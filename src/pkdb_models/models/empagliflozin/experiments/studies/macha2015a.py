"""Simulation experiment of the study Macha2015a."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.empagliflozin import MODEL
from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)


class Macha2015a(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2015a.

    Twice daily (5 mg BID) and once daily (10 mg QD) regimens of empagliflozin
    for 5 days in healthy subjects (crossover, n=16), the profiles and the
    urinary glucose excretion of day 5 (steady state). The morning doses were
    taken fasted (an overnight fast of at least 10 hr on days 1 and 5, on day
    5 2 hr before a light breakfast), the evening doses of the BID regimen 30
    min after dinner. The times are relative to the morning dose of day 5.
    """

    doses: ClassVar[dict[str, int]] = {"EMP5": 5, "EMP10": 10}
    #: the dosing interval of the regimens [hr]
    intervals: ClassVar[dict[str, int]] = {"EMP5": 12, "EMP10": 24}
    interventions: ClassVar[list[str]] = list(doses.keys())

    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "KI__UGE": "uge",
    }

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    # [kg] not reported (BMI 23.3 kg/m2, 7 men, 9 women): the reference weight
    bodyweight = EmpagliflozinSimulationExperiment.bodyweight_healthy

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if label.startswith("uge"):
                    # the urine of day 5 is collected from the morning dose:
                    # the value at 0 hr is the empty collection, not a
                    # measurement (the model resets the urinary glucose at the
                    # dose, the value before the reset is the urine of day 4)
                    df_label = df_label[df_label.time > 0]
                dsets[label] = DataSet.from_df(df_label, self.ureg)
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention, dose in self.doses.items():
            mg = Q(dose, "mg")
            interval = self.intervals[intervention]
            # the doses of days 1-5 at 0 hr (morning, fasted) and for the BID
            # regimen at 12 hr (evening, after dinner) of every day; the urinary
            # glucose is collected from the morning dose of day 5 (96 hr)
            morning = [24 * k for k in range(1, 4)]
            evening = [24 * k + 12 for k in range(5)] if interval == 12 else []
            simulations[f"po_{intervention}"] = Simulation(
                time_unit="hr",
                end=96 + 25,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "[KI__fpg]": Q(self.fpg, "mM"),
                    "KI__f_renal_function": Q(1.0, "dimensionless"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "PODOSE_emp": mg,
                },
                changes=[
                    Change(
                        morning,
                        {"GU__fed": Q(self.fasting_map["fasted"], "dimensionless")},
                    ),
                    Change(
                        [96],
                        {
                            "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                            "KI__glc_urine": Q(0, "mmole"),
                        },
                    ),
                    self.add_doses([*morning, 96, *evening], {"PODOSE_emp": mg}),
                    *(
                        [
                            Change(
                                evening,
                                {
                                    "GU__fed": Q(
                                        self.fasting_map["fed"], "dimensionless"
                                    ),
                                },
                            )
                        ]
                        if evening
                        else []
                    ),
                ],
                time_shift=-96,
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            for sid, name in self.info.items():
                mappings[f"fm_po_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if name == "uge" else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "uge" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        # the morning doses fasted; the evening doses of the
                        # BID regimen after dinner are simulated in the fed
                        # state
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-24, max=25),
            legend=True,
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_uge, unit=self.unit_uge, min=-2, max=80)
        # the urinary glucose of day 5
        plots[1].set_xaxis(self.label_time, unit=self.unit_time, min=-1, max=25)

        for intervention in self.interventions:
            dose = self.doses[intervention]
            label = f"{dose} mg Emp"
            color = self.dose_colors[dose]

            for k, (sid, name) in enumerate(self.info.items()):
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=label,
                    color=color,
                )

                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd" if name == "uge" else None,
                    count="count",
                    label=label,
                    color=color,
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Macha2015a, output_dir=Macha2015a.__name__)
