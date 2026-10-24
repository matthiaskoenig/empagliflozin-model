"""Simulation experiment of the study Friedrich2013."""

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


class Friedrich2013(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Friedrich2013.

    Sixteen healthy men (81.5 kg) take 50 mg empagliflozin once daily for 5 days
    (EMP50_D5), directly followed by 7 days of 50 mg empagliflozin with 5 mg
    linagliptin (EMP50_D7, the twelfth dose of empagliflozin), with water after
    an overnight fast of at least 10 hr. Plasma (Fig2) and the cumulative urinary
    excretion (Fig3, collections 0-2, 2-4, 4-8, 8-12 and 12-24 hr) after the last
    dose, without error bars.

    Linagliptin does not change the exposure of empagliflozin: AUCtau,ss ratio
    1.02 (0.97-1.07), Cmax,ss 0.88 (0.79-0.99), fe0-24 20.9 vs 20.6 %. The arm
    is simulated without an interaction factor and is validation data.
    """

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    bodyweight = 81.5  # [kg]

    interventions: ClassVar[list[str]] = [
        "EMP50_D5",
        "EMP50_D7",
    ]
    colors: ClassVar[dict[str, str]] = {
        "EMP50_D5": "tab:blue",
        "EMP50_D7": "tab:orange",
    }
    #: the number of doses of empagliflozin before the profile
    n_doses: ClassVar[dict[str, int]] = {
        "EMP50_D5": 5,
        "EMP50_D7": 12,
    }
    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole
                if label.startswith("empagliflozin_urine_"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        dose = Q(50, "mg")
        for intervention in self.interventions:
            # a dose of empagliflozin every 24 hr, 5 doses alone (D5) or 5
            # doses alone followed by 7 doses with linagliptin (D7, 12 doses);
            # the time is relative to the last dose; the urine is collected
            # after the last dose (the amount is reset at every dose)
            n_doses = self.n_doses[intervention]
            end = 24 * (n_doses - 1) + 25  # [hr]
            simulations[f"po_{intervention}"] = Simulation(
                time_unit="hr",
                end=end,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "[KI__fpg]": Q(self.fpg, "mM"),
                    "PODOSE_emp": dose,
                },
                changes=[
                    Change(
                        [24 * k for k in range(1, n_doses)],
                        {"Aurine_emp": Q(0, "mmole")},
                    ),
                    self.add_doses(
                        [24 * k for k in range(1, n_doses)], {"PODOSE_emp": dose}
                    ),
                ],
                time_shift=-24 * (n_doses - 1),
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.LINAGLIPTIN
                        if "D7" in intervention
                        else Coadministration.NONE,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2_3",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )

        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-5, max=25),
            legend=True,
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)

        for kp, (sid, name) in enumerate(self.info.items()):
            # simulation
            plots[kp].add_data(
                task="task_po_EMP50_D5",
                xid="time",
                yid=sid,
                label="50 mg Emp",
                color="black",
            )

            # data
            for intervention in self.interventions:
                label = (
                    "50 mg Emp" if intervention.endswith("D5") else "50 mg Emp + LIN"
                )

                plots[kp].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Friedrich2013, output_dir=Friedrich2013.__name__)
