"""Simulation experiment of the study Seman2013."""

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
from sbmlsim.simulation import Simulation

from pkdb_models.models.empagliflozin import MODEL
from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)


class Seman2013(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Seman2013.

    Single rising doses of empagliflozin (0.5 mg to 800 mg) and placebo in
    healthy white men after an overnight fast of >= 10 hr (part 1, median
    79 kg), and 50 mg fasted and directly after a high-fat, high-calorie
    breakfast in a crossover (part 2, Fig4, 14 men, median 82 kg). The OGTT
    group of part 1 (50 mg, 75 g glucose 1 hr after the dose) is not simulated.
    The mean plasma glucose of placebo (AUEC0-72 7210 mg*h/dl, 5.6 mM) is the
    healthy fasting glucose of the model with its mean difference to the
    fasting glucose (5.0 + 0.7 mM). The exposure rises less than proportional
    to the dose above 400 mg (400 to 800 mg: AUC x1.5, Cmax x1.3, fe 15.1 to
    11.0 %), the model is linear.
    """

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy
    #: median body weight of part 1 (single rising doses)
    bodyweight = 79  # [kg]
    #: median body weight of part 2 (food effect, Fig4)
    bodyweight_food = 82  # [kg]

    doses: ClassVar[dict[str, float]] = {
        "placebo": 0,
        "EMP0.5": 0.5,
        "EMP2.5": 2.5,
        "EMP10": 10,
        "EMP25": 25,
        "EMP50": 50,
        "EMP100": 100,
        "EMP200": 200,
        "EMP400": 400,
        "EMP800": 800,
    }
    interventions: ClassVar[list[str]] = list(doses.keys())

    fig4_conditions: ClassVar[list[str]] = ["fasted", "fed"]

    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
        "KI__UGE": "uge",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "FigS2", "Fig4"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        for intervention, dose in self.doses.items():
            simulations[self.to_sid(intervention)] = Simulation(
                time_unit="hr",
                end=75,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "[KI__fpg]": Q(self.fpg, "mM"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "PODOSE_emp": Q(dose, "mg"),
                },
            )

        # Fig4: fasted/fed
        for condition in self.fig4_conditions:
            simulations[self.to_sid(f"EMP50_{condition}")] = Simulation(
                time_unit="hr",
                end=75,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight_food, "kg"),
                    "[KI__fpg]": Q(self.fpg, "mM"),
                    "GU__fed": Q(self.fasting_map[condition], "dimensionless"),
                    "PODOSE_emp": Q(50, "mg"),
                },
            )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                # no empagliflozin data for placebo
                if (
                    name in {"empagliflozin", "empagliflozin_urine"}
                    and intervention == "placebo"
                ):
                    continue

                mappings[self.to_sid(f"fm_{name}_{intervention}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{intervention}"),
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{intervention}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA
                        if name == "empagliflozin"
                        else Tissue.URINE,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )

        # Fig4: fed and fasted fit mappings
        for condition, fasting in [("fed", Fasting.FED), ("fasted", Fasting.FASTED)]:
            mappings[self.to_sid(f"fm_empagliflozin_{condition}_EMP50")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"{condition}_EMP50"),
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_EMP50_{condition}"),
                    xid="time",
                    yid="[Cve_emp]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=fasting,
                    coadministration=Coadministration.NONE,
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        # Fig1_2
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (Healthy)",
            # one legend per row, right of the grid (the panels show the same
            # doses)
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        # semilog as Fig1A
        plots[0].set_yaxis(
            self.label_emp_plasma, unit=self.unit_emp, scale="log", min=1e-4
        )
        plots[1].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)
        plots[2].set_yaxis(self.label_uge, unit=self.unit_uge, max=100)
        # the glucose was measured to 24 hr
        plots[2].set_xaxis(self.label_time, unit=self.unit_time, min=-1, max=26)

        for kp, (sid, name) in enumerate(self.info.items()):
            for intervention in self.interventions:
                dose = self.doses[intervention]
                color = self.dose_colors[dose]
                dose_label = "Placebo" if dose == 0 else f"{dose} mg Emp"

                plots[kp].add_data(
                    task=self.to_sid(f"task_{intervention}"),
                    xid="time",
                    yid=sid,
                    label=dose_label,
                    color=color,
                )

                # no empagliflozin data for placebo
                if (
                    name in {"empagliflozin", "empagliflozin_urine"}
                    and intervention == "placebo"
                ):
                    continue

                # data (the legend entries of the plasma and the placebo
                # glucose, the other panels have the same markers)
                plots[kp].add_data(
                    dataset=self.to_sid(f"{name}_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=dose_label if kp == 0 or dose == 0 else None,
                    color=color,
                )

        # Fig4: fed vs fasted
        fig4 = Figure(
            experiment=self,
            sid="Fig4",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots4 = fig4.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots4[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots4[0].set_xaxis(self.label_time, unit=self.unit_time, max=30)

        for condition in self.fig4_conditions:
            color = self.fasting_colors[condition]
            plots4[0].add_data(
                task=self.to_sid(f"task_EMP50_{condition}"),
                xid="time",
                yid="[Cve_emp]",
                label=f"50 mg Emp ({condition})",
                color=color,
            )
            plots4[0].add_data(
                dataset=self.to_sid(f"{condition}_EMP50"),
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"50 mg Emp ({condition})",
                color=color,
            )

        return self.figures_by_sid(fig, fig4)


if __name__ == "__main__":
    MODEL.run_experiments(Seman2013, output_dir=Seman2013.__name__)
