"""Simulation experiment of the study Sarashina2013."""

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


class Sarashina2013(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Sarashina2013.

    Single rising doses of empagliflozin (1 mg to 100 mg) and placebo in healthy
    Japanese men (n=6 per dose, placebo n=12), tablets after an overnight fast
    of >= 10 hr. Body weights are the medians of the groups (Table 1). The
    group with an oral glucose tolerance test (10 mg, 75 g glucose 1 hr after
    the dose) is not simulated (the model has no glucose load). The mean
    plasma glucose of the day (AUEC0-24 2480 mg*h/dl, 5.7 mM) is the healthy
    fasting glucose of the model with its mean difference to the fasting
    glucose (5.0 + 0.7 mM).
    """

    doses: ClassVar[dict[str, float]] = {
        "placebo": 0,
        "EMP1": 1,
        "EMP5": 5,
        "EMP10": 10,
        "EMP25": 25,
        "EMP100": 100,
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "placebo": 55.7,
        "EMP1": 60.4,
        "EMP5": 66,
        "EMP10": 63.1,
        "EMP25": 66.4,
        "EMP100": 66,
    }  # [kg]
    interventions: ClassVar[list[str]] = list(doses.keys())
    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
        "KI__UGE": "uge",
    }

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Tab3A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion
                if label.startswith("empagliflozin_urine"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)

                dsets[label] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention, dose in self.doses.items():
            simulations[intervention] = Simulation(
                time_unit="hr",
                end=75,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[intervention], "kg"),
                    "[KI__fpg]": Q(self.fpg, "mM"),
                    "PODOSE_emp": Q(dose, "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                if name.startswith("empagliflozin") and intervention == "placebo":
                    continue

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
                        task=f"task_{intervention}",
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
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_2_Tab3A",
            num_rows=1,
            num_cols=4,
            name=f"{self.__class__.__name__} (Healthy)",
            # one legend per row, right of the grid (the panels show the same
            # doses)
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        # Fig1B
        plots[1].set_yaxis(
            self.label_emp_plasma, unit=self.unit_emp, scale="log", min=1e-3
        )
        plots[2].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)
        plots[3].set_yaxis(self.label_uge, unit=self.unit_uge)

        panels = [
            ("[Cve_emp]", "empagliflozin", plots[0]),
            ("[Cve_emp]", "empagliflozin", plots[1]),
            ("Aurine_emp", "empagliflozin_urine", plots[2]),
            ("KI__UGE", "uge", plots[3]),
        ]
        for sid, name, plot in panels:
            for intervention, dose in self.doses.items():
                color = self.dose_colors[dose]
                dose_label = "Placebo" if dose == 0 else f"{dose} mg Emp"

                # simulation
                plot.add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=dose_label,
                    color=color,
                )

                # no empagliflozin data for placebo
                if name.startswith("empagliflozin") and intervention == "placebo":
                    continue

                # data (the legend entries of the plasma and the glucose, the
                # urine has the same markers)
                plot.add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=None if "urine" in name else dose_label,
                    color=color,
                    linestyle="" if "urine" in name else "--",
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Sarashina2013, output_dir=Sarashina2013.__name__)
