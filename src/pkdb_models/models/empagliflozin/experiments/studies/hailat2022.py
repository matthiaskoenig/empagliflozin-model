"""Simulation experiment of the study Hailat2022."""

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


class Hailat2022(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Hailat2022.

    Two bioequivalence studies of 10 mg empagliflozin (test and reference) in
    healthy Jordanians: fasting (78.45 kg, the tablet with 240 ml of a 20 %
    glucose solution) and fed (63.26 kg, a high-fat meal 30 min before the
    dose), serum to 72 hr (Fig3, labelled nmol/l, read as µg/l by the
    curation, which matches Cmax 75 µg/l of Table 1).

    The study is an outlier, it is internally inconsistent and disagrees with
    every other study: tmax 4.8 hr with a second peak at 5-6 hr, a half-life of
    4.0-4.5 hr (the other studies 8-14 hr), AUC0-inf 527-539 µg*h/l (41 % of
    Chen2020 at the same dose), the fed curves are 3 times below the fasting
    curves while the text reports 50-74 % more exposure with food, and the
    units of Fig3 are unclear.
    """

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    #: the body weight of the cohorts [kg]
    bodyweights: ClassVar[dict[str, float]] = {
        "fasted": 78.45,
        "fed": 63.26,
    }

    interventions: ClassVar[list[str]] = [
        "fasting_test",
        "fasting_reference",
        "fed_test",
        "fed_reference",
    ]
    #: the names of the interventions in the legends
    legends: ClassVar[dict[str, str]] = {
        "fasting_test": "fasted, test",
        "fasting_reference": "fasted, reference",
        "fed_test": "fed, test",
        "fed_reference": "fed, reference",
    }

    symbols: ClassVar[dict[str, str]] = {
        "test": "s",  # square
        "reference": "D",  # diamond
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[label] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for condition, bodyweight in self.bodyweights.items():
            simulations[f"po_emp10_{condition}"] = Simulation(
                time_unit="hr",
                end=52,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(bodyweight, "kg"),
                    "[KI__fpg]": Q(self.fpg, "mM"),
                    "GU__fed": Q(self.fasting_map[condition], "dimensionless"),
                    "PODOSE_emp": Q(10, "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            condition = "fasted" if "fasting" in intervention else "fed"
            mappings[f"fm_hailat2022_{intervention}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=intervention,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_po_emp10_{condition}",
                    xid="time",
                    yid="[Cve_emp]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED
                    if "fasting" in intervention
                    else Fasting.FED,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        # simulations
        for condition in self.bodyweights:
            plots[0].add_data(
                task=f"task_po_emp10_{condition}",
                xid="time",
                yid="[Cve_emp]",
                label=f"10 mg Emp ({condition})",
                color=self.fasting_colors[condition],
            )

        # data
        for intervention in self.interventions:
            condition = "fasted" if "fasting" in intervention else "fed"
            formulation = intervention.split("_")[-1]  # "test" or "reference"
            plots[0].add_data(
                dataset=intervention,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"10 mg Emp ({self.legends[intervention]})",
                color=self.fasting_colors[condition],
                marker=self.symbols[formulation],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Hailat2022, output_dir=Hailat2022.__name__)
