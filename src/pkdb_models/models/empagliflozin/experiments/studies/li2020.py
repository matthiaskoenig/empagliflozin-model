"""Simulation experiment of the study Li2020."""

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


class Li2020(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Li2020.

    Bioequivalence of a generic (test) and the brand-named (reference) tablet of
    25 mg empagliflozin in healthy Chinese subjects, in a fasting and a fed arm
    of 24 subjects each (crossover of test and reference in each arm).

    Fasting arm: overnight fast, the dose with 240 ml water containing 20 %
    glucose and 60 ml of 20 % glucose every hour until 4 hr (about 100 g oral
    glucose). Fed arm: a high-fat, high-calorie breakfast 30 min before the
    dose. Both arms have the same tmax of 4 hr (2-6 hr) and the same Cmax and
    AUC, i.e. the glucose drink delays the absorption as the meal does: the
    fasting arm is simulated with the absorption of the fed state.

    Only the reference is mapped: test and reference are the same subjects and
    bioequivalent, and the test curves are hidden behind the reference curves
    of the figures (not digitizable).
    """

    bodyweights: ClassVar[dict[str, float]] = {
        "fasting": 58.2,
        "fed": 59.9,
    }  # [kg]
    conditions: ClassVar[list[str]] = ["fasting", "fed"]

    #: the fasting state of the study (metadata) by arm
    fasting_key: ClassVar[dict[str, str]] = {
        "fasting": "fasted",
        "fed": "fed",
    }
    #: the state of the intestine in the simulation: the oral glucose load of the
    #: fasting arm delays the absorption like a meal (tmax 4 hr in both arms)
    absorption_key: ClassVar[dict[str, str]] = {
        "fasting": "fed",
        "fed": "fed",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("empagliflozin"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[label] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for condition in self.conditions:
            simulations[f"po_emp25_{condition}"] = Simulation(
                time_unit="hr",
                end=73,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[condition], "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "KI__f_renal_function": Q(1.0, "dimensionless"),
                    "GU__fed": Q(
                        self.fasting_map[self.absorption_key[condition]],
                        "dimensionless",
                    ),
                    "PODOSE_emp": Q(25, "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for condition in self.conditions:
            mappings[f"fm_po_emp25_ref_{condition}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"empagliflozin_ref_{condition}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_po_emp25_{condition}",
                    xid="time",
                    yid="[Cve_emp]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    # the fasting arm fasted overnight, but drank about 100 g
                    # glucose with and after the dose
                    fasting=Fasting.FASTED if condition == "fasting" else Fasting.FED,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2_3",
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_plasma, unit=self.unit_emp, scale="log")

        for plot in plots:
            for condition in self.conditions:
                color = self.fasting_colors[self.fasting_key[condition]]
                label = {
                    "fasting": "25 mg Emp (fasting, glucose drink)",
                    "fed": "25 mg Emp (fed)",
                }[condition]
                plot.add_data(
                    task=f"task_po_emp25_{condition}",
                    xid="time",
                    yid="[Cve_emp]",
                    label=label,
                    color=color,
                )
                plot.add_data(
                    dataset=f"empagliflozin_ref_{condition}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=color,
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Li2020, output_dir=Li2020.__name__)
