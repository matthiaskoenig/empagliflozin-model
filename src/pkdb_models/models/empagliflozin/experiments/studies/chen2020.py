"""Simulation experiment of the study Chen2020."""

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


class Chen2020(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Chen2020.

    Two cohorts of healthy Chinese subjects (fasting 62.93 kg, fed 60.99 kg)
    take a single tablet of 10 mg empagliflozin, test (T, generic) and
    reference (R, Jardiance) in a crossover, after a fast of at least 10 hr or
    30 min after a high-fat meal (Fig2). T and R are the same subjects and
    bioequivalent. The SDs of Fig2 are the upper error bars from 4 hr on; the
    bars of the absorption phase overlap, there the SD is the CV of Cmax of
    Table 3 times the mean.
    """

    # groups are "<condition>_<formulation>", test (T) and reference (R)
    bodyweights: ClassVar[dict[str, float]] = {
        "Fasting_T": 62.93,
        "Fasting_R": 62.93,
        "Fed_T": 60.99,
        "Fed_R": 60.99,
    }  # [kg]
    groups: ClassVar[list[str]] = list(bodyweights)

    fasting_key: ClassVar[dict[str, str]] = {
        "Fasting": "fasted",
        "Fed": "fed",
    }
    # panel of the figure by condition
    panels: ClassVar[dict[str, int]] = {"Fasting": 0, "Fed": 1}

    symbols: ClassVar[dict[str, str]] = {
        "T": "s",  # square
        "R": "D",  # diamond
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if label.endswith("EMP10"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[label] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        for group, bodyweight in self.bodyweights.items():
            condition, _ = group.split("_")  # "Fasting" or "Fed"
            simulations[f"po_emp10_{group}"] = Simulation(
                time_unit="hr",
                end=73,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(bodyweight, "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "GU__fed": Q(
                        self.fasting_map[self.fasting_key[condition]], "dimensionless"
                    ),
                    "PODOSE_emp": Q(10, "mg"),
                },
            )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            mappings[f"fm_po_emp10_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{group}_EMP10",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_po_emp10_{group}",
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
                    if group.startswith("Fasting")
                    else Fasting.FED,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-0.5, max=25),
            legend=True,
        )
        for condition, kp in self.panels.items():
            plots[kp].set_title(condition)
            plots[kp].set_yaxis(
                self.label_emp_plasma, unit=self.unit_emp, min=-0.02, max=0.55
            )

        for group in self.groups:
            condition, formulation = group.split("_")  # "Fasting"/"Fed", "T"/"R"
            color = self.fasting_colors[self.fasting_key[condition]]
            marker = self.symbols[formulation]
            kp = self.panels[condition]

            # the simulations of T and R are the same
            if formulation == "R":
                plots[kp].add_data(
                    task=f"task_po_emp10_{group}",
                    xid="time",
                    yid="[Cve_emp]",
                    label="10 mg Emp",
                    color=color,
                )
            plots[kp].add_data(
                dataset=f"{group}_EMP10",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"10 mg Emp ({formulation})",
                color=color,
                marker=marker,
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Chen2020, output_dir=Chen2020.__name__)
