"""Simulation experiment of the study ElDash2021."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    DataType,
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


class ElDash2021(EmpagliflozinSimulationExperiment):
    """Simulation experiment of ElDash2021.

    Sixteen healthy Egyptian men (74.5 kg, fasting glucose 84.4 mg/dl) take a
    tablet of 10 mg empagliflozin in the evening (7 pm, phase 1) and after 7
    days in the morning (9 am, phase 2), after a fast of 12 hr, each dose
    directly followed by a meal (dinner, breakfast): the meal falls into the
    absorption (tmax 2.5 and 1.6 hr), both arms are fed. The model has no time
    of day, one simulation describes both arms. Plasma (Fig4, means and the
    individual data) and the urinary glucose excretion (Table 2, the individual
    amounts of the data of the study).

    The UGE of the curation (UGE0-24 15.4 g evening, 12.4 g morning, from the
    individual data) is 5-9 times lower than Table 2 of the paper (amounts in
    "g/dl", 69 and 39 g in 12-24 hr, called UGE0-24 in the text, and the
    "cumulative" values do not rise): the UGE is not usable.
    """

    interventions: ClassVar[list[str]] = ["morning", "evening"]
    colors: ClassVar[dict[str, str]] = {
        "morning": "black",
        "evening": "tab:blue",
    }
    subjects: ClassVar[list[str]] = [f"S{k + 1}" for k in range(16)]
    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "KI__UGE": "uge",
    }

    fpg = 84.4 / 18  # [mM] (healthy subjects)
    bodyweight = 74.5  # [kg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig4", "Fig4A", "Tab2", "Tab2A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion, mean data (Fig4, Tab2) and individual data
                # (Fig4A, Tab2A)
                if label.startswith("empagliflozin_"):
                    if fig_id in ["Fig4", "Tab2"]:
                        dset.unit_conversion("mean", 1 / self.Mr.emp)
                    elif fig_id in ["Fig4A", "Tab2A"]:
                        dset.unit_conversion("value", 1 / self.Mr.emp)

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        simulations["po_EMP10"] = Simulation(
            time_unit="hr",
            end=51,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg, "mM"),
                "GU__fed": Q(self.fasting_map["fed"], "dimensionless"),
                "PODOSE_emp": Q(10, "mg"),
            },
        )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            for sid, name in self.info.items():
                # mean
                mappings[f"po_EMP10_{name}_{intervention}"] = FitMapping(
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
                        task="task_po_EMP10",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=self._metadata(name, DataType.MEAN),
                )

                # subjects
                for subject in self.subjects:
                    mappings[f"po_EMP10_{name}_{intervention}_{subject}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"{name}_{intervention}_{subject}",
                            xid="time",
                            yid="value",
                            yid_sd=None,
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task="task_po_EMP10",
                            xid="time",
                            yid=sid,
                        ),
                        metadata=self._metadata(name, DataType.INDIVIDUAL),
                    )

        return mappings

    @staticmethod
    def _metadata(name: str, data_type: DataType) -> PKPDMappingMetaData:
        """Metadata of the mean or the individual data of an observable."""
        return PKPDMappingMetaData(
            tissue=Tissue.URINE if "uge" in name else Tissue.PLASMA,
            route=Route.PO,
            application_form=ApplicationForm.TABLET,
            dosing=Dosing.SINGLE,
            health=Health.HEALTHY,
            fasting=Fasting.FED,
            coadministration=Coadministration.NONE,
            data_type=data_type,
        )

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig4_Tab2",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[1].set_xaxis(self.label_time, unit=self.unit_time, min=-0.05, max=25)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_uge, unit=self.unit_uge)

        # simulation
        for k, sid in enumerate(self.info):
            plots[k].add_data(
                task="task_po_EMP10",
                xid="time",
                yid=sid,
                label="10 mg Emp",
                color="black",
                alpha=0.7,
            )

        # individual data (morning vs evening)
        for intervention in self.interventions:
            for k, name in enumerate(self.info.values()):
                for subject in self.subjects:
                    plots[k].add_data(
                        dataset=f"{name}_{intervention}_{subject}",
                        xid="time",
                        yid="value",
                        yid_sd=None,
                        count="count",
                        label=None,
                        color=self.colors[intervention],
                        alpha=0.5,
                        marker="",
                    )

        # mean data (morning vs evening)
        for intervention in self.interventions:
            for k, name in enumerate(self.info.values()):
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"10 mg Emp ({intervention})",
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(ElDash2021, output_dir=ElDash2021.__name__)
