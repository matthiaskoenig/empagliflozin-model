"""Simulation experiment of the study Macha2014b."""

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


class Macha2014b(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2014b.

    Single dose of 50 mg empagliflozin in subjects with normal hepatic function
    and with mild, moderate and severe hepatic impairment.
    """

    bodyweights: ClassVar[dict[str, float]] = {
        "normal": 75,
        "mild": 82.5,
        "moderate": 75.5,
        "severe": 70.5,
    }  # [kg]
    colors: ClassVar[dict[str, str]] = {
        "normal": EmpagliflozinSimulationExperiment.cirrhosis_colors["Control"],
        "mild": EmpagliflozinSimulationExperiment.cirrhosis_colors["Mild cirrhosis"],
        "moderate": EmpagliflozinSimulationExperiment.cirrhosis_colors[
            "Moderate cirrhosis"
        ],
        "severe": EmpagliflozinSimulationExperiment.cirrhosis_colors[
            "Severe cirrhosis"
        ],
    }
    f_cirrhosis: ClassVar[dict[str, float]] = {
        "normal": EmpagliflozinSimulationExperiment.cirrhosis_map["Control"],
        "mild": EmpagliflozinSimulationExperiment.cirrhosis_map["Mild cirrhosis"],
        "moderate": EmpagliflozinSimulationExperiment.cirrhosis_map[
            "Moderate cirrhosis"
        ],
        "severe": EmpagliflozinSimulationExperiment.cirrhosis_map["Severe cirrhosis"],
    }
    groups: ClassVar[list[str]] = list(bodyweights.keys())
    # observable and data of the fit mappings
    info: ClassVar[list[tuple[str, str]]] = [
        ("[Cve_emp]", "empagliflozin"),
        ("Aurine_emp", "empagliflozin_urine"),
        ("KI__UGE", "uge"),
        ("Aurine_eg", "eg_urine"),
        ("Aurine_emptot", "emptot_urine"),
    ]
    panels: ClassVar[dict[str, int]] = {
        "empagliflozin": 0,
        "empagliflozin_urine": 1,
        "eg_urine": 2,
        "emptot_urine": 3,
        "uge": 4,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig3", "Fig4", "Tab2A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for group in self.groups:
            simulations[self.to_sid(group)] = Simulation(
                time_unit="hr",
                end=100,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[group], "kg"),
                    # 1 subject with T2DM in every hepatic impairment group
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "PODOSE_emp": Q(50, "mg"),
                    "f_cirrhosis": Q(self.f_cirrhosis[group], "dimensionless"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info:
            for group in self.groups:
                mappings[self.to_sid(f"fm_{name}_{group}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd=None if name == "empagliflozin" else "mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{group}"),
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
                        health=Health.HEALTHY
                        if group == "normal"
                        else Health.HEPATIC_IMPAIRMENT,
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
            num_cols=5,
            name=f"{self.__class__.__name__} (Hepatic impairment)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)
        plots[2].set_yaxis(self.label_eg_urine, unit=self.unit_eg_urine)
        plots[3].set_yaxis(self.label_emptot_urine, unit=self.unit_emptot_urine)
        plots[4].set_yaxis(self.label_uge, unit=self.unit_uge)

        for sid, name in self.info:
            kp = self.panels[name]
            for group in self.groups:
                # simulation
                plots[kp].add_data(
                    task=self.to_sid(f"task_{group}"),
                    xid="time",
                    yid=sid,
                    label=f"50 mg Emp ({group})",
                    color=self.colors[group],
                )

                # data
                plots[kp].add_data(
                    dataset=self.to_sid(f"{name}_{group}"),
                    xid="time",
                    yid="mean",
                    yid_sd=None if name == "empagliflozin" else "mean_sd",
                    count="count",
                    label=f"50 mg Emp ({group})",
                    color=self.colors[group],
                    linestyle="--" if kp not in {2, 3} else "",
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Macha2014b, output_dir=Macha2014b.__name__)
