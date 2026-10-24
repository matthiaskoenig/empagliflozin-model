"""Simulation experiment of the study Sarashina2014."""

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


class Sarashina2014(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Sarashina2014.

    Single dose of 25 mg empagliflozin (tablet, after an overnight fast) in
    Japanese patients with T2DM and normal renal function or mild, moderate
    and severe renal impairment (n=8 per group, Japanese eGFR >= 90, 60-<90,
    30-<60 and 15-<30 ml/min/1.73 m2). The renal function of a group is its
    mean eGFR over 101 ml/min (the eGFR of Fig2A grouped by the thresholds:
    109.1, 75.3, 46.8 and 20.8 ml/min/1.73 m2; the group means of the AUC of
    the scatter match Table II); the eGFR is not scaled to the body surface
    area of the patients. Body weights are the medians of Table I.
    """

    #: fasting plasma glucose [mM]: not reported, the mean glucose estimated
    #: from the HbA1c of Table I (ADAG, 28.7 * HbA1c - 46.7 mg/dl: 8.5, 9.1,
    #: 9.2 and 8.2 mM) minus the difference of the mean and the fasting glucose
    #: of the model (`KI__PG_delta`, 0.7 mM)
    fpgs: ClassVar[dict[str, float]] = {
        "normal": 7.8,
        "mild": 8.4,
        "moderate": 8.5,
        "severe": 7.5,
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "normal": 56.9,
        "mild": 64.65,
        "moderate": 69.2,
        "severe": 61.35,
    }  # [kg]

    #: mean eGFR of the groups [ml/min/1.73 m2] (Fig2A)
    egfrs: ClassVar[dict[str, float]] = {
        "normal": 109.1,
        "mild": 75.3,
        "moderate": 46.8,
        "severe": 20.8,
    }
    renal_functions: ClassVar[dict[str, float]] = {
        group: egfr / EmpagliflozinSimulationExperiment.gfr_healthy
        for group, egfr in egfrs.items()
    }
    colors: ClassVar[dict[str, str]] = {
        "normal": EmpagliflozinSimulationExperiment.renal_colors[
            "Normal renal function"
        ],
        "mild": EmpagliflozinSimulationExperiment.renal_colors["Mild renal impairment"],
        "moderate": EmpagliflozinSimulationExperiment.renal_colors[
            "Moderate renal impairment"
        ],
        "severe": EmpagliflozinSimulationExperiment.renal_colors[
            "Severe renal impairment"
        ],
    }
    groups: ClassVar[list[str]] = list(bodyweights.keys())

    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
        "KI__UGE": "uge",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig3", "Tab2A"]:
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
        for group in self.groups:
            simulations[group] = Simulation(
                time_unit="hr",
                end=100,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[group], "kg"),
                    "[KI__fpg]": Q(self.fpgs[group], "mM"),
                    "KI__f_renal_function": Q(
                        self.renal_functions[group], "dimensionless"
                    ),
                    "PODOSE_emp": Q(25, "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for group in self.groups:
                mappings[f"fm_{name}_{group}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{group}",
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
                        health=Health.T2DM
                        if group == "normal"
                        else Health.T2DM_RENAL_IMPAIRMENT,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_3_Tab2A",
            num_rows=1,
            num_cols=4,
            name=f"{self.__class__.__name__} (T2DM, renal impairment)",
            # one legend per row, right of the grid (the panels show the same
            # groups)
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        # the terminal phase (the semilog insert of Fig1A)
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
            for group in self.groups:
                label = f"25 mg Emp ({group})"
                # simulation
                plot.add_data(
                    task=f"task_{group}",
                    xid="time",
                    yid=sid,
                    label=label,
                    color=self.colors[group],
                )
                # data (the legend entries of the plasma and the glucose, the
                # urine has the same markers)
                plot.add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=None if "urine" in name else label,
                    color=self.colors[group],
                    linestyle="" if "urine" in name else "--",
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Sarashina2014, output_dir=Sarashina2014.__name__)
