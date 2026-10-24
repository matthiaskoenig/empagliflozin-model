"""Simulation experiment of the study Macha2014f."""

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


class Macha2014f(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2014f.

    Single dose of 50 mg empagliflozin after an overnight fast in subjects with
    normal renal function and with mild, moderate and severe renal impairment
    (eGFR > 90, 60-89, 30-59, < 30 ml/min/1.73 m2, MDRD) and with end stage
    renal disease on dialysis (ESRD). The groups with normal renal function and
    with mild and moderate impairment have T2DM, 4 of 8 subjects of the severe
    group and none of the ESRD group.
    """

    # [mM] fasting plasma glucose of the groups with T2DM (assumed, not
    # reported; the severe group, half T2DM, is simulated with it as well);
    # the ESRD group without diabetes has the healthy reference
    fpg_t2dm = 10
    fpgs: ClassVar[dict[str, float]] = {
        "normal": fpg_t2dm,
        "mild": fpg_t2dm,
        "moderate": fpg_t2dm,
        "severe": fpg_t2dm,
        "esrd": EmpagliflozinSimulationExperiment.fpg_healthy,
    }
    # [kg] median weight (Table 1)
    bodyweight: ClassVar[dict[str, float]] = {
        "normal": 81.7,
        "mild": 70.2,
        "moderate": 78.6,
        "severe": 75.1,
        "esrd": 77.25,
    }
    # [ml/min/1.73 m2] mean eGFR of the groups (individual eGFR of Figure 2);
    # the renal function is eGFR/101 (GFR of healthy kidneys at 75 kg)
    egfr: ClassVar[dict[str, float]] = {
        "normal": 108,
        "mild": 75,
        "moderate": 46,
        "severe": 23,
        "esrd": 8,
    }
    renal_functions: ClassVar[dict[str, float]] = {
        group: gfr / EmpagliflozinSimulationExperiment.gfr_healthy
        for group, gfr in egfr.items()
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
        "esrd": EmpagliflozinSimulationExperiment.renal_colors[
            "End stage renal disease"
        ],
    }
    health: ClassVar[dict[str, Health]] = {
        "normal": Health.T2DM,
        "mild": Health.T2DM_RENAL_IMPAIRMENT,
        "moderate": Health.T2DM_RENAL_IMPAIRMENT,
        # 4 of 8 subjects with T2DM
        "severe": Health.T2DM_RENAL_IMPAIRMENT,
        # dialysis, no T2DM; the timing of the dialysis is not reported and the
        # dialysis is not modelled
        "esrd": Health.RENAL_IMPAIRMENT,
    }
    groups: ClassVar[list[str]] = list(renal_functions.keys())
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
                    "BW": Q(self.bodyweight[group], "kg"),
                    "[KI__fpg]": Q(self.fpgs[group], "mM"),
                    "KI__f_renal_function": Q(
                        self.renal_functions[group], "dimensionless"
                    ),  # [0, 1]  <=> [0, 100] gfr
                    "PODOSE_emp": Q(50, "mg"),
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
                        health=self.health[group],
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
            num_cols=3,
            name=f"{self.__class__.__name__} (T2DM, Renal Impairment, ESRD)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)
        plots[2].set_yaxis(self.label_uge, unit=self.unit_uge)

        for kp, (sid, name) in enumerate(self.info.items()):
            for group in self.groups:
                # simulation
                plots[kp].add_data(
                    task=f"task_{group}",
                    xid="time",
                    yid=sid,
                    label=f"50 mg Emp ({group})",
                    color=self.colors[group],
                )
                # data
                plots[kp].add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"50 mg Emp ({group})",
                    color=self.colors[group],
                    linestyle="" if "urine" in name else "--",
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Macha2014f, output_dir=Macha2014f.__name__)
