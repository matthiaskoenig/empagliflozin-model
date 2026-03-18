from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.experiments.base_experiment import EmpagliflozinSimulationExperiment
from pkdb_models.models.empagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, EmpagliflozinMappingMetaData, Coadministration
from pkdb_models.models.empagliflozin.helpers import run_experiments


class Seman2013(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Seman2013."""

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy
    bodyweight = 79  # [kg]
    gfr = EmpagliflozinSimulationExperiment.gfr_healthy

    doses = {
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
    interventions = list(doses.keys())

    fig4_conditions = ["fasted", "fed"]

    info = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
        "KI__UGE": "uge",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "FigS2", "Fig4"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[f"{label}"] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for intervention, dose in self.doses.items():
            tcsims[intervention] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=75 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "[KI__fpg]": Q_(self.fpg, "mM"),
                        "GU__f_absorption": Q_(self.fasting_map["NR"], "dimensionless"),
                        "PODOSE_emp": Q_(dose, "mg"),
                    },
                )]
            )

        # Fig4: fasted/fed
        for condition in self.fig4_conditions:
            tcsims[f"EMP50_{condition}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=75 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "[KI__fpg]": Q_(self.fpg, "mM"),
                        "GU__f_absorption": Q_(self.fasting_map[condition], "dimensionless"),
                        "PODOSE_emp": Q_(50, "mg"),
                    },
                )]
            )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for kp, sid in enumerate(self.info):
            name = self.info[sid]
            for intervention in self.interventions:

                if name == "empagliflozin" and intervention == "placebo":
                    continue
                if name == "empagliflozin_urine" and intervention == "placebo":
                    continue

                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_{intervention}", xid="time", yid=sid,
                    ),
                    metadata=EmpagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA if name == "empagliflozin" else Tissue.URINE,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration=Coadministration.NONE,
                    ),
                )

        # Fig4: fed and fasted fit mappings
        for condition, fasting in [("fed", Fasting.FED), ("fasted", Fasting.FASTED)]:
            mappings[f"fm_empagliflozin_{condition}_EMP50"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{condition}_EMP50",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_EMP50_{condition}", xid="time", yid="[Cve_emp]",
                ),
                metadata=EmpagliflozinMappingMetaData(
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

    def figures(self) -> Dict[str, Figure]:
        # Fig1_2
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        Figure.legend_fontsize = 9
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)
        plots[2].set_yaxis(self.label_uge, unit=self.unit_uge)
        for kp in [1, 2]:
            plots[kp].xaxis.max = 25

        for kp, sid in enumerate(self.info):
            name = self.info[sid]
            for intervention in self.interventions:
                dose = self.doses[intervention]
                color = self.dose_colors[dose]
                dose_label = "Placebo" if dose == 0 else f"{dose} mg Emp"

                plots[kp].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=dose_label,
                    color=color,
                )

                if name == "empagliflozin" and intervention == "placebo":
                    continue
                if name == "empagliflozin_urine" and intervention == "placebo":
                    continue

                plots[kp].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=dose_label,
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
        plots4 = fig4.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots4[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots4[0].xaxis.max = 30

        for condition in self.fig4_conditions:
            color = self.fasting_colors[condition]
            plots4[0].add_data(
                task=f"task_EMP50_{condition}",
                xid="time",
                yid="[Cve_emp]",
                label=f"50 mg Emp ({condition})",
                color=color,
            )
            plots4[0].add_data(
                dataset=f"{condition}_EMP50",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"50 mg Emp ({condition})",
                color=color,
            )

        return {
            fig.sid: fig,
            fig4.sid: fig4,
        }


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import RESULTS_PATH_SIMULATION
    run_experiments(Seman2013, output_dir=RESULTS_PATH_SIMULATION)