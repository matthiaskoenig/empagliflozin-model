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

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects, assuming 5 mM)
    bodyweight = 79  # [kg]
    gfr = EmpagliflozinSimulationExperiment.gfr_healthy  # [ml/min] (healthy subjects, assuming 100 ml/min)

    doses = {
        "placebo": 0,
        "EMP0.5": 0.5,
        "EMP2.5": 2.5,
        "EMP10": 10,
        "EMP25": 25,
        "EMP50": 50,
        # "OGTT50": 50,  not using this subset
        "EMP100": 100,
        "EMP200": 200,
        "EMP400": 400,
        "EMP800": 800,
    }
    interventions = list(doses.keys())

    # colors = {
    #     "placebo": "black",
    #     "EMP0.5": "#1f77b4",
    #     "EMP2.5": "#ff7f0e",
    #     "EMP10": "#2ca02c",
    #     "EMP25": "#d62728",
    #     "EMP50": "#9467bd",
    #     # "OGTT50": "tab:green",  not using this subset
    #     "EMP100": "#8c564b",
    #     "EMP200": "#e377c2",
    #     "EMP400": "#7f7f7f",
    #     "EMP800": "#bcbd22"
    # }

    # [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

    info = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
        "KI__UGE": "uge",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
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
                        # "KI__f_renal_function": Q_(self.gfr / 100, "dimensionless"),  # [0, 1]  <=> [0, 100] gfr
                        "PODOSE_emp": Q_(dose, "mg"),
                    },
                )]
            )

        # console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for kp, sid in enumerate(self.info):
            name = self.info[sid]
            for intervention in self.interventions:

                if name == "empagliflozin" and intervention == "placebo":
                    continue

                if name == "empagliflozin_urine":
                    # FIXME: add missing urinary emp data
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
        return mappings

    def figures(self) -> Dict[str, Figure]:
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
        for k in range(3):
            plots[k].xaxis.max = 30

        for kp, sid in enumerate(self.info):
            name = self.info[sid]

            for intervention in self.interventions:
                dose = self.doses[intervention]
                color = self.dose_colors[dose]
                dose_label = "Placebo" if dose == 0 else f"{dose} mg Emp"

                # simulation
                plots[kp].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=dose_label,
                    color=color,
                )

                if name == "empagliflozin" and intervention == "placebo":
                    continue

                if name == "empagliflozin_urine":
                    # FIXME: add missing urinary emp data
                    continue

                # data
                plots[kp].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=dose_label,
                    color=color,
                )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Seman2013, output_dir=Seman2013.__name__)
