from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)
from pkdb_models.models.empagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, EmpagliflozinMappingMetaData, Coadministration

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.helpers import run_experiments


class Sarashina2013(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Sarashina2013."""
    doses = {
        "placebo": 0,
        "EMP1": 1,
        "EMP5": 5,
        "EMP10": 10,
        "EMP25": 25,
        "EMP100": 100,
        # "OGTT10": 10, not using this subset
    }
    bodyweights = {
        "placebo": 55.7,
        "EMP1": 60.4,
        "EMP5": 66,
        "EMP10": 63.1,
        "EMP25": 66.4,
        "EMP100": 66,
        # "OGTT10": 56.8,
    }  # [kg]
    interventions = list(doses.keys())
    colors = {
        "placebo": "black",
        "EMP1": "tab:blue",
        "EMP5": "tab:green",
        "EMP10": "tab:red",
        "EMP25": "tab:purple",
        "EMP100": "darkgrey",
        # "OGTT10": "grey"
    }
    info = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
        "KI__UGE": "uge",
    }

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects, assuming 5 mM)
    bodyweight = 78  # [kg]
    gfr = EmpagliflozinSimulationExperiment.gfr_healthy  # [ml/min] (healthy subjects, assuming 100 ml/min)

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Tab3A"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion
                if label.startswith("empagliflozin_urine"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)

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

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}
        for kp, sid in enumerate(self.info):
            name = self.info[sid]
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

        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_2_Tab3A",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        Figure.legend_fontsize = 10
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)
        plots[2].set_yaxis(self.label_uge, unit=self.unit_uge)

        for kp, sid in enumerate(self.info):
            name = self.info[sid]

            for intervention in self.interventions:
                dose = self.doses[intervention]
                color = self.dose_colors[dose]
                dose_label = "Placebo" if dose == 0 else f"{dose} mg Emp"

                # simulation (no legend entry)
                plots[kp].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=dose_label,
                    color=color,
                )

                if name.startswith("empagliflozin") and intervention == "placebo":
                    continue

                # data
                plots[kp].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=dose_label,
                    color=color,
                    linestyle="" if "urine" in name else "--",
                )

        return {
            fig.sid: fig,
        }



if __name__ == "__main__":
    run_experiments(Sarashina2013, output_dir=Sarashina2013.__name__)
