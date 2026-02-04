from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)
from pkdb_models.models.empagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, EmpagliflozinMappingMetaData

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.helpers import run_experiments


class Chen2015a(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Chen2015a."""

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects, assuming 5 mM)
    bodyweight = 79  # [kg]
    gfr = EmpagliflozinSimulationExperiment.gfr_healthy  # [ml/min] (healthy subjects, assuming 100 ml/min)

    info = {
        # Fig2
        "[Cve_emptot]": "EMPTOT_plasma",
        "[Cve_emp]": "EMP_plasma",
        "[Cve_eg]": "EG_plasma",
        # Fig1
        "Afeces_emptot": "EMPTOT_feces",
        "Afeces_emp": "EMP_feces",
        "Afeces_eg": "EG_feces",
        "Aurine_emptot": "EMPTOT_urine",
        "Aurine_emp": "EMP_urine",
        "Aurine_eg": "EG_urine",
    }
    colors = {
        "[Cve_emptot]": "black",
        "[Cve_emp]": "tab:blue",
        "[Cve_eg]": "tab:orange",
        # Fig1
        "Afeces_emptot": "black",
        "Afeces_emp": "tab:blue",
        "Afeces_eg": "tab:orange",
        "Aurine_emptot": "black",
        "Aurine_emp":  "tab:blue",
        "Aurine_eg": "tab:orange",
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

        tcsims[f"po_C14EMP"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=180 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "[KI__fpg]": Q_(self.fpg, "mM"),
                    "KI__f_renal_function": Q_(self.gfr/100, "dimensionless"),   # [0, 1]  <=> [0, 100] gfr
                    "PODOSE_emp": Q_(50, "mg"),  # radioactive empagliflozin dose and nonradioactive (49.1 + 0.9)
                },
            )]
        )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}
        for sid, name in self.info.items():
            if name.endswith("_plasma"):
                tissue = Tissue.PLASMA
            elif name.endswith("_urine"):
                tissue = Tissue.URINE
            elif name.endswith("_feces"):
                tissue = Tissue.FECES

            mappings[f"fm_po_emp_{name}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{name}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_C14EMP", xid="time", yid=f"{sid}",
                ),
                metadata=EmpagliflozinMappingMetaData(
                    tissue=tissue,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR,
                ),
            )
            # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_rows=1,
            num_cols=4,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(label="Plasma", unit=self.unit_emp)
        plots[1].set_yaxis(self.label_eg, unit=self.unit_eg)
        plots[2].set_yaxis(label="Urine", unit=self.unit_emp_urine)
        plots[3].set_yaxis(label="Feces", unit=self.unit_emp_feces)
        plots[1].yaxis.min = -0.05
        plots[1].yaxis.max = 0.3

        for sid, name in self.info.items():
            if name.endswith("_plasma") and "EG" not in name:
                kp = 0
            if name.endswith("_plasma") and "EG" in name:
                kp = 1
            elif name.endswith("_urine"):
                kp = 2
            elif name.endswith("_feces"):
                kp = 3

            label_map = {
                "emptot": "EMPTOT",
                "_eg": "EG",
            }

            # determine species
            species = "EMP"
            for k, v in label_map.items():
                if k in sid.lower() or k in name.lower():
                    species = v
                    break

            label = f"50 mg Emp ({species})"

            # simulation
            plots[kp].add_data(
                task="task_po_C14EMP",
                xid="time",
                yid=sid,
                label=label,
                color=self.colors[sid],
            )

            # data
            plots[kp].add_data(
                dataset=name,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label,
                color=self.colors[sid],
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Chen2015a, output_dir=Chen2015a.__name__)
