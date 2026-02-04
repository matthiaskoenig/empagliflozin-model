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


class Zhao2015(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Zhao2015."""
    bodyweights = {
        "MULTI10": 74.0,  # [kg]
        "MULTI25": 67.0,  # [kg]
    }
    fpgs = {
        "MULTI10": 166.67/180,  # [mM]
        "MULTI25": 160.56/180,  # [mM]
    }
    doses = {
        "MULTI10": 10,
        "MULTI25": 25,
    }
    # colors = {
    #     "MULTI10": "tab:blue",
    #     "MULTI25": "tab:orange",
    # }
    gfr = EmpagliflozinSimulationExperiment.gfr_healthy  # [ml/min] (healthy subjects, assuming 100 ml/min)
    interventions = list(bodyweights.keys())
    info = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab2A"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for intervention, dose in self.doses.items():
            tc0 = Timecourse(
                start=0,
                end=48 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweights[intervention], "kg"),
                    "[KI__fpg]": Q_(self.fpgs[intervention], "mM"),
                    # "KI__f_renal_function": Q_(self.gfr/100, "dimensionless"),
                    "PODOSE_emp": Q_(dose, "mg"),
                    "Aurine_emp": Q_(0, "mmole"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_emp": Q_(dose, "mg"),
                    "Aurine_emp": Q_(0, "mmole"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=75 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_emp": Q_(dose, "mg"),
                    "Aurine_emp": Q_(0, "mmole"),
                },
            )
            tcsims[f"{intervention}"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(6)] + [tc2],
                # time_offset=-8*24*60,
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for intervention in self.interventions:
            for k, sid in enumerate(self.info):
                name = self.info[sid]
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
                        tissue=Tissue.URINE if "urine" in intervention else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.T2DM,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (T2DM)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)

        for intervention in self.interventions:
            dose = self.doses[intervention]
            color = self.dose_colors[dose]
            dose_label = f"{dose} mg Emp"

            for k, sid in enumerate(self.info):
                name = self.info[sid]

                # simulation
                plots[k].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=dose_label,
                    color=color,
                )

                # data
                plots[k].add_data(
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
    run_experiments(Zhao2015, output_dir=Zhao2015.__name__)
