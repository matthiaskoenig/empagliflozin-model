from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from sbmlutils.console import console

from pkdb_models.models.empagliflozin.experiments.base_experiment import EmpagliflozinSimulationExperiment
from pkdb_models.models.empagliflozin.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, Coadministration, EmpagliflozinMappingMetaData
)
from pkdb_models.models.empagliflozin.helpers import run_experiments


class Macha2015b(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2015b."""

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    # fpg = 7  # [mM] (healthy subjects, application with food)
    interventions = ["emp50", "emp50, piog45"]

    colors = {
        "emp50": "black",
        "emp50, piog45": "tab:blue",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[label] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        dose = 50  # [mg]

        # Single dosing
        for intervention in self.interventions:
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "[KI__fpg]": Q_(self.fpg_healthy, "mM"),  # healthy reference value
                    "KI__f_renal_function": Q_(1.0, "dimensionless"), # healthy
                    "PODOSE_emp": Q_(dose, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_emp": Q_(dose, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=30 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_emp": Q_(dose, "mg"),
                },
            )

            tcsims[f"po_emp50_{intervention}"] = TimecourseSim(
               [tc0] + [tc1 for _ in range(3)] + [tc2],
               time_offset=-4*24*60
            )
            # console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
            mappings = {}
            for intervention in self.interventions:
                mappings[f"fm_po_emp50_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"empagliflozin_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_emp50_{intervention}", xid="time", yid=f"[Cve_emp]",
                    ),
                    metadata=EmpagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,  # overnight fast
                        coadministration=Coadministration.NONE if intervention == "emp50" else Coadministration.PIOGLITAZONE,
                    ),
                )
            # console.print(mappings)
            return mappings

    def figures(self) -> Dict[str, Figure]:
            fig = Figure(
                experiment=self,
                sid="Fig1",
                name=f"{self.__class__.__name__}",
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time, min=-24), legend=True)
            plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

            # simulation
            plots[0].add_data(
                task=f"task_po_emp50_emp50",
                xid="time",
                yid=f"[Cve_emp]",
                label=f"50mg Emp",
                color="black",
            )

            label_map = {
                "emp50": "50 mg Emp",
                "emp50, piog45": "50 mg Emp + Piog 45 mg",
            }

            for intervention in self.interventions:
                # data
                plots[0].add_data(
                    dataset=f"empagliflozin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label_map[intervention],
                    color=self.colors[f"{intervention}"],
                )

            return {
                fig.sid: fig,
            }

if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import RESULTS_PATH_SIMULATION
    run_experiments(Macha2015b, output_dir=RESULTS_PATH_SIMULATION)