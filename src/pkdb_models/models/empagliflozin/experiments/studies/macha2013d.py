from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.experiments.base_experiment import EmpagliflozinSimulationExperiment
from pkdb_models.models.empagliflozin.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, Coadministration, EmpagliflozinMappingMetaData
)
from pkdb_models.models.empagliflozin.helpers import run_experiments


class Macha2013d(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2013d."""

    #groups = ["emp_verapamil", "emp_ramipril"]

    bodyweights = {
          "emp25_single": 71,
          "emp25_single, ver120_single": 71,
          "emp25_multi": 69,
          "emp25_multi, ram5_multi": 69
           }
    #bodyweights = {
    #   "emp_verapamil": 71,
    #   "emp_ramipril": 69}  # [kg]

    colors = {
        "emp25_single": "black",
        "emp25_multi": "tab:blue",
        "emp25_single, ver120_single": "tab:orange",
        "emp25_multi, ram5_multi": "tab:green"
    }

    interventions = ["emp25_single", "emp25_single, ver120_single", "emp25_multi", "emp25_multi, ram5_multi"]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                #if label.startswith("empagliflozin"):
                   #dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[label] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        # Single dosing
        for intervention in self.interventions:
            tcsims[f"po_emp25_{intervention}"] = TimecourseSim(
                Timecourse(
                    start=0,
                    end=73 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweights[intervention], "kg"),
                        "[KI__fpg]": Q_(self.fpg_healthy, "mM"),  # healthy reference value
                        "KI__f_renal_function": Q_(1.0, "dimensionless"),  # healthy
                        # FIXME: fasting
                        "PODOSE_emp": Q_(25, "mg"),
                    },
                )
            )

        # console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
            mappings = {}
            for intervention in self.interventions:
                if "single" in intervention:
                    mappings[f"fm_po_emp25_{intervention}"] = FitMapping(
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
                            self, task=f"task_po_emp25_{intervention}", xid="time", yid=f"[Cve_emp]",
                        ),
                        metadata=EmpagliflozinMappingMetaData(
                            tissue=Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY,
                            fasting=Fasting.FASTED,  # overnight fast
                            coadministration=Coadministration.NONE
                            #coadministration=Coadministration.NONE if "ver120" not in intervention else Coadministration.VERAPAMIL,
                        ),
                    )

                else:
                    mappings[f"fm_po_emp25_{intervention}"] = FitMapping(
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
                            self, task=f"task_po_emp25_{intervention}", xid="time", yid=f"[Cve_emp]",
                        ),
                        metadata=EmpagliflozinMappingMetaData(
                            tissue=Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.MULTIPLE,
                            health=Health.HEALTHY,
                            fasting=Fasting.FASTED,  # overnight fast
                            coadministration=Coadministration.NONE
                            # coadministration=Coadministration.NONE if "ram5" not in intervention else Coadministration.RAMIPRIL,
                        ),
                    )

            # console.print(mappings)
            return mappings


    def figures(self) -> Dict[str, Figure]:
            fig = Figure(
                experiment=self,
                sid="Fig1",
                name=f"{self.__class__.__name__} (Healthy)",
            )
            Figure.legend_fontsize = 10
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
            plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

            for intervention in self.interventions:
                # clean intervention label for legend
                suffix = intervention.replace("emp25_", "").strip()
                suffix = suffix.lstrip(", ").strip()
                label = "25 mg Emp" if suffix == "" else f"25 mg Emp ({suffix})"

                # simulation
                plots[0].add_data(
                    task=f"task_po_emp25_{intervention}",
                    xid="time",
                    yid=f"[Cve_emp]",
                    label=label,
                    color=self.colors[f"{intervention}"],
                )
                # data
                plots[0].add_data(
                    dataset=f"empagliflozin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=self.colors[f"{intervention}"],
                )

            return {
                fig.sid: fig,
            }


if __name__ == "__main__":
    run_experiments(Macha2013d, output_dir=Macha2013d.__name__)
