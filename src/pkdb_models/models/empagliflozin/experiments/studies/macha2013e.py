from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)
from pkdb_models.models.empagliflozin.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, Coadministration, EmpagliflozinMappingMetaData
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.helpers import run_experiments


class Macha2013e(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2013e."""

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    # fpg = 7  # [mM] (healthy subjects, application with food)
    interventions = ["emp", "emp, war"]

    colors = {
        "emp": "black",
        "emp, war": "tab:blue",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                #if label.startswith("empagliflozin"):
                  # dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[label] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        # Single dosing
        for intervention in self.interventions:
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    #"BW": Q_(self.bodyweight, "kg"),
                    "[KI__fpg]": Q_(self.fpg_healthy, "mM"),  # healthy reference value
                    "KI__f_renal_function": Q_(1.0, "dimensionless"), # healthy
                    # healthy reference value
                    # FIXME: fasting
                    "PODOSE_emp": Q_(25, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_emp": Q_(25, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=170 * 60,  # [min]
                steps=2000,
                changes={
                    "PODOSE_emp": Q_(25, "mg"),
                },
            )

            if intervention == "emp":
                n_repeats = 3
            elif intervention == "emp, war":
                n_repeats = 5

            tcsims[f"po_emp25_{intervention}"] = TimecourseSim(
               [tc0] + [tc1 for _ in range(n_repeats)] + [tc2],
               time_offset=-(n_repeats+1)*24*60
            )
            # console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
            mappings = {}
            for intervention in self.interventions:
                mappings[f"fm_po_emp25_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"empagliflozin_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
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
                        coadministration=Coadministration.NONE if intervention == "emp" else Coadministration.WARFARIN,
                    ),
                )
            # console.print(mappings)
            return mappings

    def figures(self) -> Dict[str, Figure]:
            fig = Figure(
                experiment=self,
                sid="Fig3",
                name=f"{self.__class__.__name__} (Healthy)",
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
            plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

            plots[0].add_data(
                task=f"task_po_emp25_emp",
                xid="time",
                yid=f"[Cve_emp]",
                label=f"25 mg Emp",
                color="black",
            )

            label_map = {
                "emp": "25 mg Emp",
                "emp, war": "25 mg Emp + War",
            }

            for intervention in self.interventions:
                plots[0].add_data(
                    dataset=f"empagliflozin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=label_map[intervention],
                    color=self.colors[intervention],
                )

            return {
                fig.sid: fig,
            }


if __name__ == "__main__":
    run_experiments(Macha2013e, output_dir=Macha2013e.__name__)
