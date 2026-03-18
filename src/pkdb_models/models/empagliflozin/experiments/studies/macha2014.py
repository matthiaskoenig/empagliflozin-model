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


class Macha2014(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2014."""

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    interventions = ["EMP25", "EMP25, GEM", "EMP10_R", "EMP10_P", "EMP10, RIF", "EMP10, PROB"]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
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

        # Single dosing
        for intervention in self.interventions:
            if "25" in intervention:
                tcsims[f"po_emp25_{intervention}"] = TimecourseSim(
                    Timecourse(
                        start=0,
                        end=73 * 60,  # [min]
                        steps=500,
                        changes={
                            **self.default_changes(),
                            "[KI__fpg]": Q_(self.fpg_healthy, "mM"),  # healthy reference value
                            "KI__f_renal_function": Q_(1.0, "dimensionless"),  # healthy
                            "PODOSE_emp": Q_(25, "mg"),
                        },
                    )
                )

            else:
                tcsims[f"po_emp10_{intervention}"] = TimecourseSim(
                    Timecourse(
                        start=0,
                        end=73 * 60,  # [min]
                        steps=500,
                        changes={
                            **self.default_changes(),
                            # "BW": Q_(self.bodyweights[intervention], "kg"),
                            "[KI__fpg]": Q_(self.fpg_healthy, "mM"),  # healthy reference value
                            "KI__f_renal_function": Q_(1.0, "dimensionless"),  # healthy
                            # FIXME: fasting
                            "PODOSE_emp": Q_(10, "mg"),
                        },
                    )
                )

        # console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}
        for intervention in self.interventions:
            if "25" in intervention:
                mappings[f"fm_emp25_{intervention}"] = FitMapping(
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
                        fasting=Fasting.FASTED)
                )
            else:
                mappings[f"fm_po_emp10_{intervention}"] = FitMapping(
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
                        self, task=f"task_po_emp10_{intervention}", xid="time", yid=f"[Cve_emp]",
                    ),
                    metadata=EmpagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED
                    ),
                )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__} (Healthy)",
            num_cols=2,
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time, min=-24), legend=True)

        # Set y-axis for both plots
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        plots[0].xaxis.min = -1
        plots[1].xaxis.min = -1

        # 10mg simulation
        plots[0].add_data(
            task=f"task_po_emp10_EMP10_R",
            xid="time",
            yid=f"[Cve_emp]",
            label=f"10 mg Emp",
            color="black",
        )

        plots[0].add_data(
            dataset=f"empagliflozin_EMP10_R",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=f"10 mg Emp",
            color="black",
        )

        plots[0].add_data(
            dataset=f"empagliflozin_EMP10_P",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=f"10 mg Emp",
            color="black",
            marker="D",
        )

        plots[0].add_data(
            dataset=f"empagliflozin_EMP10, RIF",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=f"10 mg Emp + Rif",
            color="tab:blue",
        )

        plots[0].add_data(
            dataset=f"empagliflozin_EMP10, PROB",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=f"10 mg Emp + Prob",
            color="tab:orange",
        )

        # 25mg simulation
        plots[1].add_data(
            task=f"task_po_emp25_EMP25",
            xid="time",
            yid=f"[Cve_emp]",
            label=f"25 mg Emp",
            color="black",
        )

        plots[1].add_data(
            dataset=f"empagliflozin_EMP25",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=f"25 mg Emp",
            color="black",
        )

        plots[1].add_data(
            dataset=f"empagliflozin_EMP25, GEM",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=f"25 mg Emp + Gem",
            color="tab:blue",
        )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import RESULTS_PATH_SIMULATION
    run_experiments(Macha2014, output_dir=RESULTS_PATH_SIMULATION)