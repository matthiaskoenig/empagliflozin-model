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


class Chen2020(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Chen2020."""

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy, actual value not reported)
    bodyweights = {
        "Fasting_T": 62.93,
        "Fasting_R": 62.93,
        "Fed_T": 60.99,
        "Fed_R": 60.99,
    }  # [kg]

    colors = {
        "Fasting_T": "black",
        "Fasting_R": "tab:blue",
        "Fed_T": "tab:green",
        "Fed_R": "tab:orange",
    }
    groups = list(bodyweights.keys())

    gfr = EmpagliflozinSimulationExperiment.gfr_healthy  # [ml/min] (Healthy)

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.endswith("EMP10"):
                   dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        # Single dosing
        for group, bodyweight in self.bodyweights.items():
            tcsims[f"po_emp10_{group}"] = TimecourseSim(
                [Timecourse(
                start=0,
                end=73 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(bodyweight, "kg"),
                    "[KI__fpg]": Q_(self.fpg_healthy, "mM"),  # healthy reference value
                    # "KI__f_renal_function": Q_(self.gfr, "dimensionless")
                    "PODOSE_emp": Q_(10, "mg"),
                },
            )]
            )


        # console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            mappings[f"fm_po_emp10_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{group}_EMP10",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_emp10_{group}", xid="time", yid=f"[Cve_emp]",
                ),
                metadata=EmpagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED if group.startswith("Fasting") else Fasting.FED,
                    coadministration=Coadministration.NONE,
                ),
            )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp, min=-0.05)

        for group in self.groups:
            # simulation
            plots[0].add_data(
                task=f"task_po_emp10_{group}",
                xid="time",
                yid=f"[Cve_emp]",
                label=f"10 mg Emp ({group})",
                color=self.colors[group],
            )
            # data
            plots[0].add_data(
                dataset=f"{group}_EMP10",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"10 mg Emp ({group})",
                color=self.colors[group]
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Chen2020, output_dir=Chen2020.__name__)
