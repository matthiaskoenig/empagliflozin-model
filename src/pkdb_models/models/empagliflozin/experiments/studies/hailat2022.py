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


class Hailat2022(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Friedrich2013."""

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects, assuming 5 mM)
    bodyweight = 79.5  # [kg]
    gfr = EmpagliflozinSimulationExperiment.gfr_healthy  # [ml/min] (healthy subjects, assuming 100 ml/min)

    interventions = [
        "fasting_test",
        "fasting_reference",
        "fed_test",
        "fed_reference",
    ]
    colors = {
        "fasting_test": "tab:blue",
        "fasting_reference": "tab:orange",
        "fed_test": "tab:purple",
        "fed_reference": "tab:green",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                # dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        tcsims[f"po_emp10"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=52 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "[KI__fpg]": Q_(self.fpg, "mM"),
                    # "KI__f_renal_function": Q_(self.gfr/100, "dimensionless"),
                    "PODOSE_emp": Q_(10, "mg"),
                },
            )]
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for intervention in self.interventions:
            mappings[f"fm_{intervention}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_emp10", xid="time", yid="[Cve_emp]",
                ),
                metadata=EmpagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED if "fasting" in intervention else Fasting.FED,
                    coadministration=Coadministration.NONE,
                    outlier=True  # probably incorrect units; Cmax 75 µg/l = 75/450.909 nmole/l = 0.166 nmol/l
                ),
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            name=f"{self.__class__.__name__} (Healthy)",
        )

        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        # simulation
        plots[0].add_data(
            task=f"task_po_emp10",
            xid="time",
            yid="[Cve_emp]",
            label="10 mg Emp",
            color="black",
        )

        for intervention in self.interventions:
            # data
            plots[0].add_data(
                dataset=intervention,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"10 mg Emp ({intervention})",
                color=self.colors[intervention],
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Hailat2022, output_dir=Hailat2022.__name__)
