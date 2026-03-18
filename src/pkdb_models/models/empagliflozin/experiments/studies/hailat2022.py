from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)
from pkdb_models.models.empagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, EmpagliflozinMappingMetaData, Coadministration

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.helpers import run_experiments


class Hailat2022(EmpagliflozinSimulationExperiment):
    """
    Simulation experiment of Hailat2022.
    Outlier.
    Probably wrong units.
    Most likely µg/l instead of nmol/l (has been updated accordingly).
    """

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy
    bodyweight = 79.5  # [kg]
    gfr = EmpagliflozinSimulationExperiment.gfr_healthy

    interventions = [
        "fasting_test",
        "fasting_reference",
        "fed_test",
        "fed_reference",
    ]

    symbols = {
        "test": "s",       # square
        "reference": "D",  # diamond
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[f"{label}"] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for condition in ["fasted", "fed"]:
            tcsims[f"po_emp10_{condition}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=52 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "[KI__fpg]": Q_(self.fpg, "mM"),
                        "GU__f_absorption": Q_(self.fasting_map[condition], "dimensionless"),
                        "PODOSE_emp": Q_(10, "mg"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            condition = "fasted" if "fasting" in intervention else "fed"
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
                    self, task=f"task_po_emp10_{condition}", xid="time", yid="[Cve_emp]",
                ),
                metadata=EmpagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED if "fasting" in intervention else Fasting.FED,
                    coadministration=Coadministration.NONE,
                    outlier=True,
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

        # simulations
        for condition in ["fasted", "fed"]:
            plots[0].add_data(
                task=f"task_po_emp10_{condition}",
                xid="time",
                yid="[Cve_emp]",
                label=f"10 mg Emp ({condition})",
                color=self.fasting_colors[condition],
            )

        # data
        for intervention in self.interventions:
            condition = "fasted" if "fasting" in intervention else "fed"
            formulation = intervention.split("_")[-1]  # "test" or "reference"
            plots[0].add_data(
                dataset=intervention,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"10 mg Emp ({intervention})",
                color=self.fasting_colors[condition],
                marker=self.symbols[formulation],
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import RESULTS_PATH_SIMULATION
    run_experiments(Hailat2022, output_dir=RESULTS_PATH_SIMULATION)