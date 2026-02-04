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


class Kim2023(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Kim2023."""

    colors = {
        "EP25": "black",
        "EV5, EP25": "tab:blue",
    }
    interventions = list(colors.keys())

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("empagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)

                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        tc0 = Timecourse(
            start=0,
            end=24 * 60,  # [min]
            steps=500,
            changes={
                **self.default_changes(),
                # "BW": Q_(self.bodyweight, "kg"),
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
            end=25 * 60,  # [min]
            steps=500,
            changes={
                "PODOSE_emp": Q_(25, "mg"),
            },
        )

        tcsims[f"po_emp25"] = TimecourseSim(
            [tc0] + [tc1 for _ in range(3)] + [tc2],
            time_offset=-4 * 24 * 60,
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}
        for intervention in self.interventions:
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
                    self, task=f"task_po_emp25", xid="time", yid=f"[Cve_emp]",
                ),
                metadata=EmpagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.EVOGLIPTIN if "EV5" in intervention else Coadministration.NONE
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
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[0].xaxis.min = -5
        plots[0].xaxis.max = 25

        label_map = {
            "EP25": "25 mg Emp",
            "EV5, EP25": "25 mg Emp + Evo",
        }

        # simulation
        plots[0].add_data(
            task="task_po_emp25",
            xid="time",
            yid="[Cve_emp]",
            label="25 mg Emp",
            color="black",
        )

        # data
        for intervention in self.interventions:
            plots[0].add_data(
                dataset=f"empagliflozin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label_map[intervention],
                color=self.colors[intervention],
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Kim2023, output_dir=Kim2023.__name__)
