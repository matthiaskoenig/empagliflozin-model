from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.experiments.base_experiment import EmpagliflozinSimulationExperiment
from pkdb_models.models.empagliflozin.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, Coadministration, EmpagliflozinMappingMetaData
)
from pkdb_models.models.empagliflozin.helpers import run_experiments


class Macha2013d(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2013d."""

    bodyweight_single = 71  # [kg]
    bodyweight_multi = 69   # [kg]

    single_interventions = ["emp25_single", "emp25_single, ver120_single"]
    multi_interventions = ["emp25_multi", "emp25_multi, ram5_multi"]

    colors = {
        "emp25_single": "black",
        "emp25_single, ver120_single": "tab:blue",
        "emp25_multi": "black",
        "emp25_multi, ram5_multi": "tab:orange",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        df = load_pkdb_dataframe("Macha2013_Fig1", data_path=self.data_path)
        for label, df_label in df.groupby("label"):
            dsets[label] = DataSet.from_df(df_label, self.ureg)
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_

        tc_single = TimecourseSim(
            Timecourse(
                start=0,
                end=73 * 60,
                steps=1000,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight_single, "kg"),
                    "[KI__fpg]": Q_(self.fpg_healthy, "mM"),
                    "KI__f_renal_function": Q_(1.0, "dimensionless"),
                    "PODOSE_emp": Q_(25, "mg"),
                },
            )
        )

        bw = Q_(self.bodyweight_multi, "kg")
        fpg = Q_(self.fpg_healthy, "mM")
        tc0 = Timecourse(
            start=0,
            end=24 * 60,
            steps=1000,
            changes={
                **self.default_changes(),
                "BW": bw,
                "[KI__fpg]": fpg,
                "KI__f_renal_function": Q_(1.0, "dimensionless"),
                "PODOSE_emp": Q_(25, "mg"),
            },
        )
        tc1 = Timecourse(
            start=0,
            end=24 * 60,
            steps=1000,
            changes={
                "BW": bw,
                "[KI__fpg]": fpg,
                "PODOSE_emp": Q_(25, "mg"),
            },
        )
        tc2 = Timecourse(
            start=0,
            end=73 * 60,
            steps=1000,
            changes={
                "BW": bw,
                "[KI__fpg]": fpg,
                "PODOSE_emp": Q_(25, "mg"),
            },
        )
        tc_multi = TimecourseSim(
            [tc0] + [tc1 for _ in range(3)] + [tc2],
            time_offset=-4 * 24 * 60,
        )

        return {
            "po_emp25_single": tc_single,
            "po_emp25_multi": tc_multi,
        }

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        all_interventions = [
            (self.single_interventions, "task_po_emp25_single", Dosing.SINGLE),
            (self.multi_interventions,  "task_po_emp25_multi",  Dosing.MULTIPLE),
        ]
        for interventions, task, dosing in all_interventions:
            for intervention in interventions:
                if "ver120" in intervention:
                    coadmin = Coadministration.VERAPAMIL
                elif "ram5" in intervention:
                    coadmin = Coadministration.RAMIPRIL
                else:
                    coadmin = Coadministration.NONE

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
                        self, task=task, xid="time", yid="[Cve_emp]",
                    ),
                    metadata=EmpagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=dosing,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=coadmin,
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        Figure.legend_fontsize = 10
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        plots[0].add_data(
            task="task_po_emp25_single",
            xid="time",
            yid="[Cve_emp]",
            label="25 mg Emp",
            color="black",
        )
        for intervention in self.single_interventions:
            label = "25 mg Emp" if "ver120" not in intervention else "25 mg Emp + Ver"
            plots[0].add_data(
                dataset=f"empagliflozin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label,
                color=self.colors[intervention],
            )

        plots[1].add_data(
            task="task_po_emp25_multi",
            xid="time",
            yid="[Cve_emp]",
            label="25 mg Emp QD",
            color="black",
        )
        for intervention in self.multi_interventions:
            label = "25 mg Emp QD" if "ram5" not in intervention else "25 mg Emp QD + Ram"
            plots[1].add_data(
                dataset=f"empagliflozin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label,
                color=self.colors[intervention],
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import RESULTS_PATH_SIMULATION
    run_experiments(Macha2013d, output_dir=RESULTS_PATH_SIMULATION)