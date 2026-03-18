from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.experiments.base_experiment import EmpagliflozinSimulationExperiment
from pkdb_models.models.empagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, EmpagliflozinMappingMetaData, Coadministration
from pkdb_models.models.empagliflozin.helpers import run_experiments


class Laffel2018(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Laffel2018."""

    groups = {
        "EMP5":  {"dose": 5,  "bw": 90.0,  "fpg": 8.5, "gfr": 178.3},
        "EMP10": {"dose": 10, "bw": 111.0, "fpg": 8.6, "gfr": 162.3},
        "EMP25": {"dose": 25, "bw": 91.1,  "fpg": 6.4, "gfr": 157.3},
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # if label.startswith("EMP"):
                #    dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[f"{label}"] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for group, params in self.groups.items():
            tcsims[f"po_{group}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(params["bw"], "kg"),
                        "[KI__fpg]": Q_(params["fpg"], "mM"),
                        "KI__f_renal_function": Q_(params["gfr"] / 100, "dimensionless"),
                        "PODOSE_emp": Q_(params["dose"], "mg"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group, params in self.groups.items():
            mappings[f"fm_po_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{group}_concentration",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_{group}", xid="time", yid="[Cve_emp]",
                ),
                metadata=EmpagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.T2DM,
                    fasting=Fasting.FED,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__} (T2DM, Pediatric)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        for group, params in self.groups.items():
            dose = params["dose"]
            color = self.dose_colors[dose]
            label = f"{dose} mg Emp"

            # simulation
            plots[0].add_data(
                task=f"task_po_{group}",
                xid="time",
                yid="[Cve_emp]",
                label=label,
                color=color,
            )
            # data
            plots[0].add_data(
                dataset=f"{group}_concentration",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label,
                color=color,
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import RESULTS_PATH_SIMULATION
    run_experiments(Laffel2018, output_dir=RESULTS_PATH_SIMULATION)