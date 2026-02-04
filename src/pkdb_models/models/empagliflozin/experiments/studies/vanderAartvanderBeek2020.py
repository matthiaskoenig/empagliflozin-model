from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.experiments.base_experiment import EmpagliflozinSimulationExperiment
from pkdb_models.models.empagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, EmpagliflozinMappingMetaData
from pkdb_models.models.empagliflozin.helpers import run_experiments


class vanderAartvanderBeek2020(EmpagliflozinSimulationExperiment):
    """Simulation experiment of vanderAartvanderBeek2020.

    T2DM
    """
    fpg = 10  # [mM]
    subjects = [f"S{k+1}" for k in range(6)]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("empagliflozin_"):
                    dset.unit_conversion("value", 1 / self.Mr.emp)

                dsets[f"{label}"] = dset

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
                "[KI__fpg]": Q_(self.fpg, "mM"),
                "PODOSE_emp": Q_(10, "mg"),
            },
        )
        tc1 = Timecourse(
            start=0,
            end=24 * 60,  # [min]
            steps=500,
            changes={
                "PODOSE_emp": Q_(10, "mg"),
            },
        )

        tcsims[f"po_emp10"] = TimecourseSim(
            [tc0] + [tc1 for _ in range(9)],
            time_offset=-9 * 24 * 60,
        )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}
        for subject in self.subjects:
            mappings[f"fm_emp10_{subject}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"empagliflozin_EMP_{subject}",
                    xid="time",
                    yid="value",
                    yid_sd=None,
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_emp10", xid="time", yid=f"[Cve_emp]",
                ),
                metadata=EmpagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.T2DM,
                    fasting=Fasting.NR,
                ),
            )

        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig3",
            name=f"{self.__class__.__name__} (T2DM)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[0].xaxis.min = -3
        plots[0].xaxis.max = 24

        # simulation
        plots[0].add_data(
            task=f"task_po_emp10",
            xid="time",
            yid=f"[Cve_emp]",
            label="10 mg Emp",
            color="black",
        )
        # data
        for subject in self.subjects:
            plots[0].add_data(
                dataset=f"empagliflozin_EMP_{subject}",
                xid="time",
                yid="value",
                yid_sd=None,
                count="count",
                label=f"10 mg Emp ({subject})",
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(vanderAartvanderBeek2020, output_dir=vanderAartvanderBeek2020.__name__)
