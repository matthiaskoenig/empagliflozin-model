from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.experiments.base_experiment import EmpagliflozinSimulationExperiment
from pkdb_models.models.empagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, EmpagliflozinMappingMetaData, Coadministration
from pkdb_models.models.empagliflozin.helpers import run_experiments


class Ayoub2017(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Ayoub2017."""

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects, assuming 5 mM)
    bodyweight = 77.8  # [kg]
    # gfr = EmpagliflozinSimulationExperiment.gfr_reference  # [ml/min] (healthy subjects, assuming 100 ml/min)

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig6"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                # if label.startswith("empagliflozin_"):
                #    dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        tcsims[f"po_emp25"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=13 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "[KI__fpg]": Q_(self.fpg, "mM"),
                    # "KI__f_renal_function": Q_(self.gfr, "dimensionless")
                    "PODOSE_emp": Q_(25, "mg"),
                },
            )]
        )
        # console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        mappings[f"fm_po_emp25"] = FitMapping(
            self,
            reference=FitData(
                self,
                dataset=f"empagliflozin_EMP25",
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
                dosing=Dosing.SINGLE,
                health=Health.HEALTHY,
                fasting=Fasting.NR,
                coadministration=Coadministration.NONE,
            ),
        )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig6",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        # simulation
        plots[0].add_data(
            task=f"task_po_emp25",
            xid="time",
            yid=f"[Cve_emp]",
            label=f"25 mg Emp",
            color="black",
        )
        # data
        plots[0].add_data(
            dataset=f"empagliflozin_EMP25",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=f"25 mg Emp",
            color="black",
        )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Ayoub2017, output_dir=Ayoub2017.__name__)
