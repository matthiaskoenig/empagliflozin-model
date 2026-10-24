"""Simulation experiment of the study Ayoub2017."""

from typing import override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.empagliflozin import MODEL
from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)


class Ayoub2017(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Ayoub2017.

    Six healthy Egyptian men (77.8 kg) take one tablet of 25 mg empagliflozin
    fasted, plasma to 12 hr (Fig6). The error bars of Fig6 are asymmetric, the
    SD of the data is half the distance between the lower and the upper bar
    (Cmax 576 +- 86 nmol/l in Table 1).
    """

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    bodyweight = 77.8  # [kg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig6"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        simulations["po_emp25"] = Simulation(
            time_unit="hr",
            end=13,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg, "mM"),
                "PODOSE_emp": Q(25, "mg"),
            },
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        mappings["fm_po_emp25"] = FitMapping(
            self,
            reference=FitData(
                self,
                dataset="empagliflozin_EMP25",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
            ),
            observable=FitData(
                self,
                task="task_po_emp25",
                xid="time",
                yid="[Cve_emp]",
            ),
            metadata=PKPDMappingMetaData(
                tissue=Tissue.PLASMA,
                route=Route.PO,
                application_form=ApplicationForm.TABLET,
                dosing=Dosing.SINGLE,
                health=Health.HEALTHY,
                fasting=Fasting.FASTED,
                coadministration=Coadministration.NONE,
            ),
        )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig6",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        # simulation
        plots[0].add_data(
            task="task_po_emp25",
            xid="time",
            yid="[Cve_emp]",
            label="25 mg Emp",
            color="black",
        )
        # data
        plots[0].add_data(
            dataset="empagliflozin_EMP25",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="25 mg Emp",
            color="black",
        )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Ayoub2017, output_dir=Ayoub2017.__name__)
