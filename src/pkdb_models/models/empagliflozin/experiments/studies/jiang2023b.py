"""Simulation experiment of the study Jiang2023b."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
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


class Jiang2023b(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Jiang2023b."""

    #: the test formulation CKD-370 (empagliflozin L-proline cocrystal) and the
    #: reference Jardiance, bioequivalent, in the same 27 subjects (crossover)
    interventions: ClassVar[list[str]] = ["test", "ref"]

    colors: ClassVar[dict[str, str]] = {
        "test": "tab:blue",
        "ref": "black",
    }

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    bodyweight = 69.4  # [kg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("empagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        simulations["po_emp25"] = Simulation(
            time_unit="hr",
            end=50,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg, "mM"),
                "KI__f_renal_function": Q(1.0, "dimensionless"),  # healthy
                "PODOSE_emp": Q(25, "mg"),
            },
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
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
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        # linear (A) and semi-log (B) as in the paper
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(
            self.label_emp_plasma, unit=self.unit_emp, scale="log", min=0.002
        )

        for plot in plots:
            # simulation
            plot.add_data(
                task="task_po_emp25",
                xid="time",
                yid="[Cve_emp]",
                label="25 mg Emp",
                color="black",
            )
            # data
            for intervention in self.interventions:
                plot.add_data(
                    dataset=f"empagliflozin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"25 mg Emp ({intervention})",
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Jiang2023b, output_dir=Jiang2023b.__name__)
