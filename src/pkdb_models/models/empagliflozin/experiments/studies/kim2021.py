"""Simulation experiment of the study Kim2021."""

from typing import ClassVar, override

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


class Kim2021(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Kim2021."""

    bodyweight = 73.1  # [kg]
    # empagliflozin (E) alone and with lobeglitazone (L)
    interventions: ClassVar[list[str]] = ["E", "E, L"]

    colors: ClassVar[dict[str, str]] = {
        "E": "black",
        "E, L": "tab:blue",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.startswith("empagliflozin"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        dose = Q(25, "mg")

        for intervention in self.interventions:
            # a dose every 24 hr for 5 days, the time is relative to the last
            # dose
            simulations[self.to_sid(f"po_emp25_{intervention}")] = Simulation(
                time_unit="hr",
                end=126,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),  # healthy reference value
                    "KI__f_renal_function": Q(1.0, "dimensionless"),  # healthy
                    "PODOSE_emp": dose,
                },
                changes=[self.add_doses([24, 48, 72, 96], {"PODOSE_emp": dose})],
                time_shift=-96,
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            mappings[self.to_sid(f"fm_po_emp25_{intervention}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"empagliflozin_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_emp25_{intervention}"),
                    xid="time",
                    yid="[Cve_emp]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,  # overnight fast
                    coadministration=Coadministration.NONE
                    if intervention == "E"
                    else Coadministration.LOBEGLITAZONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-24), legend=True
        )
        # linear and semi-log (the inset of the paper)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(
            self.label_emp_plasma, unit=self.unit_emp, scale="log", min=0.02
        )

        label_map = {
            "E": "25 mg Emp",
            "E, L": "25 mg Emp + Lob",
        }

        for plot in plots:
            # simulation, the same for both arms (no interaction)
            plot.add_data(
                task="task_po_emp25_E",
                xid="time",
                yid="[Cve_emp]",
                label="25 mg Emp",
                color="black",
            )

            # data
            for intervention in self.interventions:
                plot.add_data(
                    dataset=self.to_sid(f"empagliflozin_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label_map[intervention],
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kim2021, output_dir=Kim2021.__name__)
