"""Simulation experiment of the study Macha2015b."""

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


class Macha2015b(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2015b.

    Interaction of multiple doses of empagliflozin with pioglitazone in healthy
    subjects (study 1 of the paper, 20 white men, median 86.5 and 81.5 kg in the
    two sequences, Table I). Empagliflozin 50 mg qd for 5 days (A), directly
    followed by empagliflozin 50 mg with pioglitazone 45 mg qd for 7 days (B);
    the profiles are of the last dose of each period. Both arms are simulated
    as 5 doses: with a half-life of 9-12 hr the steady state is reached after
    5 doses, the 7 further days of B do not change it. Pioglitazone does not
    change the exposure of empagliflozin (GMR AUC 1.00, Cmax 0.93, Table II):
    the arm B is the simulation of A without an interaction factor, an
    interaction validation. Fasted overnight for A and C; B not stated,
    assumed fasted as A.
    """

    #: median body weight of the 20 subjects (86.5 and 81.5 kg, Table I)
    bodyweight = 84.0  # [kg]

    interventions: ClassVar[list[str]] = ["emp50", "emp50, piog45"]

    colors: ClassVar[dict[str, str]] = {
        "emp50": "black",
        "emp50, piog45": "tab:blue",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig1")
        for label, df_label in df.groupby("label"):
            label = str(label)
            dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        dose = Q(50, "mg")

        for intervention in self.interventions:
            # a dose every 24 hr for 5 days, the time is relative to the last
            # dose; pioglitazone without an interaction factor (no effect)
            simulations[self.to_sid(f"po_emp50_{intervention}")] = Simulation(
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
            mappings[self.to_sid(f"fm_po_emp50_{intervention}")] = FitMapping(
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
                    task=self.to_sid(f"task_po_emp50_{intervention}"),
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
                    if intervention == "emp50"
                    else Coadministration.PIOGLITAZONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-24), legend=True
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        # simulation
        plots[0].add_data(
            task="task_po_emp50_emp50",
            xid="time",
            yid="[Cve_emp]",
            label="50 mg Emp",
            color="black",
        )

        label_map = {
            "emp50": "50 mg Emp",
            "emp50, piog45": "50 mg Emp + Piog 45 mg",
        }

        for intervention in self.interventions:
            # data
            plots[0].add_data(
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
    MODEL.run_experiments(Macha2015b, output_dir=Macha2015b.__name__)
