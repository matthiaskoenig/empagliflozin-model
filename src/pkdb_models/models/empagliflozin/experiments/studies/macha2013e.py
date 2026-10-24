"""Simulation experiment of the study Macha2013e."""

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


class Macha2013e(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2013e.

    Interaction of multiple doses of empagliflozin with warfarin in healthy
    subjects.
    """

    interventions: ClassVar[list[str]] = ["emp", "emp, war"]

    colors: ClassVar[dict[str, str]] = {
        "emp": "black",
        "emp, war": "tab:blue",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig2")
        for label, df_label in df.groupby("label"):
            label = str(label)
            dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)
        return dsets

    #: the daily doses before the profile and the doses of the arm: arm A
    #: (empagliflozin alone) has its profile after the fifth and last dose
    #: (day 5); arm B continues empagliflozin on days 6-12 with warfarin on
    #: day 6, its profile follows the sixth dose and the troughs reach day 13
    #: (24 hr after the twelfth dose). The times of the data are relative to
    #: the dose of the profile.
    doses: ClassVar[dict[str, tuple[int, int]]] = {
        "emp": (4, 5),
        "emp, war": (5, 12),
    }

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        dose = Q(25, "mg")

        for intervention in self.interventions:
            n_before, n_doses = self.doses[intervention]
            simulations[self.to_sid(f"po_emp25_{intervention}")] = Simulation(
                time_unit="hr",
                end=24 * n_doses + 1,
                preinit_changes={
                    **self.default_changes(),
                    # the weight is not reported (BMI 25.4 kg/m2): the
                    # reference of 75 kg
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),  # healthy reference value
                    "KI__f_renal_function": Q(1.0, "dimensionless"),  # healthy
                    "PODOSE_emp": dose,
                },
                changes=[
                    self.add_doses(
                        [24 * k for k in range(1, n_doses)], {"PODOSE_emp": dose}
                    )
                ],
                time_shift=-24 * n_before,
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
                    yid_sd=None,
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
                    # the publication reports no fasting
                    fasting=Fasting.NR,
                    coadministration=Coadministration.NONE
                    if intervention == "emp"
                    else Coadministration.WARFARIN,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        # the time is relative to the dose of the profile (day 5 or 6)
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-24, max=170),
            legend=True,
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        label_map = {
            "emp": "25 mg Emp",
            "emp, war": "25 mg Emp + War",
        }

        for intervention in self.interventions:
            plots[0].add_data(
                task=self.to_sid(f"task_po_emp25_{intervention}"),
                xid="time",
                yid="[Cve_emp]",
                label=label_map[intervention],
                color=self.colors[intervention],
            )
            plots[0].add_data(
                dataset=self.to_sid(f"empagliflozin_{intervention}"),
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=label_map[intervention],
                color=self.colors[intervention],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Macha2013e, output_dir=Macha2013e.__name__)
