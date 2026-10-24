"""Simulation experiment of the study Macha2013d."""

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


class Macha2013d(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2013d.

    Interaction of empagliflozin with verapamil (single dose) and ramipril
    (multiple dose) in healthy subjects. The data is curated in the ramipril
    study Macha2013.
    """

    bodyweight_single = 71  # [kg]
    bodyweight_multi = 69  # [kg]

    single_interventions: ClassVar[list[str]] = [
        "emp25_single",
        "emp25_single, ver120_single",
    ]
    multi_interventions: ClassVar[list[str]] = [
        "emp25_multi",
        "emp25_multi, ram5_multi",
    ]

    colors: ClassVar[dict[str, str]] = {
        "emp25_single": "black",
        "emp25_single, ver120_single": "tab:blue",
        "emp25_multi": "black",
        "emp25_multi, ram5_multi": "tab:orange",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        # the data of the study id Macha2013, not the id of the experiment
        df = self.load_dataframe("Fig1", study="Macha2013")
        dsets = {}
        for label, df_label in df.groupby("label"):
            label = str(label)
            dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tc_single = Simulation(
            time_unit="hr",
            end=73,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight_single, "kg"),
                "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                "KI__f_renal_function": Q(1.0, "dimensionless"),
                "PODOSE_emp": Q(25, "mg"),
            },
        )

        # a dose every 24 hr for 5 days, the time is relative to the last dose
        dose = Q(25, "mg")
        tc_multi = Simulation(
            time_unit="hr",
            end=169,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight_multi, "kg"),
                "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                "KI__f_renal_function": Q(1.0, "dimensionless"),
                "PODOSE_emp": dose,
            },
            changes=[self.add_doses([24, 48, 72, 96], {"PODOSE_emp": dose})],
            time_shift=-96,
        )

        return {
            "po_emp25_single": tc_single,
            "po_emp25_multi": tc_multi,
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}

        all_interventions = [
            (self.single_interventions, "task_po_emp25_single", Dosing.SINGLE),
            (self.multi_interventions, "task_po_emp25_multi", Dosing.MULTIPLE),
        ]
        for interventions, task, dosing in all_interventions:
            for intervention in interventions:
                if "ver120" in intervention:
                    coadmin = Coadministration.VERAPAMIL
                elif "ram5" in intervention:
                    coadmin = Coadministration.RAMIPRIL
                else:
                    coadmin = Coadministration.NONE

                mappings[self.to_sid(f"fm_po_emp25_{intervention}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"empagliflozin_{intervention}"),
                        xid="time",
                        yid="mean",
                        yid_sd=None,  # the figure has no SD
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(task),
                        xid="time",
                        yid="[Cve_emp]",
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        # the publication names neither the formulation nor
                        # the fasting state
                        application_form=ApplicationForm.NR,
                        dosing=dosing,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration=coadmin,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        # the profile of day 5, the time is relative to the last dose
        plots[1].set_xaxis(self.label_time, unit=self.unit_time, min=-24, max=73)

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
                dataset=self.to_sid(f"empagliflozin_{intervention}"),
                xid="time",
                yid="mean",
                yid_sd=None,  # the figure has no SD
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
            label = (
                "25 mg Emp QD" if "ram5" not in intervention else "25 mg Emp QD + Ram"
            )
            plots[1].add_data(
                dataset=self.to_sid(f"empagliflozin_{intervention}"),
                xid="time",
                yid="mean",
                yid_sd=None,  # the figure has no SD
                count="count",
                label=label,
                color=self.colors[intervention],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Macha2013d, output_dir=Macha2013d.__name__)
