"""Simulation experiment of the study vanderAartvanderBeek2020."""

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


class vanderAartvanderBeek2020(EmpagliflozinSimulationExperiment):
    """Simulation experiment of vanderAartvanderBeek2020.

    Individual plasma concentrations of six patients with T2DM and elevated
    albuminuria after the first dose of 10 mg empagliflozin (section 2.7 of
    the method paper: samples at 0, 15, 30, 45, 60, 90, 120, 240 and 360 min
    and 24 hr after the first dose of a once daily treatment). Body weight,
    renal function, glucose and food are not reported: the reference body
    weight and the healthy renal function of the model are used. The 24 hr
    values of S1 (138 µg/l) and S6 (85 µg/l) lie above their values at 4 to
    6 hr, implausible for a single dose of 10 mg (the other subjects 5-15
    µg/l); the values are kept as drawn in Fig3.
    """

    #: fasting plasma glucose, assumed (T2DM, not reported)
    fpg = 10  # [mM]
    subjects: ClassVar[list[str]] = [f"S{k + 1}" for k in range(6)]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig3")
        for label, df_label in df.groupby("label"):
            label = str(label)
            dset = DataSet.from_df(df_label, self.ureg)
            # unit conversion
            if label.startswith("empagliflozin_"):
                dset.unit_conversion("value", 1 / self.Mr.emp)

            dsets[label] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        # the first dose of the treatment
        simulations["po_emp10"] = Simulation(
            time_unit="hr",
            end=30,
            preinit_changes={
                **self.default_changes(),
                "[KI__fpg]": Q(self.fpg, "mM"),
                "PODOSE_emp": Q(10, "mg"),
            },
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
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
                    self,
                    task="task_po_emp10",
                    xid="time",
                    yid="[Cve_emp]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.T2DM,
                    fasting=Fasting.NR,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            name=f"{self.__class__.__name__} (T2DM)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-1, max=26),
            legend=True,
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)

        # simulation
        plots[0].add_data(
            task="task_po_emp10",
            xid="time",
            yid="[Cve_emp]",
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

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(
        vanderAartvanderBeek2020, output_dir=vanderAartvanderBeek2020.__name__
    )
