"""Simulation experiment of the study Laffel2018."""

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


class Laffel2018(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Laffel2018 (pediatric T2DM)."""

    #: dose [mg], body weight [kg], fasting plasma glucose [mM] and eGFR
    #: [ml/min/1.73 m2] of the groups (Table 1); the renal function is the
    #: eGFR normalized to the body surface over the GFR of healthy adults (the
    #: GFR of the model scales with the body weight, the absolute GFR of the
    #: adolescents of 90-111 kg would be about 210-260 ml/min)
    groups: ClassVar[dict[str, dict[str, float]]] = {
        "EMP5": {"dose": 5, "bw": 90.0, "fpg": 8.5, "gfr": 178.3},
        "EMP10": {"dose": 10, "bw": 111.0, "fpg": 8.6, "gfr": 162.3},
        "EMP25": {"dose": 25, "bw": 91.1, "fpg": 6.4, "gfr": 157.3},
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig1")
        for label, df_label in df.groupby("label"):
            label = str(label)
            dsets[label] = DataSet.from_df(df_label, self.ureg)
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for group, params in self.groups.items():
            simulations[f"po_{group}"] = Simulation(
                time_unit="hr",
                end=50,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(params["bw"], "kg"),
                    "[KI__fpg]": Q(params["fpg"], "mM"),
                    "KI__f_renal_function": Q(
                        params["gfr"] / self.gfr_healthy, "dimensionless"
                    ),
                    "PODOSE_emp": Q(params["dose"], "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            mappings[f"fm_laffel2018_po_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{group}_concentration",
                    xid="time",
                    yid="mean",
                    yid_sd=None,  # the figure has no SD
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_po_{group}",
                    xid="time",
                    yid="[Cve_emp]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.T2DM,
                    # dosed in the morning before breakfast
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=2,
            name=f"{self.__class__.__name__} (T2DM, Pediatric)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        # linear (a) and semi-log (b) as in the paper
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(
            self.label_emp_plasma, unit=self.unit_emp, scale="log", min=1e-3
        )

        for plot in plots:
            for group, params in self.groups.items():
                dose = params["dose"]
                color = self.dose_colors[dose]
                label = f"{dose} mg Emp"

                # simulation
                plot.add_data(
                    task=f"task_po_{group}",
                    xid="time",
                    yid="[Cve_emp]",
                    label=label,
                    color=color,
                )
                # data
                plot.add_data(
                    dataset=f"{group}_concentration",
                    xid="time",
                    yid="mean",
                    yid_sd=None,  # the figure has no SD
                    count="count",
                    label=label,
                    color=color,
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Laffel2018, output_dir=Laffel2018.__name__)
