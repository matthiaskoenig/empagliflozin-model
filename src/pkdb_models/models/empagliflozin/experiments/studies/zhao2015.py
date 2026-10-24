"""Simulation experiment of the study Zhao2015."""

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
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.empagliflozin import MODEL
from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)


class Zhao2015(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Zhao2015.

    Single and multiple doses of 10 mg and 25 mg empagliflozin once daily in
    Chinese patients with T2DM (n=9 per dose, eGFR >= 80 ml/min/1.73 m2,
    background metformin, sulfonylurea or acarbose without an interaction).
    A single dose on day 1 (profile and urine 0-24 hr, `SING`), then once
    daily on days 3 to 9 (7 doses at 48 to 192 hr, profile and urine 192-216
    hr at steady state, `MULTI`), fasted overnight; one simulation per dose.
    Body weight and fasting plasma glucose are the means of Table I. The
    exposure is 1.4-1.9 times that of white patients (as in the Japanese
    studies).
    """

    bodyweights: ClassVar[dict[str, float]] = {
        "MULTI10": 74.0,  # [kg]
        "MULTI25": 67.0,  # [kg]
    }
    #: fasting plasma glucose, 166.67 and 160.56 mg/dl (Table I)
    fpgs: ClassVar[dict[str, float]] = {
        "MULTI10": 166.67 / 18,  # [mM]
        "MULTI25": 160.56 / 18,  # [mM]
    }
    doses: ClassVar[dict[str, int]] = {
        "MULTI10": 10,
        "MULTI25": 25,
    }
    interventions: ClassVar[list[str]] = list(bodyweights.keys())
    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
    }
    #: the data of day 1 (single dose) and of day 9 (steady state)
    days: ClassVar[dict[str, Dosing]] = {
        "SING": Dosing.SINGLE,
        "MULTI": Dosing.MULTIPLE,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab2A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dsets[label] = DataSet.from_df(df_label, self.ureg)
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        for intervention, dose in self.doses.items():
            # doses at 0 hr and every 24 hr from 48 to 192 hr, the urine is
            # collected from dose to dose
            simulations[intervention] = Simulation(
                time_unit="hr",
                end=267,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[intervention], "kg"),
                    "[KI__fpg]": Q(self.fpgs[intervention], "mM"),
                    "PODOSE_emp": Q(dose, "mg"),
                    "Aurine_emp": Q(0, "mmole"),
                },
                changes=[
                    Change(
                        [48 + 24 * k for k in range(7)], {"Aurine_emp": Q(0, "mmole")}
                    ),
                    self.add_doses(
                        [48 + 24 * k for k in range(7)], {"PODOSE_emp": Q(dose, "mg")}
                    ),
                ],
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}

        for intervention, dose in self.doses.items():
            for sid, name in self.info.items():
                for day, dosing in self.days.items():
                    mappings[f"fm_{name}_{day}{dose}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"{name}_{day}{dose}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=f"task_{intervention}",
                            xid="time",
                            yid=sid,
                        ),
                        metadata=PKPDMappingMetaData(
                            tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=dosing,
                            health=Health.T2DM,
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
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (T2DM)",
            # one legend per row, right of the grid (the panels show the same
            # doses)
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        # plasma of day 1 (single dose) and day 9 (steady state), urine
        plots[0].set_xaxis(self.label_time, unit=self.unit_time, min=-2, max=50)
        plots[1].set_xaxis(self.label_time, unit=self.unit_time, min=190, max=266)
        plots[2].set_xaxis(self.label_time, unit=self.unit_time, max=225)
        for kp in [0, 1]:
            plots[kp].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[2].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)

        panels = [
            ("[Cve_emp]", "empagliflozin", ["SING"], plots[0]),
            ("[Cve_emp]", "empagliflozin", ["MULTI"], plots[1]),
            ("Aurine_emp", "empagliflozin_urine", ["SING", "MULTI"], plots[2]),
        ]
        for intervention, dose in self.doses.items():
            color = self.dose_colors[dose]
            dose_label = f"{dose} mg Emp"
            for sid, name, days, plot in panels:
                # simulation
                plot.add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=dose_label,
                    color=color,
                )
                # data (one legend entry per dose, the urine has the same
                # markers)
                for kd, day in enumerate(days):
                    plot.add_data(
                        dataset=f"{name}_{day}{dose}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=dose_label if kd == 0 and sid == "[Cve_emp]" else None,
                        color=color,
                        linestyle="" if "urine" in name else "--",
                    )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Zhao2015, output_dir=Zhao2015.__name__)
