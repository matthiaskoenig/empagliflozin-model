"""Simulation experiment of the study Brand2012."""

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


class Brand2012(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Brand2012.

    Sixteen healthy men (78 kg) take 50 mg empagliflozin once daily for 5 days
    (EMP50_D5), directly followed by 5 days of 50 mg empagliflozin with 100 mg
    sitagliptin (EMP50_D10, the tenth dose of empagliflozin). The profile days
    start after a fast of at least 10 hr. Plasma (Fig2, semi-log panel b with
    the pre-dose sample) and urine (Fig3) after the last dose; the error bars
    are SD (2 hr CV 21-31 %, Cmax,ss CV 23.8 % in Table 1).

    Sitagliptin is cleared by the kidneys (OAT3, P-gp) and does not change the
    exposure of empagliflozin: AUCtau,ss ratio 1.10 (1.04-1.17), Cmax,ss 1.08
    (0.97-1.19), fe0-24 17.1 vs 19.3 %. The arm is simulated without an
    interaction factor and is validation data.
    """

    interventions: ClassVar[list[str]] = ["EMP50_D5", "EMP50_D10"]
    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "Aurine_emp": "empagliflozin_urine",
    }

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    bodyweight = 78  # [kg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion
                if label.startswith("empagliflozin_urine_"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        dose = Q(50, "mg")
        for intervention in self.interventions:
            # a dose every 24 hr for 5 or 10 days, the time is relative to the
            # last dose; the urine is collected after the last dose (the amount
            # is reset at every dose)
            n_doses = 5 if intervention.endswith("D5") else 10
            end = 24 * (n_doses - 1) + 25  # [hr]
            simulations[f"po_{intervention}"] = Simulation(
                time_unit="hr",
                end=end,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "[KI__fpg]": Q(self.fpg, "mM"),
                    "PODOSE_emp": dose,
                },
                changes=[
                    Change(
                        [24 * k for k in range(1, n_doses)],
                        {"Aurine_emp": Q(0, "mmole")},
                    ),
                    self.add_doses(
                        [24 * k for k in range(1, n_doses)], {"PODOSE_emp": dose}
                    ),
                ],
                time_shift=-24 * (n_doses - 1),
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.SITAGLIPTIN
                        if "D10" in intervention
                        else Coadministration.NONE,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2_3",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-5, max=25),
            legend=True,
        )
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)

        for kp, (sid, name) in enumerate(self.info.items()):
            # simulation: plotted once
            plots[kp].add_data(
                task="task_po_EMP50_D5",
                xid="time",
                yid=sid,
                label="50 mg Emp",
                color="black",
            )

            # data
            for intervention in self.interventions:
                is_emp_only = intervention.endswith("D5")

                plots[kp].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label="50 mg Emp" if is_emp_only else "50 mg Emp + SIT",
                    color="black" if is_emp_only else "tab:blue",
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Brand2012, output_dir=Brand2012.__name__)
