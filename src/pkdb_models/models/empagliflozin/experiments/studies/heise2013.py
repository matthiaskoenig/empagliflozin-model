"""Simulation experiment of the study Heise2013."""

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


class Heise2013(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Heise2013."""

    #: fasting plasma glucose (FPG) at baseline, Table 1 [mM]
    fpgs: ClassVar[dict[str, float]] = {
        "Placebo": 156.9 / 18,
        "EMP2.5": 144.1 / 18,
        "EMP10": 150.1 / 18,
        "EMP25": 142.2 / 18,
        "EMP100": 164.3 / 18,
    }
    #: weighted mean daily glucose (MDG) at baseline (day -2), Table 7 [mM];
    #: the MDG falls by 1.4-2.1 mM until day 8, the model keeps the baseline
    mdgs: ClassVar[dict[str, float]] = {
        "Placebo": 169 / 18,
        "EMP2.5": 161 / 18,
        "EMP10": 160 / 18,
        "EMP25": 159 / 18,
        "EMP100": 185 / 18,
    }
    #: creatinine clearance (Cockcroft-Gault) at baseline, Table 1 [ml/min];
    #: with the total body weight (BMI 27.6-32.9) it overestimates the GFR
    crcls: ClassVar[dict[str, float]] = {
        "Placebo": 103.8,
        "EMP2.5": 125.3,
        "EMP10": 78.0,
        "EMP25": 146.8,
        "EMP100": 117.2,
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "Placebo": 100.6,
        "EMP2.5": 100.1,
        "EMP10": 101.7,
        "EMP25": 91.1,
        "EMP100": 85.4,
    }  # [kg]
    doses: ClassVar[dict[str, float]] = {
        "Placebo": 0,
        "EMP2.5": 2.5,
        "EMP10": 10,
        "EMP25": 25,
        "EMP100": 100,
    }  # [mg]
    info: ClassVar[dict[str, str]] = {
        "[Cve_emp]": "empagliflozin",
        "KI__UGE": "uge",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        for intervention, dose in self.doses.items():
            # doses at 0 hr and every 24 hr from 48 to 192 hr (day 2 without
            # a dose), the urinary glucose is collected over 24 hr from every
            # dose (the data of days 1, 8 and 9) and reset after the last one
            simulations[self.to_sid(intervention)] = Simulation(
                time_unit="hr",
                end=292,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[intervention], "kg"),
                    "[KI__fpg]": Q(self.fpgs[intervention], "mM"),
                    "KI__PG_delta": Q(
                        self.mdgs[intervention] - self.fpgs[intervention], "mM"
                    ),
                    "KI__f_renal_function": Q(
                        self.crcls[intervention] / self.gfr_healthy, "dimensionless"
                    ),
                    "PODOSE_emp": Q(dose, "mg"),
                    "KI__glc_urine": Q(0, "mmole"),
                },
                changes=[
                    Change([24, 216], {"KI__glc_urine": Q(0, "mmole")}),
                    Change(
                        [48 + 24 * k for k in range(7)],
                        {"KI__glc_urine": Q(0, "mmole")},
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
                if name == "empagliflozin" and dose == 0:
                    continue

                mappings[self.to_sid(f"fm_{name}_{intervention}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{intervention}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if name == "uge" else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{intervention}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if name == "uge" else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.T2DM,
                        fasting=Fasting.FASTED,  # overnight fast of 10 hr
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (T2DM)",
        )
        # the legends inside the panels cover data
        fig.settings = {**fig.settings, "legend_position": "outside"}
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        # plasma (linear and semi-log as in the paper) and the glucose excretion
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(
            self.label_emp_plasma, unit=self.unit_emp, scale="log", min=1e-3
        )
        plots[2].set_xaxis(self.label_time, unit=self.unit_time, min=-0.05, max=225)
        plots[2].set_yaxis(self.label_uge, unit=self.unit_uge, max=180)
        panels = {"empagliflozin": plots[:2], "uge": plots[2:]}

        for intervention, dose in self.doses.items():
            for sid, name in self.info.items():
                dose_label = "Placebo" if dose == 0 else f"{dose} mg Emp"
                color = self.dose_colors[dose]
                for plot in panels[name]:
                    # simulation
                    plot.add_data(
                        task=self.to_sid(f"task_{intervention}"),
                        xid="time",
                        yid=sid,
                        label=dose_label,
                        color=color,
                    )

                    # data
                    if name == "empagliflozin" and dose == 0:
                        continue

                    plot.add_data(
                        dataset=self.to_sid(f"{name}_{intervention}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if name == "uge" else None,
                        count="count",
                        label=dose_label,
                        color=color,
                        linestyle="" if name == "uge" else "--",
                    )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Heise2013, output_dir=Heise2013.__name__)
