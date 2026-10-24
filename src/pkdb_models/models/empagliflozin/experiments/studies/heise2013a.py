"""Simulation experiment of the study Heise2013a."""

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


class Heise2013a(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Heise2013a."""

    #: fasting plasma glucose (FPG) at baseline, Table 1 [mM]
    fpgs: ClassVar[dict[str, float]] = {
        "Placebo": 153.9 / 18,
        "EMP10": 186.2 / 18,
        "EMP25": 167.5 / 18,
        "EMP100": 150.0 / 18,
    }
    #: mean daily glucose (MDG) at baseline, Table 1 [mM]; below the FPG in
    #: the 10 mg group
    mdgs: ClassVar[dict[str, float]] = {
        "Placebo": 152.1 / 18,
        "EMP10": 164.8 / 18,
        "EMP25": 166.3 / 18,
        "EMP100": 149.4 / 18,
    }
    #: creatinine clearance (Cockcroft-Gault) at baseline, Table 1 [ml/min];
    #: with the total body weight (BMI 29.5-30.6) it overestimates the GFR
    crcls: ClassVar[dict[str, float]] = {
        "Placebo": 133.3,
        "EMP10": 114.0,
        "EMP25": 130.2,
        "EMP100": 122.6,
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "Placebo": 89.3,  # [kg]
        "EMP10": 89.6,  # [kg]
        "EMP25": 89.6,  # [kg]
        "EMP100": 92.4,  # [kg]
    }
    doses: ClassVar[dict[str, float]] = {
        "Placebo": 0,  # [mg]
        "EMP10": 10,  # [mg]
        "EMP25": 25,  # [mg]
        "EMP100": 100,  # [mg]
    }
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
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        for intervention, dose in self.doses.items():
            # a dose every 24 hr for 28 days, the urinary glucose is collected
            # from dose to dose and reset after the last one
            simulations[intervention] = Simulation(
                time_unit="hr",
                end=748,
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
                    Change([672], {"KI__glc_urine": Q(0, "mmole")}),
                    Change(
                        [24 * k for k in range(1, 28)],
                        {"KI__glc_urine": Q(0, "mmole")},
                    ),
                    self.add_doses(
                        [24 * k for k in range(1, 28)], {"PODOSE_emp": Q(dose, "mg")}
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

                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if name == "uge" else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if name == "uge" else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.T2DM,
                        fasting=Fasting.FASTED,  # overnight fast
                        coadministration=Coadministration.NONE,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        # time window of the figures: x limits of the plasma (linear and
        # semi-log as in the paper) and of the urine plot, lower limit of the
        # linear plasma concentration
        xlimits = {
            "all": [(-10, 750), (-10, 750), (-10, 680)],
            "start": [(-1, 27), (-1, 27), (-1, 27)],
            "end": [(645, 750), (645, 750), (645, 675)],
        }
        ymins = {"all": -0.2, "start": -0.1, "end": -0.1}

        figures = []
        for subplot, xlims in xlimits.items():
            fig = Figure(
                experiment=self,
                sid=f"Fig1_2_{subplot}",
                num_rows=1,
                num_cols=3,
                name=f"{self.__class__.__name__} (T2DM)",
            )
            # the legends inside the panels cover data
            fig.settings = {**fig.settings, "legend_position": "outside"}
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )
            for kp, (xmin, xmax) in enumerate(xlims):
                plots[kp].set_xaxis(
                    self.label_time, unit=self.unit_time, min=xmin, max=xmax
                )
            plots[0].set_yaxis(
                self.label_emp_plasma, unit=self.unit_emp, min=ymins[subplot]
            )
            plots[1].set_yaxis(
                self.label_emp_plasma, unit=self.unit_emp, scale="log", min=1e-3
            )
            plots[2].set_yaxis(self.label_uge, unit=self.unit_uge, max=150)
            panels = {"empagliflozin": plots[:2], "uge": plots[2:]}

            for intervention, dose in self.doses.items():
                dose_label = "Placebo" if dose == 0 else f"{dose} mg Emp"
                color = self.dose_colors[dose]

                for sid, name in self.info.items():
                    for plot in panels[name]:
                        # simulation
                        plot.add_data(
                            task=f"task_{intervention}",
                            xid="time",
                            yid=sid,
                            label=dose_label,
                            color=color,
                        )

                        # data
                        if name == "empagliflozin" and dose == 0:
                            continue

                        # the plasma data have a SD only at some points
                        plot.add_data(
                            dataset=f"{name}_{intervention}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd" if name == "uge" else None,
                            count="count",
                            label=dose_label,
                            color=color,
                            linestyle="" if name == "uge" else "--",
                        )

            figures.append(fig)

        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(Heise2013a, output_dir=Heise2013a.__name__)
