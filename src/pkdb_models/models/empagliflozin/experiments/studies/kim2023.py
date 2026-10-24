"""Simulation experiment of the study Kim2023."""

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


class Kim2023(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Kim2023."""

    # empagliflozin (EP25) alone and with evogliptin (EV5)
    colors: ClassVar[dict[str, str]] = {
        "EP25": "black",
        "EV5, EP25": "tab:blue",
    }
    interventions: ClassVar[list[str]] = list(colors)
    #: body weight of the empagliflozin arm (group 1), Table S1 [kg]
    bodyweight = 69.5

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "TabS2"]:
            df = self.load_dataframe(fig_id)
            if fig_id == "TabS2":
                # the 24 hr urine of day 5 (1 day after the 5th dose); the rows
                # at 0 days are not data
                df = df[df["time"] > 0]
            for label, df_label in df.groupby("label"):
                label = str(label)
                if not label.startswith(("empagliflozin_", "uge_")):
                    continue
                if label.startswith("uge_") and "EP25" not in label:
                    # the dapagliflozin arm
                    continue
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("empagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.emp)

                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        dose = Q(25, "mg")
        # a dose every 24 hr for 5 days, the time is relative to the last dose;
        # the urinary glucose is collected from dose to dose, the data are the
        # 24 hr after the last dose (the OGTT of day 4 is not simulated)
        simulations["po_emp25"] = Simulation(
            time_unit="hr",
            end=121,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg_healthy, "mM"),  # healthy reference value
                "KI__f_renal_function": Q(1.0, "dimensionless"),  # healthy
                "PODOSE_emp": dose,
            },
            changes=[
                Change([24, 48, 72, 96], {"KI__glc_urine": Q(0, "mmole")}),
                self.add_doses([24, 48, 72, 96], {"PODOSE_emp": dose}),
            ],
            time_shift=-96,
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        # plasma and the glucose excretion of day 5: dataset, observable and
        # the prefix of the key
        info = {
            "empagliflozin": ("[Cve_emp]", "fm_emp25"),
            "uge": ("KI__UGE", "fm_uge"),
        }
        for intervention in self.interventions:
            for name, (sid, prefix) in info.items():
                mappings[self.to_sid(f"{prefix}_{intervention}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{intervention}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task="task_po_emp25",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if name == "uge" else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.EVOGLIPTIN
                        if "EV5" in intervention
                        else Coadministration.NONE,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=3,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-5, max=25),
            legend=True,
        )
        # plasma (linear and semi-log as in the paper) and the glucose
        # excretion of day 5 (Table S2)
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(
            self.label_emp_plasma, unit=self.unit_emp, scale="log", min=0.02
        )
        plots[2].set_xaxis(self.label_time, unit=self.unit_time, min=-1, max=25)
        plots[2].set_yaxis(self.label_uge, unit=self.unit_uge)
        panels = {
            "empagliflozin": ("[Cve_emp]", plots[:2]),
            "uge": ("KI__UGE", plots[2:]),
        }

        label_map = {
            "EP25": "25 mg Emp",
            "EV5, EP25": "25 mg Emp + Evo",
        }

        for name, (sid, panel_plots) in panels.items():
            for plot in panel_plots:
                # simulation, the same for both arms (no interaction)
                plot.add_data(
                    task="task_po_emp25",
                    xid="time",
                    yid=sid,
                    label="25 mg Emp",
                    color="black",
                )

                # data
                for intervention in self.interventions:
                    plot.add_data(
                        dataset=self.to_sid(f"{name}_{intervention}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=label_map[intervention],
                        color=self.colors[intervention],
                        linestyle="" if name == "uge" else "--",
                    )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kim2023, output_dir=Kim2023.__name__)
