"""Simulation experiment of the study Heise2015."""

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


class Heise2015(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Heise2015."""

    #: fasting plasma glucose of T2DM patients on metformin [mM], not reported;
    #: it changes only the glucose excretion, which is not mapped
    fpg = 10
    #: body weight [kg], not reported: the mean BMI 27.1 kg/m2 at an assumed
    #: height of 1.69 m
    bodyweight = 77.8

    dset_ids: ClassVar[list[str]] = [
        "emp25",
        "emp25,hctz25",
        "emp25,tor5",
    ]
    colors: ClassVar[dict[str, str]] = {
        "emp25": "black",
        "emp25,hctz25": "tab:blue",
        "emp25,tor5": "tab:orange",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        dose = Q(25, "mg")
        # a dose every 24 hr for 5 days, the time is relative to the last dose
        simulations["po_emp25"] = Simulation(
            time_unit="hr",
            end=121,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg, "mM"),
                "PODOSE_emp": dose,
            },
            changes=[self.add_doses([24, 48, 72, 96], {"PODOSE_emp": dose})],
            time_shift=-96,
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for dset in self.dset_ids:
            coadministration = Coadministration.NONE
            if "hctz25" in dset:
                coadministration = Coadministration.HCTZ
            elif "tor5" in dset:
                coadministration = Coadministration.TORASEMIDE

            mappings[self.to_sid(f"fm_po_emp25_{dset}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(dset),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_emp25",
                    xid="time",
                    yid="[Cve_emp]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.T2DM,
                    fasting=Fasting.FASTED,  # overnight fast of 10 hr
                    coadministration=coadministration,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_cols=2,
            name=f"{self.__class__.__name__} (T2DM)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-5, max=25),
            legend=True,
        )
        # linear and semi-log as in the paper
        plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[1].set_yaxis(
            self.label_emp_plasma, unit=self.unit_emp, scale="log", min=0.02
        )
        coadmin_labels = {
            "emp25": None,
            "emp25,hctz25": "HCTZ",
            "emp25,tor5": "Tor",
        }

        for plot in plots:
            # simulation
            plot.add_data(
                task="task_po_emp25",
                xid="time",
                yid="[Cve_emp]",
                label="25 mg Emp",
                color="black",
            )

            for dset in self.dset_ids:
                coadmin = coadmin_labels[dset]
                label = "25 mg Emp" if coadmin is None else f"25 mg Emp + {coadmin}"

                plot.add_data(
                    dataset=self.to_sid(dset),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=self.colors[dset],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Heise2015, output_dir=Heise2015.__name__)
