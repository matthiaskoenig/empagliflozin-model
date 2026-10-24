"""Simulation experiment of the study Macha2014."""

from typing import ClassVar, NamedTuple, override

import pandas as pd
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


class Arm(NamedTuple):
    """An arm of the study."""

    #: dose of empagliflozin [mg]
    dose: int
    #: the coadministered drug, a key of `interactions` of the base, or None
    drug: str | None
    #: the fasting state of the study, a key of `fasting_map`
    fasting: str
    #: the simulation of the arm
    simulation: str


class Macha2014(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Macha2014.

    Two open-label crossover studies of single oral doses of empagliflozin in
    healthy subjects (Figure 2, Tables I and III):

    - gemfibrozil study (n=18): 25 mg alone and 25 mg on day 3 of gemfibrozil
      600 mg BID (days 1-5), both 0.5 hr after a standardized breakfast which
      followed an overnight fast (Figure 2A);
    - rifampicin/probenecid study (n=18): 10 mg alone, 10 mg with a single dose
      of rifampicin 600 mg and 10 mg on day 2 of probenecid 500 mg BID (days
      1-4, n=16); fasting not reported. The 10 mg alone arm is drawn twice, in
      Figure 2B (with rifampicin, `EMP10_R`) and 2C (with probenecid,
      `EMP10_P`). Urinary excretion 0-72 hr of the arms alone and with
      probenecid (Table III).

    The interactions are simulated with the factors of the base experiment
    (`interaction_changes`): gemfibrozil inhibits OATP1B1 and OAT3, rifampicin
    OATP1B1/1B3 and probenecid OAT3 (observed AUC ratios 1.59, 1.35, 1.53).
    """

    #: the arms by the label of the data
    arms: ClassVar[dict[str, Arm]] = {
        "EMP25": Arm(25, None, "fed", "po_emp25"),
        "EMP25_GEM": Arm(25, "gemfibrozil", "fed", "po_emp25_gem"),
        "EMP10_R": Arm(10, None, "NR", "po_emp10"),
        "EMP10_P": Arm(10, None, "NR", "po_emp10"),
        "EMP10_RIF": Arm(10, "rifampicin", "NR", "po_emp10_rif"),
        "EMP10_PROB": Arm(10, "probenecid", "NR", "po_emp10_prob"),
    }
    #: the arms with the urinary excretion of Table III, by the intervention of
    #: the table
    urine_arms: ClassVar[dict[str, str]] = {
        "EMP10": "EMP10_R",
        "EMP10, PROB": "EMP10_PROB",
    }

    coadministrations: ClassVar[dict[str | None, Coadministration]] = {
        None: Coadministration.NONE,
        "gemfibrozil": Coadministration.GEMFIBROZIL,
        "rifampicin": Coadministration.RIFAMPICIN,
        "probenecid": Coadministration.PROBENECID,
    }
    fastings: ClassVar[dict[str, Fasting]] = {
        "fed": Fasting.FED,
        "NR": Fasting.NR,
    }
    colors: ClassVar[dict[str | None, str]] = {
        None: "black",
        "gemfibrozil": "tab:blue",
        "rifampicin": "tab:blue",
        "probenecid": "tab:orange",
    }
    arm_labels: ClassVar[dict[str, str]] = {
        "EMP25": "25 mg Emp",
        "EMP25_GEM": "25 mg Emp + Gem",
        "EMP10_R": "10 mg Emp",
        "EMP10_P": "10 mg Emp",
        "EMP10_RIF": "10 mg Emp + Rif",
        "EMP10_PROB": "10 mg Emp + Prob",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig2")
        for label, df_label in df.groupby("label"):
            label = str(label)
            dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)

        # the fraction of the dose excreted unchanged in the urine 0-72 hr
        # (Table III, mean (CV)) as the amount of empagliflozin
        df = self.load_dataframe("Tab3")
        df = df[df.measurement_type == "recovery"]
        for intervention, df_int in df.groupby("intervention"):
            arm = self.urine_arms[str(intervention)]
            dose = self.arms[arm].dose
            fe = df_int.iloc[0]
            df_urine = pd.DataFrame(
                {
                    "time": [fe.time],
                    "time_unit": [fe.time_unit],
                    "count": [fe["count"]],
                    "mean": [fe["mean"] / 100 * dose],
                    "sd": [fe["cv"] / 100 * fe["mean"] / 100 * dose],
                    "unit": ["mg"],
                }
            )
            dset = DataSet.from_df(df_urine, self.ureg)
            dset.unit_conversion("mean", 1 / self.Mr.emp)
            dsets[self.to_sid(f"empagliflozin_urine_{arm}")] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for arm in self.arms.values():
            if arm.simulation in simulations:
                continue
            simulations[arm.simulation] = Simulation(
                time_unit="hr",
                end=73,
                preinit_changes={
                    **self.default_changes(),
                    # weights not reported (BMI 23.5 and 23.0 kg/m2)
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),  # healthy reference value
                    "KI__f_renal_function": Q(1.0, "dimensionless"),  # healthy
                    "GU__fed": Q(self.fasting_map[arm.fasting], "dimensionless"),
                    **(
                        self.interaction_changes(arm.drug)
                        if arm.drug is not None
                        else {}
                    ),
                    "PODOSE_emp": Q(arm.dose, "mg"),
                },
            )
        return simulations

    def _metadata(self, arm: Arm, tissue: Tissue) -> PKPDMappingMetaData:
        """Get the metadata of the mappings of an arm.

        Args:
            arm: the arm.
            tissue: the tissue of the data.

        Returns:
            The metadata of the mapping.
        """
        return PKPDMappingMetaData(
            tissue=tissue,
            route=Route.PO,
            application_form=ApplicationForm.TABLET,
            dosing=Dosing.SINGLE,
            health=Health.HEALTHY,
            fasting=self.fastings[arm.fasting],
            coadministration=self.coadministrations[arm.drug],
        )

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for label, arm in self.arms.items():
            mappings[self.to_sid(f"fm_po_{label}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"empagliflozin_{label}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_{arm.simulation}",
                    xid="time",
                    yid="[Cve_emp]",
                ),
                metadata=self._metadata(arm, Tissue.PLASMA),
            )
        for label in self.urine_arms.values():
            arm = self.arms[label]
            mappings[self.to_sid(f"fm_po_{label}_urine")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"empagliflozin_urine_{label}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_{arm.simulation}",
                    xid="time",
                    yid="Aurine_emp",
                ),
                metadata=self._metadata(arm, Tissue.URINE),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__} (Healthy)",
            num_cols=4,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-1, max=75),
            legend=True,
        )
        for plot in plots[:3]:
            plot.set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
        plots[3].set_yaxis(self.label_emp_urine, unit=self.unit_emp_urine)

        # A: gemfibrozil, B: rifampicin, C: probenecid, D: urine
        panels = {
            0: ["EMP25", "EMP25_GEM"],
            1: ["EMP10_R", "EMP10_RIF"],
            2: ["EMP10_P", "EMP10_PROB"],
        }
        for kp, labels in panels.items():
            for label in labels:
                arm = self.arms[label]
                plots[kp].add_data(
                    task=f"task_{arm.simulation}",
                    xid="time",
                    yid="[Cve_emp]",
                    label=self.arm_labels[label],
                    color=self.colors[arm.drug],
                )
                plots[kp].add_data(
                    dataset=self.to_sid(f"empagliflozin_{label}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.arm_labels[label],
                    color=self.colors[arm.drug],
                )
        for label in self.urine_arms.values():
            arm = self.arms[label]
            plots[3].add_data(
                task=f"task_{arm.simulation}",
                xid="time",
                yid="Aurine_emp",
                label=self.arm_labels[label],
                color=self.colors[arm.drug],
            )
            plots[3].add_data(
                dataset=self.to_sid(f"empagliflozin_urine_{label}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.arm_labels[label],
                color=self.colors[arm.drug],
                linestyle="",
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Macha2014, output_dir=Macha2014.__name__)
