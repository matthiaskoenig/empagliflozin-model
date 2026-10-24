"""Simulation experiment of the study Chen2015a."""

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


class Chen2015a(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Chen2015a.

    Mass balance of eight healthy men (79 kg): a single oral solution of 50 mg
    empagliflozin with about 100 µCi [14C]-empagliflozin after an overnight
    fast of 10 hr (no food until 4 hr after the dose). Plasma (Fig2), urine
    and feces (Fig1) of empagliflozin (LC-MS/MS) and of the total
    radioactivity, the metabolites by radiochromatography (Table 2).

    The glucuronide in plasma is the sum of the glucuronides M626/1, M626/2 and
    M626/3 of Table 2 (2, 6 and 12 hr). The curve "EG_plasma" of the curation
    is the plasma radioactivity minus the empagliflozin of the LC-MS/MS assay
    (zero up to 1 hr, where the assays disagree): it underestimates the
    glucuronides of Table 2 about twofold. The glucuronide in urine is the
    urinary radioactivity minus empagliflozin (26.5 % of the dose at 96 hr, the
    glucuronides of Table 2 are 23.1 %, the rest are oxidative metabolites).
    No glucuronide is found in the feces: the fecal curves are the fecal
    radioactivity times the fraction of empagliflozin (82.9 %) of the pooled
    feces, "EG_feces" (the remaining 17.1 %, oxidative metabolites and
    unidentified radioactivity) is not mapped.
    """

    fpg = EmpagliflozinSimulationExperiment.fpg_healthy  # [mM] (healthy subjects)
    bodyweight = 79  # [kg]

    #: the glucuronides of empagliflozin in Table 2
    glucuronides: ClassVar[list[str]] = ["C14-M626/1", "C14-M626/2", "C14-M626/3"]

    # datasets are named "<species>_<tissue>"
    info: ClassVar[dict[str, str]] = {
        # Fig2
        "[Cve_emptot]": "EMPTOT_plasma",
        "[Cve_emp]": "EMP_plasma",
        "[Cve_eg]": "EG_plasma",
        # Fig1
        "Afeces_emptot": "EMPTOT_feces",
        "Afeces_emp": "EMP_feces",
        "Aurine_emptot": "EMPTOT_urine",
        "Aurine_emp": "EMP_urine",
        "Aurine_eg": "EG_urine",
    }
    colors: ClassVar[dict[str, str]] = {
        "[Cve_emptot]": "black",
        "[Cve_emp]": "tab:blue",
        "[Cve_eg]": "tab:orange",
        # Fig1
        "Afeces_emptot": "black",
        "Afeces_emp": "tab:blue",
        "Aurine_emptot": "black",
        "Aurine_emp": "tab:blue",
        "Aurine_eg": "tab:orange",
    }
    tissues: ClassVar[dict[str, Tissue]] = {
        "plasma": Tissue.PLASMA,
        "urine": Tissue.URINE,
        "feces": Tissue.FECES,
    }
    # panel of the figure by tissue
    panels: ClassVar[dict[str, int]] = {"plasma": 0, "urine": 1, "feces": 2}

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[label] = dset

        # the glucuronides in plasma, the sum of Table 2
        df = self.load_dataframe("Tab2")
        df = df[df.substance.isin(self.glucuronides)]
        df_eg = df[df.substance == self.glucuronides[0]].copy()
        df_eg["mean"] = df.groupby("time")["mean"].sum().loc[df_eg.time].to_numpy()
        df_eg["substance"] = "empagliflozin-glucuronide"
        df_eg["label"] = "EG_plasma_Tab2"
        dsets["EG_plasma_Tab2"] = DataSet.from_df(df_eg, self.ureg)

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        simulations["po_C14EMP"] = Simulation(
            time_unit="hr",
            end=180,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg, "mM"),
                # 49.1 mg empagliflozin plus 0.9 mg [14C]-empagliflozin
                "PODOSE_emp": Q(50, "mg"),
            },
        )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        datasets = {name: sid for sid, name in self.info.items()}
        # the glucuronides of Table 2
        datasets["EG_plasma_Tab2"] = "[Cve_eg]"
        for name, sid in datasets.items():
            tissue = name.split("_")[1]
            mappings[f"fm_po_emp_{name}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=name,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_C14EMP",
                    xid="time",
                    yid=sid,
                ),
                metadata=PKPDMappingMetaData(
                    tissue=self.tissues[tissue],
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
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
            num_cols=3,  # plasma | urine | feces
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_xaxis(self.label_time, unit=self.unit_time, min=-0.05, max=25)
        plots[1].set_xaxis(self.label_time, unit=self.unit_time, min=-0.05, max=150)
        plots[0].set_yaxis(label="Plasma", unit=self.unit_emp, min=-0.05)
        plots[1].set_yaxis(label="Urine", unit=self.unit_emp_urine)
        plots[2].set_yaxis(label="Feces", unit=self.unit_emp_feces)

        # the glucuronides of Table 2
        plots[0].add_data(
            dataset="EG_plasma_Tab2",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="50 mg Emp (EG, Tab2)",
            color=self.colors["[Cve_eg]"],
            marker="s",
        )

        for sid, name in self.info.items():
            species, tissue = name.split("_")
            kp = self.panels[tissue]
            label = f"50 mg Emp ({species})"

            # simulation
            plots[kp].add_data(
                task="task_po_C14EMP",
                xid="time",
                yid=sid,
                label=label,
                color=self.colors[sid],
            )

            # data
            plots[kp].add_data(
                dataset=name,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                # the difference of radioactivity and empagliflozin, an outlier
                label="50 mg Emp (EG, difference)" if name == "EG_plasma" else label,
                color=self.colors[sid],
                marker="x" if name == "EG_plasma" else "o",
            )

        fig.settings = {**fig.settings, "legend_position": "outside"}
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Chen2015a, output_dir=Chen2015a.__name__)
