from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)
from pkdb_models.models.empagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, EmpagliflozinMappingMetaData, Coadministration

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.helpers import run_experiments


class Heise2013a(EmpagliflozinSimulationExperiment):
    """Simulation experiment of Heise2013a."""
    fpgs = {
        # using MDG values
        "Placebo": 152.1/18,  # 153.9/18,  # [mM]
        "EMP10": 164.8/18,  #  186.2/18,  # [mM]
        "EMP25": 166.3/18,  #  167.5/18,  # [mM]
        "EMP100": 149.4/18,  #  149.4/18,  # [mM]
    }
    bodyweights = {
        "Placebo": 89.3,  # [kg]
        "EMP10": 89.6,  # [kg]
        "EMP25": 89.6,  # [kg]
        "EMP100": 92.4,  # [kg]
    }
    doses = {
        "Placebo": 0,  # [mg]
        "EMP10": 10,  # [mg]
        "EMP25": 25,  # [mg]
        "EMP100": 100,  # [mg]
    }
    colors = {
        "Placebo": "black",
        "EMP10": "tab:red",
        "EMP25": "tab:orange",
        "EMP100": "tab:green",
    }
    info = {
        "[Cve_emp]": "empagliflozin",
        "KI__UGE": "uge",
    }
    interventions = list(bodyweights.keys())

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for intervention, dose in self.doses.items():
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweights[intervention], "kg"),
                    "[KI__fpg]": Q_(self.fpgs[intervention], "mM"),
                    # "KI__f_renal_function": Q_(self.gfr/100, "dimensionless"),
                    "PODOSE_emp": Q_(dose, "mg"),
                    "KI__glc_urine": Q_(0, "mmole"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_emp": Q_(dose, "mg"),
                    "KI__glc_urine": Q_(0, "mmole"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=100 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_emp": Q_(dose, "mg"),
                    "KI__glc_urine": Q_(0, "mmole"),
                },
            )
            tcsims[intervention] = TimecourseSim(
                [tc0] + [tc1 for _ in range(26)] + [tc2],
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for intervention, dose in self.doses.items():
            for k, sid in enumerate(self.info):
                name = self.info[sid]
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
                        self, task=f"task_{intervention}", xid="time", yid=sid,
                    ),
                    metadata=EmpagliflozinMappingMetaData(
                        tissue=Tissue.URINE if name == "uge" else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.T2DM,
                        fasting=Fasting.NR,
                        coadministration=Coadministration.NONE,
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        figures = {}
        subplots = ["all", "start", "end"]
        for subplot in subplots:
            fig = Figure(
                experiment=self,
                sid=f"Fig1_2_{subplot}",
                num_rows=1,
                num_cols=2,
                name=f"{self.__class__.__name__} (T2DM)",
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
            plots[0].set_yaxis(self.label_emp_plasma, unit=self.unit_emp)
            plots[1].set_yaxis(self.label_uge, unit=self.unit_uge)

            if subplot == "start":
                for kp in [0, 1]:
                    plots[kp].xaxis.min = -3
                    plots[kp].xaxis.max = 30
            elif subplot == "end":
                for kp in [0, 1]:
                    plots[kp].xaxis.min = 640
                    plots[kp].xaxis.max = 750

            for intervention, dose in self.doses.items():
                dose_label = "Placebo" if dose == 0 else f"{dose} mg Emp"
                color = self.dose_colors[dose]

                for k, sid in enumerate(self.info):
                    name = self.info[sid]

                    # simulation
                    plots[k].add_data(
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                        label=dose_label,
                        color=color,
                    )

                    # data
                    if name == "empagliflozin" and dose == 0:
                        continue

                    plots[k].add_data(
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=dose_label,
                        color=color,
                        linestyle="" if "uge" in name else "--",
                    )

            figures[fig.sid] = fig

        return figures


if __name__ == "__main__":
    run_experiments(Heise2013a, output_dir=Heise2013a.__name__)
