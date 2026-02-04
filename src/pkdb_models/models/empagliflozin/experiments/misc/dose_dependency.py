from typing import Dict

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)
from pkdb_models.models.empagliflozin.helpers import run_experiments


class DoseDependencyExperiment(EmpagliflozinSimulationExperiment):
    """Test dose dependency of Empagliflozin PO application."""

    doses = [0, 2.5, 5, 10, 50, 100, 200, 800]  # [mg]
    glucoses = [5, 6, 7, 8, 9, 10, 11]  # [mM]

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for dose in self.doses:
            tcsims[f"emp_dose_{dose}"] = TimecourseSim(
                Timecourse(
                    start=0,
                    end=25 * 60,  # [min]
                    steps=1000,
                    changes={
                        **self.default_changes(),
                        "PODOSE_emp": Q_(dose, "mg"),
                    },
                )
            )

        for glc in self.glucoses:
            tcsims[f"emp_glucose_{glc}"] = TimecourseSim(
                Timecourse(
                    start=0,
                    end=25 * 60,  # [min]
                    steps=5000,
                    changes={
                        **self.default_changes(),
                        "PODOSE_emp": Q_(25, "mg"),
                        "[KI__fpg]": Q_(glc, "mM"),  # FPG
                    },
                )
            )

        return tcsims

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def figure_pk(self) -> Dict[str, Figure]:
        figures = {}
        for key in ["dose", "glucose"]:
            values = self.doses if key == "dose" else self.glucoses

            fig = Figure(
                experiment=self,
                sid=f"Fig_{key}_dependency_pk",
                name="Dose dependency of empagliflozin (PK)" if key == "dose" else "Glucose dependency (PK)",
                num_rows=3,
                num_cols=3,
            )
            plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)

            sids = [
                # Plasma
                "[Cve_emp]",
                "[Cve_eg]",
                "[Cve_emptot]",
                # Urine
                "Aurine_emp",
                "Aurine_eg",
                "Aurine_emptot",
                # Feces
                "Afeces_emp",
                "Afeces_eg",
                "Afeces_emptot",
            ]

            for ksid, sid in enumerate(sids):
                if not sid:
                    continue
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for value in values:
                    color = self.dose_colors[value] if key == "dose" else self.glucose_colors[value]
                    label = (
                        f"{value} mg Emp PO"
                        if key == "dose"
                        else f"FPG {value} mM"
                    )
                    plots[ksid].add_data(
                        task=f"task_emp_{key}_{value}",
                        xid="time",
                        yid=sid,
                        label=label,
                        color=color,
                    )

            figures[fig.sid] = fig

        return figures

    def figure_pd(self) -> Dict[str, Figure]:
        figures = {}
        for key in ["dose", "glucose"]:
            values = self.doses if key == "dose" else self.glucoses

            fig = Figure(
                experiment=self,
                sid=f"Fig_{key}_dependency_pd",
                name="Dose dependency of empagliflozin (PD)" if key == "dose" else "Glucose dependency (PD)",
                num_rows=2,
                num_cols=3,
            )
            plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)

            sids = ["KI__UGE", "KI__RTG", "[KI__fpg]"]
            for ksid, sid in enumerate(sids):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            # time courses
            for ksid, sid in enumerate(sids):
                for value in values:
                    color = self.dose_colors[value] if key == "dose" else self.glucose_colors[value]
                    label = (
                        f"{value} mg Emp PO"
                        if key == "dose"
                        else f"FPG {value} mM"
                    )
                    plots[ksid].add_data(
                        task=f"task_emp_{key}_{value}",
                        xid="time",
                        yid=sid,
                        label=label,
                        color=color,
                    )

            # concentration-response panels
            for ksid in range(3, 6):
                plots[ksid].set_xaxis(label=self.label_emp, unit=self.unit_emp)

            plots[3].set_yaxis(label="Rate glucose excretion", unit="mmole/min")
            plots[4].set_yaxis(label=self.label_rtg, unit=self.unit_rtg)
            plots[5].set_yaxis(label=self.label_uge, unit=self.unit_uge)

            for value in values:
                color = self.dose_colors[value] if key == "dose" else self.glucose_colors[value]
                label = (
                    f"{value} mg Emp PO"
                    if key == "dose"
                    else f"FPG {value} mM"
                )

                plots[3].add_data(
                    task=f"task_emp_{key}_{value}",
                    xid="[Cve_emp]",
                    yid="KI__GLCEX",
                    label=label,
                    color=color,
                )
                plots[4].add_data(
                    task=f"task_emp_{key}_{value}",
                    xid="[Cve_emp]",
                    yid="KI__RTG",
                    label=label,
                    color=color,
                )
                plots[5].add_data(
                    task=f"task_emp_{key}_{value}",
                    xid="[Cve_emp]",
                    yid="KI__UGE",
                    label=label,
                    color=color,
                )

            figures[fig.sid] = fig

        return figures


if __name__ == "__main__":
    run_experiments(DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__)
