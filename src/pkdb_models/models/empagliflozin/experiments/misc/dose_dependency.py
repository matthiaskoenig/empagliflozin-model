"""Dependency of the pharmacokinetics and pharmacodynamics on the dose and glucose."""

from typing import ClassVar, override

from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.empagliflozin import MODEL
from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)


class DoseDependencyExperiment(EmpagliflozinSimulationExperiment):
    """Test dose dependency of Empagliflozin PO application."""

    doses: ClassVar[list[float]] = [0, 2.5, 5, 10, 50, 100, 200, 800]  # [mg]
    glucoses: ClassVar[list[int]] = [5, 6, 7, 8, 9, 10, 11]  # [mM]

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        for dose in self.doses:
            tcsims[self.to_sid(f"emp_dose_{dose}")] = Simulation(
                time_unit="min",
                end=25 * 60 * 7,
                preinit_changes={
                    **self.default_changes(),
                    "PODOSE_emp": Q(dose, "mg"),
                },
            )

        for glc in self.glucoses:
            tcsims[self.to_sid(f"emp_glucose_{glc}")] = Simulation(
                time_unit="min",
                end=25 * 60 * 7,
                preinit_changes={
                    **self.default_changes(),
                    "PODOSE_emp": Q(25, "mg"),
                    "[KI__fpg]": Q(glc, "mM"),  # FPG
                },
            )

        return tcsims

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def _curves(self, key: str) -> list[tuple[float, str, str]]:
        """Values of the doses or glucose levels with the color and label of the curves.

        Args:
            key: `dose` or `glucose`.
        """
        if key == "dose":
            return [
                (dose, self.dose_colors[dose], f"{dose} mg Emp PO")
                for dose in self.doses
            ]
        return [
            (glc, self.glucose_colors[glc], f"FPG {glc} mM") for glc in self.glucoses
        ]

    def figure_pk(self) -> dict[str, Figure]:
        """Timecourses of empagliflozin and its glucuronide in plasma, urine and feces."""
        figures = []
        for key in ["dose", "glucose"]:
            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_{key}_dependency_pk"),
                name="Dose dependency of empagliflozin (PK)"
                if key == "dose"
                else "Glucose dependency (PK) - 25 mg EMP PO",
                num_rows=3,
                num_cols=3,
                settings={"legend_scope": "row"},
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )

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
                "Afeces_emptot",
                # renal excretion rate
                "KI__EMPEX",
            ]

            for ksid, sid in enumerate(sids):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for value, color, label in self._curves(key):
                    plots[ksid].add_data(
                        task=self.to_sid(f"task_emp_{key}_{value}"),
                        xid="time",
                        yid=sid,
                        label=label,
                        color=color,
                    )

            figures.append(fig)

        return self.figures_by_sid(*figures)

    def figure_pd(self) -> dict[str, Figure]:
        """Timecourses of the glucose excretion and its dependency on empagliflozin."""
        figures = []
        for key in ["dose", "glucose"]:
            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_{key}_dependency_pd"),
                name="Dose dependency of empagliflozin (PD)"
                if key == "dose"
                else "Glucose dependency (PD) - 25 mg EMP PO",
                num_rows=2,
                num_cols=3,
                settings={"legend_scope": "row"},
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )

            sids = ["KI__UGE", "KI__RTG", "[KI__fpg]"]
            for ksid, sid in enumerate(sids):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            # time courses
            for ksid, sid in enumerate(sids):
                for value, color, label in self._curves(key):
                    plots[ksid].add_data(
                        task=self.to_sid(f"task_emp_{key}_{value}"),
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

            for value, color, label in self._curves(key):
                plots[3].add_data(
                    task=self.to_sid(f"task_emp_{key}_{value}"),
                    xid="[Cve_emp]",
                    yid="KI__GLCEX",
                    label=label,
                    color=color,
                )
                plots[4].add_data(
                    task=self.to_sid(f"task_emp_{key}_{value}"),
                    xid="[Cve_emp]",
                    yid="KI__RTG",
                    label=label,
                    color=color,
                )
                plots[5].add_data(
                    task=self.to_sid(f"task_emp_{key}_{value}"),
                    xid="[Cve_emp]",
                    yid="KI__UGE",
                    label=label,
                    color=color,
                )

            figures.append(fig)

        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(
        DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__
    )
