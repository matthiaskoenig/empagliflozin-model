"""Run all simulation experiments."""
import shutil
from typing import List
from pathlib import Path
from sbmlutils.console import console

from pkdb_models.models.empagliflozin.experiments.scans.scan_parameter import EmpagliflozinParameterScan
from pkdb_models.models.empagliflozin.helpers import run_experiments
from pkdb_models.models.empagliflozin.experiments.studies import *
from pkdb_models.models.empagliflozin.experiments.misc import *
import pkdb_models.models.empagliflozin as empagliflozin
from sbmlutils import log
from sbmlsim.plot import Figure

logger = log.get_logger(__name__)

EXPERIMENTS = {
    "studies": [
        Ayoub2017,
        Brand2012,
        Chen2015a,
        Chen2020,
        ElDash2021,
        Friedrich2013,
        # Hailat2022,  # Outlier
        Heise2013,
        Heise2013a,
        Heise2015,
        Jiang2023b,
        Kim2021,
        Kim2023,
        Li2020,
        # Macha2013d, # wip
        # Macha2013e, # wip
        # Macha2014 # wip
        Macha2014b,
        Macha2014f,
        Macha2015a,
        # Macha2015b, # wip
        Sarashina2013,
        Sarashina2014,
        Seman2013, # wip urine data
        vanderAartvanderBeek2020,
        Zhao2015,
    ],
    "pharmacodynamics": [
        ElDash2021,
        Heise2013,
        Heise2013a,
        Macha2014b,
        Macha2015a,
        Sarashina2013,
        Seman2013
    ],
    "dose_dependency": [
        Heise2013,
        Heise2013a,
        Macha2015a,
        Sarashina2013,
        Seman2013,
        Zhao2015
    ],
    "hepatic_impairment": [
        Macha2014b,
    ],
    "renal_impairment": [
        Macha2014f,
        Sarashina2014,
    ],
    "misc": [
        DoseDependencyExperiment,
    ],
    "scan": [
        EmpagliflozinParameterScan
    ]
}

EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] # + EXPERIMENTS["scan"]


def run_simulation_experiments(
        selected: str = None,
        experiment_classes: List = None,
        output_dir: Path = None
) -> None:
    """Run empagliflozin simulation experiments."""

    Figure.fig_dpi = 300
    Figure.legend_fontsize = 10

    # Determine which experiments to run
    if experiment_classes is not None:
        experiments_to_run = experiment_classes
        if output_dir is None:
            output_dir = empagliflozin.RESULTS_PATH_SIMULATION / "custom_selection"
    elif selected:
        # Using the 'selected' parameter
        if selected not in EXPERIMENTS:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(EXPERIMENTS.keys())}[/red]"
            )
            console.rule(style="red bold")
            return
        experiments_to_run = EXPERIMENTS[selected]
        if output_dir is None:
            output_dir = empagliflozin.RESULTS_PATH_SIMULATION / selected
    else:
        console.print("\n[red bold]Error: No experiments specified![/red bold]")
        console.print("[yellow]Use selected='all' or selected='studies' or provide experiment_classes=[...][/yellow]\n")
        return

    # Run the experiments
    run_experiments(experiment_classes=experiments_to_run, output_dir=output_dir)

    # Collect figures into one folder
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        if f.parent == figures_dir:
            continue
        try:
            shutil.copy2(f, figures_dir / f.name)
        except Exception as err:
            print(f"file {f.name} in {f.parent} fails, skipping. Error: {err}")
    console.print(f"Figures copied to: file://{figures_dir}", style="info")


if __name__ == "__main__":
    """
    # Run experiments

    # selected = "all"
    # selected = "misc"
    # selected = "studies"
    # selected = "pharmacodynamics"
    # selected = "dose_dependency"
    # selected = "hepatic_impairment"
    # selected = "renal_impairment"
    # selected = "scan"
    """

    run_simulation_experiments(selected="all")
