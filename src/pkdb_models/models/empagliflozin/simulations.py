"""The simulation experiments of the empagliflozin model and their groups.

They are the registry of `MODEL`, the experiments its commands simulate.
"""

from sbmlsim.experiment import SimulationExperiment

from pkdb_models.models.empagliflozin import MODEL
from pkdb_models.models.empagliflozin.experiments.misc import (
    DoseDependencyExperiment,
)
from pkdb_models.models.empagliflozin.experiments.scans.scan_parameter import (
    EmpagliflozinParameterScan,
)
from pkdb_models.models.empagliflozin.experiments.studies import (
    Ayoub2017,
    Brand2012,
    Chen2015a,
    Chen2020,
    ElDash2021,
    Friedrich2013,
    Hailat2022,
    Heise2013,
    Heise2013a,
    Heise2015,
    Jiang2023b,
    Kim2021,
    Kim2023,
    Laffel2018,
    Li2020,
    Macha2013d,
    Macha2013e,
    Macha2014,
    Macha2014b,
    Macha2014f,
    Macha2015a,
    Macha2015b,
    Sarashina2013,
    Sarashina2014,
    Seman2013,
    Zhao2015,
    vanderAartvanderBeek2020,
)

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {
    "studies": [
        Ayoub2017,
        Brand2012,
        Chen2015a,
        Chen2020,
        ElDash2021,
        Friedrich2013,
        Hailat2022,  # outlier (probably units)
        Heise2013,
        Heise2013a,
        Heise2015,
        Jiang2023b,
        Kim2021,
        Kim2023,
        Laffel2018,  # pediatric cohort
        Li2020,
        Macha2013d,  # the empagliflozin groups of the ramipril study Macha2013
        Macha2013e,
        Macha2014,
        Macha2014b,
        Macha2014f,
        Macha2015a,
        Macha2015b,
        Sarashina2013,
        Sarashina2014,
        Seman2013,
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
        Seman2013,
    ],
    "dose_dependency": [
        Heise2013,
        Heise2013a,
        Macha2015a,
        Sarashina2013,
        Seman2013,
        Zhao2015,
    ],
    "hepatic_impairment": [
        Macha2014b,
    ],
    "renal_impairment": [
        Macha2014f,
        Sarashina2014,
    ],
    "food": [
        Chen2020,
        Hailat2022,
        Li2020,
        Seman2013,  # Figure 4
    ],
    "misc": [
        DoseDependencyExperiment,
    ],
    "scans": [
        EmpagliflozinParameterScan,
    ],
}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] + EXPERIMENTS["scans"]


if __name__ == "__main__":
    MODEL.run_simulation_experiments(EXPERIMENTS["all"], MODEL.simulation_dir / "all")
