"""Parameter fit problems for empagliflozin."""
from typing import Dict, List
from sbmlsim.fit.helpers import f_fitexp, filter_empty, filter_outlier
from sbmlutils.console import console
from sbmlutils.log import get_logger

from sbmlsim.fit import FitExperiment, FitMapping

from pkdb_models.models.empagliflozin import EMPAGLIFLOZIN_PATH, DATA_PATHS
from pkdb_models.models.empagliflozin.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, EmpagliflozinMappingMetaData, Coadministration
)
from pkdb_models.models.empagliflozin.experiments.studies import *


logger = get_logger(__name__)


# --- Experiment classes ---
experiment_classes = [
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
    Li2020,
    Macha2014b,
    Macha2014f,
    Macha2015a,
    Sarashina2013,
    Sarashina2014,
    Seman2013,
    Zhao2015,
]


# --- Filters ---
def filter_control(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Return control experiments/mappings."""

    metadata: EmpagliflozinMappingMetaData = fit_mapping.metadata

    # only PO and IV (no SL, MU, RE)
    if metadata.route not in {Route.PO, Route.IV}:
        return False

    # filter coadminstration
    if metadata.coadministration != Coadministration.NONE:
        return False

    # filter health (no renal, cardiac impairment, ...)
    if metadata.health not in {Health.HEALTHY, Health.T2DM, Health.HYPERTENSION}:
        return False

    # filter multiple dosing (only single dosing)
    # if metadata.dosing == Dosing.MULTIPLE:
    #     return False

    # only fasted subjects
    if metadata.fasting not in {Fasting.FASTED, Fasting.NR}:
        return False

    # remove outliers
    if metadata.outlier is True:
        return False

    return True


def filter_pk(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only pharmacokinetics data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if yid not in {
        "Cve_emp",
        "Cve_eg",
        "Cve_emptot",
        "Aurine_emp",
        "Aurine_eg",
        "Aurine_emptot",
        "Afeces_emp",
        "Afeces_eg",
        "Afeces_emptot",
    }:
        return False
    return True


def filter_pd(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only pharmacodynamics data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if yid not in {
        "KI__RTG",
        "KI__UGE"
        }:
        return False
    return True


# --- Fit experiments ---
def f_fitexp_all() -> Dict[str, List[FitExperiment]]:
    """Control data."""
    return f_fitexp(experiment_classes, metadata_filters=filter_empty,
        base_path=EMPAGLIFLOZIN_PATH,
        data_path=DATA_PATHS,
    )


def f_fitexp_control() -> Dict[str, List[FitExperiment]]:
    """Control data."""
    return f_fitexp(experiment_classes, metadata_filters=filter_control,
        base_path=EMPAGLIFLOZIN_PATH,
        data_path=DATA_PATHS,
    )


def f_fitexp_pk() -> Dict[str, List[FitExperiment]]:
    """Control data."""
    return f_fitexp(experiment_classes, metadata_filters=[filter_control, filter_pk],
        base_path=EMPAGLIFLOZIN_PATH,
        data_path=DATA_PATHS,
    )


def f_fitexp_pd() -> Dict[str, List[FitExperiment]]:
    """Control data."""
    return f_fitexp(experiment_classes, metadata_filters=[filter_control, filter_pd],
        base_path=EMPAGLIFLOZIN_PATH,
        data_path=DATA_PATHS,
    )


if __name__ == "__main__":
    """Test construction of FitExperiments."""

    for f in [
        f_fitexp_all,
        f_fitexp_control,
        # f_fitexp_pk,
        # f_fitexp_pd,
    ]:
        console.rule(style="white")
        console.print(f"{f.__name__}")
        fitexp = f()
