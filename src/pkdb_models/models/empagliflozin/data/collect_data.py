"""Collect the data of the empagliflozin studies from the PK-DB study folders.

The studies are copied from a `pkdb_data` checkout (`PKDB_DATA`, else
`../pkdb_data`, see `pkpdsim.data.studies_dir`) into the data folder of the
model, see `pkpdsim.data.collect_tsv_files`.
Besides the empagliflozin studies, studies of other substances have
empagliflozin groups: Kim2023 and vanderAartvanderBeek2020 (dapagliflozin),
Heise2015 (hydrochlorothiazide) and Macha2013 (ramipril), the data of the
experiment Macha2013d.
"""

from pathlib import Path

from pkpdsim.data import collect_tsv_files
from pkpdsim.data import studies_dir as pkdb_studies_dir

from pkdb_models.models.empagliflozin import (
    DATA_PATH_DAPAGLIFLOZIN,
    DATA_PATH_EMPAGLIFLOZIN,
    DATA_PATH_HCTZ,
    DATA_PATH_RAMIPRIL,
)

#: the studies of other substances which the empagliflozin experiments use
DAPAGLIFLOZIN_STUDIES = {"Kim2023", "vanderAartvanderBeek2020"}
HYDROCHLOROTHIAZIDE_STUDIES = {"Heise2015"}
RAMIPRIL_STUDIES = {"Macha2013"}


def collect_data(studies_dir: Path | None = None) -> None:
    """Copy the data of the empagliflozin studies into the data folder of the model.

    Args:
        studies_dir: the `studies` folder of pkdb_data, by default the one of
            the checkout `pkpdsim.data.studies_dir` finds.
    """
    if studies_dir is None:
        studies_dir = pkdb_studies_dir()
    collect_tsv_files(
        source_dir=studies_dir / "empagliflozin", target_dir=DATA_PATH_EMPAGLIFLOZIN
    )
    collect_tsv_files(
        source_dir=studies_dir / "dapagliflozin",
        target_dir=DATA_PATH_DAPAGLIFLOZIN,
        filter_study=lambda study: study in DAPAGLIFLOZIN_STUDIES,
    )
    collect_tsv_files(
        source_dir=studies_dir / "hydrochlorothiazide",
        target_dir=DATA_PATH_HCTZ,
        filter_study=lambda study: study in HYDROCHLOROTHIAZIDE_STUDIES,
    )
    collect_tsv_files(
        source_dir=studies_dir / "ramipril",
        target_dir=DATA_PATH_RAMIPRIL,
        filter_study=lambda study: study in RAMIPRIL_STUDIES,
    )


if __name__ == "__main__":
    collect_data()
