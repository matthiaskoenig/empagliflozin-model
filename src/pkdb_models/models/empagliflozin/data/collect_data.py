from pathlib import Path
from pkdb_models.models.data import collect_tsv_files

def collect_empagliflozin_data():
    common_parent: Path = Path(__file__).parents[5]
    source_dir = common_parent / "pkdb_data" / "studies" / "empagliflozin"
    target_dir = Path(__file__).parent / "empagliflozin"
    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)

    # Collect Kim2023 (dapagliflozin)
    def is_Kim2023(study_name) -> bool:
        return study_name == "Kim2023"
    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "dapagliflozin",
        target_dir=Path(__file__).parent / "dapagliflozin",
        filter_study=is_Kim2023,
    )

    # Collect Heise2015 (hydrochlorothiazide)
    def is_Heise2015(study_name) -> bool:
        return study_name == "Heise2015"
    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "hydrochlorothiazide",
        target_dir=Path(__file__).parent / "hydrochlorothiazide",
        filter_study=is_Heise2015,
    )

    # Collect vanderAartvanderBeek2020 (dapagliflozin)
    def is_vanderAartvanderBeek2020(study_name) -> bool:
        return study_name == "vanderAartvanderBeek2020"
    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "dapagliflozin",
        target_dir=Path(__file__).parent / "dapagliflozin",
        filter_study=is_vanderAartvanderBeek2020,
    )

if __name__ == "__main__":
    collect_empagliflozin_data()

