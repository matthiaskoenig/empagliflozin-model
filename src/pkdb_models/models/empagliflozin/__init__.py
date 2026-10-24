"""Empagliflozin: PBPK/PD model of empagliflozin and the urinary glucose excretion.

The paths of the model and its model object `MODEL`, which the commands
`run_simulation`, `run_fit`, `run_report` and `run_identifiability` of pkpdsim
run by its package `pkdb_models.models.empagliflozin`.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings
from pkpdsim.release.settings import Author, ReleaseSettings

from pkdb_models.models.empagliflozin.tolerances import ABSOLUTE_TOLERANCE

EMPAGLIFLOZIN_PATH = Path(__file__).parent

MODEL_BASE_PATH = EMPAGLIFLOZIN_PATH / "models" / "output" / "models"
MODEL_PATH = MODEL_BASE_PATH / "empagliflozin_body_flat.xml"

# the subset of the PK-DB studies the model uses, see `data/collect_data.py`
DATA_PATH_BASE = EMPAGLIFLOZIN_PATH / "data"
DATA_PATH_EMPAGLIFLOZIN = DATA_PATH_BASE / "empagliflozin"
DATA_PATH_DAPAGLIFLOZIN = DATA_PATH_BASE / "dapagliflozin"
DATA_PATH_HCTZ = DATA_PATH_BASE / "hydrochlorothiazide"
DATA_PATH_RAMIPRIL = DATA_PATH_BASE / "ramipril"
DATA_PATHS = [
    DATA_PATH_EMPAGLIFLOZIN,
    DATA_PATH_DAPAGLIFLOZIN,
    DATA_PATH_HCTZ,
    DATA_PATH_RAMIPRIL,
]

MODEL = PKPDModel(
    name="empagliflozin",
    package=__name__,
    title="Empagliflozin PBPK/PD model",
    path=EMPAGLIFLOZIN_PATH,
    model_path=MODEL_PATH,
    data_paths=tuple(DATA_PATHS),
    archive_path=MODEL_BASE_PATH.parent / "empagliflozin_model.omex",
    factory=".models.factory:create_empagliflozin_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
    simulator=SimulatorSettings(
        absolute_tolerance=ABSOLUTE_TOLERANCE, relative_tolerance=1e-10
    ),
    figures=FigureSettings(formats=("svg",), legend_fontsize=9, axes_labelsize=16),
    fits=("PK_20261019_201013__13ddc", "PD_20261019_212349__c28ac"),
    fit_identifiability=True,
    # the release of the model repository
    release=ReleaseSettings(
        title=(
            "Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model"
            " of empagliflozin"
        ),
        repository="matthiaskoenig/empagliflozin-model",
        doi="10.5281/zenodo.18480677",
        authors=(
            Author(
                "Alejandro, Jennesse",
                "0009-0009-3917-4482",
                "Brown University, Providence, Rhode Island",
            ),
            Author(
                "Elias, Michelle",
                "0009-0005-7959-3164",
                "Humboldt-University Berlin, Institute for Theoretical Biology, Berlin",
            ),
            Author(
                "Babaeva, Mariia",
                "0009-0007-6659-7638",
                "Charité – Universitätsmedizin Berlin, Charité Medical Faculty, Berlin",
            ),
            Author(
                "König, Matthias",
                "0000-0003-1725-179X",
                (
                    "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver; University of Stuttgart, Institute of Structural Mechanics and Dynamics in Aerospace Engineering; University of Southern Denmark, Clinical Pharmacology, Pharmacy and Environmental Medicine, Department of Public Health; University Hospital Schleswig-Holstein, Campus Lübeck, First Department of Medicine, Lübeck, Germany"
                ),
            ),
        ),
        keywords=(
            "model",
            "empagliflozin",
            "SGLT2 inhibitor",
            "urinary glucose excretion",
            "type 2 diabetes",
            "pharmacokinetics",
            "pharmacodynamics",
            "PBPK/PD",
            "physiologically based pharmacokinetic/ pharmacodynamic model",
            "SBML",
        ),
        docker_image="matthiaskoenig/empagliflozin",
        funding=(
            "Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B)."
        ),
        description=(
            "This repository provides the empagliflozin physiologically based"
            " pharmacokinetic/pharmacodynamic (PBPK/PD) model: the absorption,"
            " the tissue distribution, the hepatic and renal glucuronidation and"
            " the renal excretion of empagliflozin, and the inhibition of SGLT2,"
            " which lowers the renal threshold of glucose and increases the"
            " urinary glucose excretion."
        ),
        first_year=2024,
    ),
)
"""The empagliflozin model: its commands, see `pkpdsim.model`."""
