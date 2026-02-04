from pathlib import Path

EMPAGLIFLOZIN_PATH = Path(__file__).parent

MODEL_BASE_PATH = EMPAGLIFLOZIN_PATH / "models" / "results" / "models"
MODEL_PATH = MODEL_BASE_PATH / "empagliflozin_body_flat.xml"

RESULTS_PATH = EMPAGLIFLOZIN_PATH / "results"
RESULTS_PATH_SIMULATION = RESULTS_PATH / "simulation"
RESULTS_PATH_FIT = RESULTS_PATH / "fit"

# DATA_PATH_BASE = EMPAGLIFLOZIN_PATH.parents[3] / "pkdb_data" / "studies"
DATA_PATH_BASE = EMPAGLIFLOZIN_PATH / "data"

DATA_PATH_EMPAGLIFLOZIN = DATA_PATH_BASE / "empagliflozin"
DATA_PATH_DAPAGLIFLOZIN = DATA_PATH_BASE / "dapagliflozin"
DATA_PATH_HCTZ = DATA_PATH_BASE / "hydrochlorothiazide"

DATA_PATHS = [
     DATA_PATH_EMPAGLIFLOZIN,
     DATA_PATH_DAPAGLIFLOZIN,
     DATA_PATH_HCTZ,
]
