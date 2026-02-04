"""
Reusable functionality for multiple simulation experiments.
"""
from collections import namedtuple
from typing import Dict
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from pkdb_models.models.empagliflozin import MODEL_PATH
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task
from pkdb_models.models.empagliflozin.empagliflozin_pk import calculate_empagliflozin_pk

# Constants for conversion
MolecularWeights = namedtuple("MolecularWeights", "emp eg")


class EmpagliflozinSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 20}
    scan_font = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 9
    suptitle_font_size = 25

    # healthy reference values
    fpg_healthy = 6.5  # [mM] (healthy, mean daily glucose)
    gfr_healthy = 100.0  # [ml/min] (healthy reference)
    bodyweight_healthy = 75  # [kg] (healthy reference)


    # labels
    label_time = "Time"
    label_emp = "Empagliflozin"
    label_eg = "Empagliflozin-glucuronide\n"
    label_emptot = "Empagliflozin Total\n"

    label_glc = "Glucose"
    label_uge = "UGE"
    label_rtg = "RTG"

    label_emp_plasma = label_emp + " Plasma\n"
    label_emp_urine = label_emp + " Urine\n"
    label_eg_urine = label_eg + " Urine"
    label_emptot_urine = label_emptot + " Urine"

    label_glc_urine = label_glc + " Urine"

    label_emp_feces = label_emp + " Feces\n"
    label_eg_feces = label_eg + " Feces"
    label_emptot_feces = label_emptot + " Feces"

    labels: Dict[str, str] = {
        "time": "time",

        "[Cve_emp]": label_emp,
        "[Cve_eg]": label_eg,
        "[Cve_emptot]": label_emptot,

        "Aurine_emp": label_emp_urine,
        "Aurine_eg": label_eg_urine,
        "Aurine_emptot": label_emptot_urine,

        "Afeces_emp": label_emp_feces,
        "Afeces_eg": label_eg_feces,
        "Afeces_emptot": label_emptot_feces,

        "[KI__fpg]": label_glc,
        "KI__glc_urine": label_glc_urine,
        "KI__UGE": label_uge,
        "KI__RTG": label_rtg,
    }

    # units
    unit_time = "hr"
    unit_metabolite = "µM"
    unit_metabolite_urine = "µmole"
    unit_metabolite_feces = "µmole"
    unit_emp = unit_metabolite
    unit_eg = unit_metabolite
    unit_emptot = unit_metabolite

    unit_emp_urine = unit_metabolite_urine
    unit_eg_urine = unit_metabolite_urine
    unit_emptot_urine = unit_metabolite_urine

    unit_emp_feces = unit_metabolite_feces
    unit_eg_feces = unit_metabolite_feces
    unit_emptot_feces = unit_metabolite_feces

    unit_glc = "mM"
    unit_glc_urine = "mole"
    unit_uge = "g"
    unit_rtg = "mM"

    units: Dict[str, str] = {
        "time": unit_time,
        "[Cve_emp]": unit_emp,
        "[Cve_eg]": unit_eg,
        "[Cve_emptot]": unit_emptot,

        "Aurine_emp": unit_emp_urine,
        "Aurine_eg": unit_eg_urine,
        "Aurine_emptot": unit_emptot_urine,

        "Afeces_emp": unit_emp_feces,
        "Afeces_eg": unit_eg_feces,
        "Afeces_emptot": unit_emptot_feces,

        "[KI__fpg]": unit_glc,
        "KI__glc_urine": unit_glc_urine,
        "KI__UGE": unit_uge,
        "KI__RTG": unit_rtg,
        "KI__GLCEX": "mmol/min",
    }

    # ----------- Fasting/food -----
    fasting_map = {
        "fasted": 1.0,
        "fed": 0.9,
    }
    fasting_colors = {
        "fasted": "black",
        "fed": "tab:blue",
    }

    # ----------- Renal map --------------
    renal_map = {
        "Normal renal function": 101.0 / 101.0,     # 1.0,
        "Mild renal impairment": 69.5 / 101.0,      # 0.69
        "Moderate renal impairment": 32.5 / 101.0,  # 0.32
        "Severe renal impairment": 19.5 / 101.0,    # 0.19
        # "End stage renal disease": 10.5 / 101.0,  # 0.1
    }
    renal_colors = {
        "Normal renal function": "black",
        "Mild renal impairment": "#66c2a4",
        "Moderate renal impairment": "#2ca25f",
        "Severe renal impairment": "#006d2c",
        "End stage renal disease": "#006d5e"
    }

    # ----------- Cirrhosis map --------------
    cirrhosis_map = {
        "Control": 0,
        "Mild cirrhosis": 0.3994897959183674,      # CPT A
        "Moderate cirrhosis": 0.6979591836734694,  # CPT B
        "Severe cirrhosis": 0.8127551020408164,    # CPT C
    }
    cirrhosis_colors = {
        "Control": "black",
        "Mild cirrhosis": "#74a9cf",      # CPT A
        "Moderate cirrhosis": "#2b8cbe",  # CPT B
        "Severe cirrhosis": "#045a8d",    # CPT C
    }

    # ----------- Dose map --------------
    dose_colors = {
        0: "black",
        0.5: "#FFF44F",
        1: "#FFD700",
        2.5: "#FFA500",
        5: "#FF8C00",
        10: "#FF6B35",
        25: "#FF5722",
        50: "#F4511E",
        100: "#E53935",
        200: "#C62828",
        400: "#B71C1C",
        800: "#8B0000",
    }

    # ----------- Glucose map --------------
    glucose_colors = {
        5:  "black",
        6: "#F5B8D8",
        7: "#EFA3CB",
        8: "#E88EBE",
        9: "#E279B1",
        10: "#DB64A4",
        11: "#D44F97",
    }

    # ----------- Cmaps --------------
    @property
    def renal_cmap(self):
        return LinearSegmentedColormap.from_list(
            'renal_function',
            ["#013817", "#006d2c", "#2ca25f", "#66c2a4", "#8AEDCC", "#D5F0E7"],
            N=256
        )

    @property
    def cirrhosis_cmap(self):
        return LinearSegmentedColormap.from_list(
            'cirrhosis_severity',
            ["#CEE2F0", "#74a9cf", "#2b8cbe", "#045a8d", "#003352"],
            N=256
        )

    @property
    def doses_cmap(self):
        return LinearSegmentedColormap.from_list(
            'doses',
            ["#FFE5B4", "#FFD18A", "#FFBD66", "#FFA53A", "#FF931F",
             "#FF7A16","#FF660D", "#F4530A", "#E04407", "#CC3705"],
            N=256
        )

    pk_labels = {
        "auc": "AUCend",
        "aucinf": "AUCinf",
        "tmax": "Tmax",
        "cl": "Total Clearance",
        "cl_hepatic": "Hepatic Clearance",
        "cmax": "Cmax",
        "thalf": "Half-life",
        "kel": "kel",
        "vd": "vd",
        "Aurine_can": "Canagliflozin Urine",
    }

    pk_units = {
        "auc": "µmole/l*hr",
        "aucinf": "ng/ml*hr",
        "tmax": "hr",
        "cl": "ml/min",
        "cl_hepatic": "ml/min",
        "cmax": "ng/ml",
        "thalf": "hr",
        "kel": "1/hr",
        "vd": "l",
        "Aurine_can": "µmole",
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_PATH,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {
            # PK
            'ftissue_emp': Q_(0.3832934310804301, 'l/min'),  # [0.01 - 10]
            'Kp_emp': Q_(0.573747586082099, 'dimensionless'),  # [0.1 - 10]
            'GU__EMPABS_k': Q_(0.007439770316739214, '1/min'),  # [0.0001 - 1]
            'GU__METEXC_k': Q_(0.0008998014627372138, '1/min'),  # [1e-06 - 0.1]
            'LI__EMP2EG_Vmax': Q_(0.012490394598835162, 'mmol/min/l'),  # [0.001 - 100]
            'LI__EMP2EG_Km_emp': Q_(0.24667266695912698, 'mM'),  # [0.001 - 1]
            'LI__EGEX_k': Q_(0.014033858848302178, '1/min'),  # [0.001 - 100]
            'LI__EGBIEX_k': Q_(0.00735528424770721, '1/min'),  # [1e-05 - 1]
            'KI__EMPEX_k': Q_(0.11795401426241267, '1/min'),  # [0.0001 - 10]
            'KI__EGEX_k': Q_(1.4994487519741735, '1/min'),  # [0.0001 - 10]

            # PD
            'KI__RTG_E50': Q_(1.7063026384665662e-06, 'mM'),  # [1e-08 - 44]
            'KI__RTG_base': Q_(12.397702595815295, 'mM'),  # [9 - 14]
            'KI__RTG_max_inhibition': Q_(0.6842787312726607, 'dimensionless'),  # [0.2 - 1.0]
            'KI__RTG_m_fpg': Q_(0.8981650940449492, 'dimensionless'),  # [0.2 - 3]
        }

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return EmpagliflozinSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
                "time",

                # dosing
                "IVDOSE_emp",
                "PODOSE_emp",

                # venous
                "[Cve_emp]",
                "[Cve_eg]",
                "[Cve_emptot]",

                # urine
                "Aurine_emp",
                "Aurine_eg",
                "Aurine_emptot",

                # feces
                "Afeces_emp",
                "Afeces_eg",
                "Afeces_emptot",

                # glucose
                "KI__glc_urine",
                "KI__GLCEX",
                "KI__UGE",
                "KI__RTG",
                "[KI__fpg]",

                # cases
                'KI__f_renal_function',
                'f_cirrhosis',
            ]
        )
        return {}

    @property
    def Mr(self):
        return MolecularWeights(
            emp=self.Q_(450.909, "g/mole"),
            eg=self.Q_(627, "g/mole"),
        )

    def calculate_empagliflozin_pk(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate empagliflozin parameters for simulations (scans)"""
       pk_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_empagliflozin_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_empagliflozin_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       return pk_dfs
