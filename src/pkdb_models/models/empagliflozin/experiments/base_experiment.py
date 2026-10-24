"""Base class of the empagliflozin simulation experiments."""

from collections.abc import Iterable, Mapping
from typing import ClassVar, NamedTuple

import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.impairment import RENAL_COLORS_END_STAGE, RENAL_FUNCTION_GFR
from sbmlsim import Q
from sbmlsim.units import Quantity

from pkdb_models.models.empagliflozin import MODEL
from pkdb_models.models.empagliflozin.empagliflozin_pk import (
    calculate_empagliflozin_pk,
)


class MolecularWeights(NamedTuple):
    """Molecular weights for the conversion of the data."""

    emp: Quantity
    eg: Quantity


class EmpagliflozinSimulationExperiment(PKPDSimulationExperiment):
    """Base class of the empagliflozin simulation experiments."""

    pkpd_model = MODEL
    selections = (
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
        "Afeces_emptot",
        # glucose
        "KI__glc_urine",
        "KI__GLCEX",
        "KI__UGE",
        "KI__RTG",
        "[KI__fpg]",
        "KI__pg",
        # renal excretion rate
        "KI__EMPEX",
        # cases
        "KI__f_renal_function",
        "KI__GFR",
        "f_cirrhosis",
        "GU__fed",
        "GU__FED_f_absorption",
    )

    font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 20}
    scan_font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 9
    suptitle_font_size = 25

    # healthy reference values
    #: fasting plasma glucose of healthy subjects [mM], `KI__fpg_healthy` of
    #: the model; the filtered glucose is the mean plasma glucose of the day,
    #: `KI__fpg` + `KI__PG_delta` (0.7 mM, 5.7 mM in healthy subjects)
    fpg_healthy = 5.0
    #: GFR of healthy kidneys at 75 kg [ml/min], `GFR_base` of the model; a
    #: measured creatinine clearance or GFR of X ml/min is the renal function
    #: X/101
    gfr_healthy = 101.0
    bodyweight_healthy = 75  # [kg] (healthy reference)

    # labels
    label_time = "Time"
    label_emp = "Empagliflozin"
    label_eg = "Empagliflozin-glucuronide\n"
    label_emptot = "Empagliflozin Total\n"

    label_glc = "Glucose"
    label_pg = "Mean Glucose"
    label_uge = "UGE"
    label_rtg = "RTG"

    label_emp_plasma = label_emp + " Plasma\n"
    label_emp_urine = label_emp + " Urine\n"
    label_eg_urine = label_eg + " Urine"
    label_emptot_urine = label_emptot + " Urine"

    label_glc_urine = label_glc + " Urine"

    label_emp_feces = label_emp + " Feces\n"
    label_emptot_feces = label_emptot + " Feces"

    labels: ClassVar[Mapping[str, str]] = {
        "time": "time",
        "[Cve_emp]": label_emp,
        "[Cve_eg]": label_eg,
        "[Cve_emptot]": label_emptot,
        "Aurine_emp": label_emp_urine,
        "Aurine_eg": label_eg_urine,
        "Aurine_emptot": label_emptot_urine,
        "Afeces_emp": label_emp_feces,
        "Afeces_emptot": label_emptot_feces,
        "[KI__fpg]": label_glc,
        "KI__pg": label_pg,
        "KI__EMPEX": f"{label_emp}\nExcretion Urine",
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
    unit_emptot_feces = unit_metabolite_feces

    unit_glc = "mM"
    unit_glc_urine = "mole"
    unit_uge = "g"
    unit_rtg = "mM"

    units: ClassVar[Mapping[str, str]] = {
        "time": unit_time,
        "[Cve_emp]": unit_emp,
        "[Cve_eg]": unit_eg,
        "[Cve_emptot]": unit_emptot,
        "Aurine_emp": unit_emp_urine,
        "Aurine_eg": unit_eg_urine,
        "Aurine_emptot": unit_emptot_urine,
        "Afeces_emp": unit_emp_feces,
        "Afeces_emptot": unit_emptot_feces,
        "[KI__fpg]": unit_glc,
        "KI__pg": unit_glc,
        "KI__EMPEX": "µmol/min",
        "KI__glc_urine": unit_glc_urine,
        "KI__UGE": unit_uge,
        "KI__RTG": unit_rtg,
        "KI__GLCEX": "mmol/min",
    }

    # ----------- Fasting/food -----
    #: the fed state of the intestine (`GU__fed`): food slows the absorption
    #: (`GU__FED_f_absorption`), Cmax 0.71, AUC 0.90 [Seman2013]
    fasting_map: ClassVar[dict[str, float]] = {
        "NR": 0.0,
        "fasted": 0.0,
        "fed": 1.0,
    }

    fasting_colors: ClassVar[dict[str, str]] = {
        "NR": "black",
        "fasted": "black",
        "fed": "tab:red",
    }

    # ----------- Renal function --------------
    #: the classes of the GFR with the end stage renal disease; a measured
    #: creatinine clearance or GFR of X ml/min is the renal function X/101
    renal_map = RENAL_FUNCTION_GFR
    renal_colors = RENAL_COLORS_END_STAGE

    # ----------- Interactions --------------
    #: the activities of the transporters and enzymes of a coadministered drug
    #: (`LI__f_oatp` hepatic uptake, `KI__f_oat3` renal uptake for the tubular
    #: secretion and the renal glucuronidation, `f_ugt` glucuronidation), from
    #: the mechanism and the observed change of the exposure, never fitted;
    #: the arms are validation data [Macha2014]:
    #: - probenecid 500 mg BID inhibits OAT1/OAT3: the renal clearance falls
    #:   from 32.7 to 15.2 ml/min (the filtration), the AUC rises to 1.53
    #:   (model with the adopted fit: CLR 31.6 to 15.8 ml/min, AUC 1.59)
    #: - rifampicin 600 mg single dose inhibits OATP1B1/1B3: the AUC rises to
    #:   1.35 (model 1.35; Cmax 1.75 and the shorter half-life, 11.8 to
    #:   7.5 hr, an effect on the absorption, are not described)
    #: - gemfibrozil 600 mg BID inhibits OATP1B1 and OAT3 (gemfibrozil and its
    #:   glucuronide): the AUC rises to 1.59, Cmax 1.15 (model 1.57, 1.11)
    #: The factors are set by hand against the observed ratios with the
    #: adopted fit.
    #: Sitagliptin, linagliptin, pioglitazone, lobeglitazone, evogliptin,
    #: verapamil, ramipril, warfarin, hydrochlorothiazide and torasemide do not
    #: change the exposure of empagliflozin (AUC ratios 0.97-1.10): no factors.
    interactions: ClassVar[dict[str, dict[str, float]]] = {
        "probenecid": {"KI__f_oat3": 0.1},
        "rifampicin": {"LI__f_oatp": 0.35},
        "gemfibrozil": {"LI__f_oatp": 0.6, "KI__f_oat3": 0.5},
    }

    def interaction_changes(self, drug: str) -> dict[str, Quantity]:
        """Get the changes of the model of a coadministered drug.

        Args:
            drug: the drug, a key of `interactions`.

        Returns:
            The activities of the transporters and enzymes by parameter id.
        """
        return {
            pid: Q(value, "dimensionless")
            for pid, value in self.interactions[drug].items()
        }

    # ----------- Dose map --------------
    # from light to dark with the hue turning from yellow over red to violet:
    # the lightness (OKLab) falls by 0.055 from dose to dose, adjacent doses
    # differ by about 7 (OKLab distance x 100)
    dose_colors: ClassVar[dict[float, str]] = {
        0: "black",
        0.5: "#F4D665",
        1: "#F9B946",
        2.5: "#FA9A3A",
        5: "#F47C44",
        10: "#E76157",
        25: "#D24C6A",
        50: "#B73C7A",
        100: "#973384",
        200: "#752D88",
        400: "#542783",
        800: "#322276",
    }

    # ----------- Glucose map --------------
    glucose_colors: ClassVar[dict[int, str]] = {
        5: "black",
        6: "#F5B8D8",
        7: "#EFA3CB",
        8: "#E88EBE",
        9: "#E279B1",
        10: "#DB64A4",
        11: "#D44F97",
    }

    # ----------- Cmaps --------------
    @property
    def renal_cmap(self) -> LinearSegmentedColormap:
        """Colormap of the renal function, from impaired (dark) to normal (light)."""
        return LinearSegmentedColormap.from_list(
            "renal_function",
            ["#013817", "#006d2c", "#2ca25f", "#66c2a4", "#8AEDCC", "#D5F0E7"],
            N=256,
        )

    @property
    def cirrhosis_cmap(self) -> LinearSegmentedColormap:
        """Colormap of the degree of cirrhosis, from none (light) to severe (dark)."""
        return LinearSegmentedColormap.from_list(
            "cirrhosis_severity",
            ["#CEE2F0", "#74a9cf", "#2b8cbe", "#045a8d", "#003352"],
            N=256,
        )

    @property
    def doses_cmap(self) -> LinearSegmentedColormap:
        """Colormap of the doses, from low (light) to high (dark)."""
        return LinearSegmentedColormap.from_list(
            "doses",
            [
                "#FFE5B4",
                "#FFD18A",
                "#FFBD66",
                "#FFA53A",
                "#FF931F",
                "#FF7A16",
                "#FF660D",
                "#F4530A",
                "#E04407",
                "#CC3705",
            ],
            N=256,
        )

    # keys are the parameters of the non-compartmental analysis of pkpdutils,
    # see `empagliflozin_pk`
    pk_labels: ClassVar[Mapping[str, str]] = {
        "auc_last": "AUCend",
        "auc_inf_obs": "AUCinf",
        "tmax": "Tmax",
        "cl_f": "Apparent Clearance",
        "cl_renal": "Renal Clearance",
        "cmax": "Cmax",
        "thalf": "Half-life",
        "lambda_z": "Elimination Rate",
        "vz_f": "Apparent Volume",
        "Aurine_emp": "Empagliflozin Urine",
    }

    #: the units of the studies, which report empagliflozin in nmol/l
    pk_units: ClassVar[Mapping[str, str]] = {
        "auc_last": "nmol/l*hr",
        "auc_inf_obs": "nmol/l*hr",
        "tmax": "hr",
        "cl_f": "ml/min",
        "cl_renal": "ml/min",
        "cmax": "nmol/l",
        "thalf": "hr",
        "lambda_z": "1/hr",
        "vz_f": "l",
        "Aurine_emp": "µmole",
    }

    # molecular weights of empagliflozin and its glucuronide
    Mr: ClassVar[MolecularWeights] = MolecularWeights(
        emp=Q(450.909, "g/mole"),
        eg=Q(627, "g/mole"),
    )

    def calculate_empagliflozin_pk(
        self, scans: Iterable[str] | None = None
    ) -> dict[str, pd.DataFrame]:
        """Calculate the pharmacokinetic parameters of the scans.

        Args:
            scans: keys of the simulations, all simulations by default.

        Returns:
            The parameters of every scan by simulation key, per draw of the
            fitted parameters in a run over them (`uncertainty`).
        """
        keys = self._simulations if scans is None else scans
        across = None if self.uncertainty is None else self.uncertainty.id
        return {
            key: calculate_empagliflozin_pk(
                experiment=self, res=self.results[f"task_{key}"], across=across
            )
            for key in keys
        }
