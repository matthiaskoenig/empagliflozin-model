"""Kidney model for the SGLT2 inhibitors."""
import libsbml
import numpy as np
from sbmlutils.converters import odefac
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.empagliflozin.models import annotations
from pkdb_models.models.empagliflozin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")



mid = "empagliflozin_kidney"
version = 3

_m = Model(
    sid=mid,
    name="Model for renal empagliflozin excretion.",
    notes=f"""
    Model for renal empagliflozin excretion.
    
    **version** {version}
    
    ## Changelog
    
    **version 1**
    
    - initial model
        
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=annotations.model + [
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7203"),  # kidney
        (BQB.OCCURS_IN, "bto/BTO:0000671"),  # kidney
        (BQB.OCCURS_IN, "NCIT:C12415"),  # kidney

        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
        (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
    ]
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True
    ),
    Compartment(
        "Vki",
        value=0.3,  # 0.4 % of bodyweight
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True
    ),
    Compartment(
        "Vmem",
        value=np.nan,
        unit=U.m2,
        name="plasma membrane",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma membrane"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vurine",
        1.0,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["urine"],
    ),

]

# ---------------------------------------------------------------------------------------------------------------------
# Pharmacokinetics
# ---------------------------------------------------------------------------------------------------------------------
_m.species = [
    Species(
        "emp_ext",
        name="empagliflozin (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
        port=True
    ),
    Species(
        "emp",
        name="empagliflozin (kidney)",
        initialConcentration=0.0,
        compartment="Vki",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
        port=True
    ),
    Species(
        "emp_urine",
        name="empagliflozin (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
        port=True
    ),
    Species(
        "eg_ext",
        name="empagliflozin-glucuronide (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
        port=True
    ),
    Species(
        "eg",
        name="empagliflozin-glucuronide (kidney)",
        initialConcentration=0.0,
        compartment="Vki",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
        port=True
    ),
    Species(
        "eg_urine",
        name="empagliflozin-glucuronide (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
        port=True
    ),
]

_m.parameters.extend([
    Parameter(
        "f_renal_function",
        name="parameter for renal function",
        value=1.0,
        unit=U.dimensionless,
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""scaling factor for renal function. 1.0: normal renal function; 
        <1.0: reduced renal function
        """
    ),
])

_m.reactions = [
    Reaction(
        sid="EMPIM",
        name="empagliflozin import (EMPIM)",
        equation="emp_ext <-> emp",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "EMPIM_k",
                100.0,  # fast
                U.per_min,
                name="rate empagliflozin import",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            ),
        ],
        formula=(
            "EMPIM_k * Vki * (emp_ext - emp)",
            U.mmole_per_min
        ),
    ),
    Reaction(
        sid="EGIM",
        name="empagliflozin-glucuronide import (EGIM)",
        equation="eg_ext <-> eg",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "EGIM_k",
                100.0,  # fast
                U.per_min,
                name="rate empagliflozin-glucuronide import",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            ),
        ],
        formula=(
            "EGIM_k * Vki * (eg_ext - eg)",
            U.mmole_per_min
        ),
    ),
    Reaction(
        sid="EMP2EG",
        name="empagliflozin glucuronidation (EMP2EG) UGT",
        equation="emp -> eg",
        compartment="Vki",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "EMP2EG_Vmax",
                0.04,
                U.mmole_per_min_l,
                name="Vmax empagliflozin conversion",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "EMP2EG_Km_emp",
                0.02,
                U.mM,
                name="Km empagliflozin UGT",
                sboTerm=SBO.MICHAELIS_CONSTANT,
            ),
            Parameter(
                "f_ugt",
                1,
                U.dimensionless,
                name="scaling factor UGT activity",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""Scaling factor to vary UGT activity.
                1.0: unchanged activity; < 1.0 decreased activity; >1.0 increased activity.
                """
            )
        ],
        formula=(
            "f_ugt * EMP2EG_Vmax * Vki * emp/(emp + EMP2EG_Km_emp)"
        ),
    ),
    Reaction(
        sid="EMPEX",
        name="empagliflozin excretion (EMPEX)",
        equation="emp_ext -> emp_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "EMPEX_k",
                0.003,
                U.per_min,
                name="rate urinary excretion of empagliflozin",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * EMPEX_k * Vki * emp_ext"
        ),
    ),
    Reaction(
         sid="EGEX",
         name="empagliflozin-glucuronide excretion (EGEX)",
         equation="eg_ext -> eg_urine",
         compartment="Vki",
         sboTerm=SBO.TRANSPORT_REACTION,
         pars=[
                Parameter(
                    "EGEX_k",
                    0.010,
                    U.per_min,
                    name="rate urinary excretion of empagliflozin-glucuronide",
                    sboTerm=SBO.KINETIC_CONSTANT,
                ),
            ],
        formula=(
            "f_renal_function * EGEX_k * Vki * eg_ext"
        ),
    )
]

# ---------------------------------------------------------------------------------------------------------------------
# Pharmacodynamics
# ---------------------------------------------------------------------------------------------------------------------
_m.species.extend([
    Species(
        "fpg",
        name="fasting plasma glucose (FPG)",
        initialConcentration=5.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        boundaryCondition=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["glc"],
        notes="""
        Fasting plasma glucose (FPG).
        Model depends on the fasting plasma glucose:
        - as a proxy of the daily glucose variation in the glucose excretion
        - to set the glucose dependent RTG values (i.e. T2DM have higher SGLT2 amounts, and consequently higher RTG values)
        """
    ),
    Species(
        "glc_urine",
        name="glucose (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["glc"],
        port=True
    ),
])

_m.parameters.extend([
    Parameter(
        "Mr_glc",
        180,
        U.g_per_mole,
        name=f"Molecular weight glc [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
    ),
    Parameter(
        "cf_mg_per_g",
        1000,
        U.mg_per_g,
        name=f"Conversion factor mg per g",
    ),
    Parameter(
        "cf_ml_per_l",
        1000,
        U.ml_per_l,
        name=f"Conversion factor ml per l",
    ),
    Parameter(
        "GFR_healthy",
        100,
        U.ml_per_min,
        name=f"Glomerular filtration rate (healthy)",
    ),
    Parameter(
        "RTG_E50",
        71.9E-7,
        # [10 - 200E-6] for optimization; 32E-3/444.518 = 71.9E-6,  # [mM]   [ng/ml]/[g/mole] = [nmole/ml] = µmole/l
        U.mM,
        name="EC50 reduction in RTG",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
    Parameter(
        "RTG_gamma",
        1,  # [1 - 4] for optimization
        U.dimensionless,
        name="hill coefficient reduction in RTG",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
    Parameter(
        "RTG_base",
        12.5,  # [10.5 - 14] for optimization
        U.mM,
        name=f"Baseline RTG value",
        notes="""Typical RTG value without SGLT2 inhibitors in healthy subjects.

        This corresponds to the SGLT2 concentrations.
        """
    ),
    Parameter(
        "RTG_m_fpg",
        0.5,  # [0.2 - 1] for optimization
        U.dimensionless,
        name=f"FPG effect on RTG",
        notes="""Effect of the FPG on the change in RTG_base."""
    ),
    Parameter(
        "RTG_max_inhibition",
        0.75,  # [0 - 1] for optimization
        U.dimensionless,
        name=f"RTG maximum inhibition",
        notes="maximum inhibition of RTG via SGLT2"
    ),
    Parameter(
        "fpg_healthy",
        5,
        U.mM,
        name=f"fasting plasma glucose (healthy)",
    ),
])

_m.rules.extend([
    AssignmentRule(
        "RTG_fpg",
        "RTG_base + RTG_m_fpg * (fpg - fpg_healthy)",
        U.mM,
        name=f"RTG value (FPG)",
        notes="""
        FPG dependent base RTG value. SGLT2 is induced in T2DM [Rahmoune2005] 
        """
    ),
    AssignmentRule(
        "RTG_delta",
        "RTG_fpg * RTG_max_inhibition",  # [7 - 10]
        U.mM,
        name=f"RTG value",
        notes="""
        ΔRTG, maximum reduction in RTG due to SGLT2 inhibition.
        Normally around 7 - 10 mM. 
        """
    ),
    AssignmentRule(
        variable="RTG",
        value="RTG_fpg - RTG_delta * power(emp_ext, RTG_gamma)/ (power(RTG_E50, RTG_gamma) + power(emp_ext, RTG_gamma))",
        unit=U.mM,
        name="renal threshold glucose (RTG)",
        notes="""
        Renal threshold glucose.

        The renal threshold for glucose (RTG) is the plasma glucose concentration at which tubular reabsorption of 
        glucose begins to saturate; glucose is excreted into the urine in direct proportion 
        to the glucose concentration above this threshold.
        12.3 - 12.7 mM RTG (placebo) [Devineni2012]
        """
    ),
    AssignmentRule(
        variable="GFR",
        value="f_renal_function * GFR_healthy",
        unit=U.ml_per_min,
        name="glomerular filtration rate",
    ),
])

# Glucose excretion (UGE)

_m.reactions.extend([
    Reaction(
        sid="GLCEX",
        name="glucose excretion (GLCEX)",
        equation="fpg -> glc_urine [emp_ext]",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        #  piecewise      | x1, y1, [x2, y2,] [...] [z] | A piecewise function: if (y1), x1.  Otherwise, if (y2), x2, etc.  Otherwise, z.
        formula=(
            # [ml/min]/[ml/l] *[mmole/l] = [mmole/min]
            "piecewise(GFR/cf_ml_per_l * (fpg - RTG), fpg > RTG, 0 mmole_per_min)"
        # FIXME: no scaling with liver volume Vki? (GFR should scale with kidney volume)
        )
    ),
])

_m.rules.extend([
    AssignmentRule(
        "UGE", "glc_urine * Mr_glc/cf_mg_per_g", unit=U.gram,
        name="urinary glucose excretion (UGE)",
        notes="""
        Urinary glucose excretion is calculated from cumulative amount of glucose in urine.
        """
    )
])

model_kidney = _m

import pandas as pd


def empagliflozin_layout(dx=200, dy=200) -> pd.DataFrame:
    """Layout definition for empagliflozin kidney model."""

    delta_y = 0.5 * dy
    delta_x = 0.7 * dx

    positions = [
        # sid, x, y
        ["fpg",       -0.3 * delta_x, 0],
        ["emp_ext",    1 * delta_x,   0],
        ["eg_ext",     3 * delta_x,   0],

        ["EMPIM",      1 * delta_x,   1 * delta_y],
        ["EGIM",       3 * delta_x,   1 * delta_y],

        ["emp",        1 * delta_x,   2 * delta_y],
        ["eg",         3 * delta_x,   2 * delta_y],

        ["EMP2EG",     2 * delta_x,   1.7 * delta_y],

        ["GLCEX",     -0.3 * delta_x, 3 * delta_y],
        ["EMPEX",      1 * delta_x,   3 * delta_y],
        ["EGEX",       3 * delta_x,   3 * delta_y],

        ["glc_urine", -0.3 * delta_x, 4 * delta_y],
        ["emp_urine",  1 * delta_x,   4 * delta_y],
        ["eg_urine",   3 * delta_x,   4 * delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df


def empagliflozin_annotations(dx=200, dy=200) -> list:
    COLOR_BLOOD = "#FF796C"
    COLOR_CELL = "#FFFFFF"
    COLOR_URINE = "#FF7F0E"
    delta_y = 0.5 * dy

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }
    xpos = -0.8 * dx
    width = 3.7 * dx

    annotations = [
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=-0.5 * delta_y, width=width, height=1.5 * delta_y,
            fill_color=COLOR_BLOOD, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=delta_y, width=width, height=2 * delta_y,
            fill_color=COLOR_CELL, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=3 * delta_y, width=width, height=1.5 * delta_y,
            fill_color=COLOR_URINE, **kwargs
        ),
    ]
    return annotations


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import MODEL_BASE_PATH
    from sbmlutils import cytoscape as cyviz
    from sbmlutils.converters import odefac
    from sbmlutils.factory import create_model

    # create SBML
    results = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(sbml_path=results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=empagliflozin_layout())
    cyviz.add_annotations(annotations=empagliflozin_annotations())
