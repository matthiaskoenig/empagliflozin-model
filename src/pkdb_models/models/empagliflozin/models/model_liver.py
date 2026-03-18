"""Liver model for SGLT2 inhibitors."""
import pandas as pd
import numpy as np

from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.empagliflozin.models import annotations
from pkdb_models.models.empagliflozin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    pass


mid = "empagliflozin_liver"
version = 1

_m = Model(
    sid=mid,
    name="Model for hepatic empagliflozin metabolism.",
    notes=f"""
    Model for empagliflozin metabolism.
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=annotations.model + [
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "NCIT:C12392"),

        (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
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
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
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
        "Vapical",
        np.nan,
        name="apical membrane",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["apical"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vbi",
        1.0,
        name="bile",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["bi"],
        port=True,
    ),
    Compartment(
        "Vlumen",
        1.2825 * 0.9,  # 0.0171 [l/kg] * 75 kg * 0.9, # FIXME: calculate from whole-body
        name="intestinal lumen (inner part of intestine)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),
]

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
        "emp",
        name="empagliflozin (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
    ),
    Species(
        "eg",
        name="empagliflozin-glucuronide (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
    ),

    Species(
        "eg_bi",
        initialConcentration=0.0,
        name="empagliflozin-glucuronide (bile)",
        compartment="Vbi",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
        notes="""
        Bile EG in amount.
        """,
    ),
    Species(
        "eg_lumen",
        initialConcentration=0.0,
        name="EG (lumen)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
        port=True,
    ),
]

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
            "EMPIM_k * Vli * (emp_ext - emp)",
            U.mmole_per_min
        ),
    ),
    Reaction(
        sid="EMP2EG",
        name="empagliflozin glucuronidation (EMP2EG) UGT",
        equation="emp -> eg",
        compartment="Vli",
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
            "f_ugt * EMP2EG_Vmax * Vli * emp/(emp + EMP2EG_Km_emp)"
        ),
    ),

    Reaction(
        sid="EGEX",
        name="EG export (EGEX)",
        equation="eg <-> eg_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "EGEX_k",
                10,
                U.per_min,
                name="rate EG export",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            )
        ],
        formula=(
            "EGEX_k * Vli * (eg - eg_ext)"
        )
    ),

    Reaction(
        sid="EGBIEX",
        name="EG bile export",
        equation="eg -> eg_bi",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        pars=[
            Parameter(
                "EGBIEX_k",
                0.0001,
                U.per_min,
                name="rate for EG export in bile",
                sboTerm=SBO.KINETIC_CONSTANT,
            )
        ],
        formula=(
            "EGBIEX_k * Vli * eg",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        "EGEHC",
        name="EG enterohepatic circulation",
        equation="eg_bi -> eg_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vlumen",
        formula=("EGBIEX", U.mmole_per_min),
    ),
]

model_liver = _m

def empagliflozin_layout(dx=200, dy=200) -> pd.DataFrame:
    """Layout definition."""

    delta_y = 0.5 * dy
    delta_x = 0.7 * dx

    positions = [
        # sid, x, y

        ["eg_ext",  0 * delta_x, 0 * delta_y],
        ["emp_ext", 1.5 * delta_x, 0 * delta_y],

        ["EGEX",   0 * delta_x, 1 * delta_y],
        ["EMPIM",  1.5 * delta_x, 1 * delta_y],

        ["eg",  0 * delta_x, 2 * delta_y],
        ["emp", 1.5 * delta_x, 2 * delta_y],

        ["EMP2EG", 0.9 * delta_x, 2.8 * delta_y],

        ["EGBIEX",  0 * delta_x, 3.0 * delta_y],

        ["eg_bi", 0.7 * delta_x, 3.5 * delta_y],

        ["EGEHC", 0.7 * delta_x, 4.5 * delta_y],

        ["eg_lumen", 0.7 * delta_x, 5.2 * delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)

    return df


def empagliflozin_annotations(dx=200, dy=200) -> list:
    COLOR_PLASMA = "#FF796C"
    COLOR_LIVER = "#FFFFFF"
    COLOR_BILE = "#F5F5C6"
    COLOR_LUMEN = "#CFEFFF"

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }

    xpos = -0.5 * dx
    width = 2.0 * dx
    delta_y = 0.5 * dy

    annotations = [
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=-0.5 * delta_y, width=width, height=1.5 * delta_y,
            fill_color=COLOR_PLASMA, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=1.0 * delta_y, width=width, height=2.0 * delta_y,
            fill_color=COLOR_LIVER, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=3.0 * delta_y, width=width, height=1.5 * delta_y,
            fill_color=COLOR_BILE, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=4.5 * delta_y, width=width, height=1.3 * delta_y,
            fill_color=COLOR_LUMEN, **kwargs
        ),
    ]

    return annotations


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import MODEL_BASE_PATH
    from sbmlutils import cytoscape as cyviz

    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    cyviz.visualize_sbml(sbml_path=results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=empagliflozin_layout())
    cyviz.add_annotations(annotations=empagliflozin_annotations())
