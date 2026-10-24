"""Liver model of empagliflozin: hepatic uptake and glucuronidation."""

from sbmlutils.factory import (
    Compartment,
    Model,
    Parameter,
    Reaction,
    Species,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.empagliflozin.models import annotations, templates
from pkdb_models.models.empagliflozin.models.fits import FITS


class U(templates.U):
    """Units of the liver model."""


mid = "empagliflozin_liver"

_m = Model(
    sid=mid,
    name="Model for the hepatic metabolism of empagliflozin.",
    notes="""
    Model for the hepatic metabolism of empagliflozin.

    - Empagliflozin is taken up into the hepatocytes by OATP1B1 and OATP1B3
      and glucuronidated by UGT2B7, UGT1A3, UGT1A8 and UGT1A9 to the 2-O-, 3-O-
      and 6-O-glucuronides (EG), which are released into the plasma and
      excreted in the urine; no glucuronide is found in the feces [Chen2015a].
    - The hepatic clearance is a first-order clearance of the plasma of the
      liver: the uptake is the rate determining step of the low permeable
      empagliflozin (`f_oatp`), the glucuronidation (`f_ugt`) is far below its
      Km. It scales with the volume of the liver tissue, so with the tissue
      loss in cirrhosis.
    - The exposure rises with the Child-Pugh class (AUC ratio 1.23, 1.47, 1.75
      [Macha2014b]); with rifampicin (single dose, inhibition of OATP1B1/1B3)
      the AUC rises by 35 % [Macha2014].
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "NCIT:C12392"),
        (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
    ],
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
        port=True,
    ),
    Compartment(
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
        port=True,
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
        port=True,
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
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        "EMP2EG_k",
        **FITS.parameter("EMP2EG_k", unit=U.per_min),
        name="rate of the hepatic glucuronidation of empagliflozin",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""First-order rate of the hepatic clearance relative to the volume
        of the liver tissue; the whole-body model sets it.""",
        port=True,
    ),
    Parameter(
        "f_ugt",
        1,
        U.dimensionless,
        name="scaling factor UGT activity",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""Scaling factor of the UGT activity of the liver (and the kidney,
        the parameter of the whole-body model replaces both).
        1.0: normal activity; <1.0 decreased activity; >1.0 increased activity.
        """,
        port=True,
    ),
    Parameter(
        "f_oatp",
        1,
        U.dimensionless,
        name="scaling factor hepatic uptake (OATP1B1/1B3)",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""Scaling factor of the hepatic uptake of empagliflozin by
        OATP1B1 and OATP1B3, the rate determining step of the hepatic
        clearance. 1.0: normal activity; <1.0 inhibited (rifampicin,
        gemfibrozil).
        """,
    ),
]

_m.reactions = [
    Reaction(
        sid="EMP2EG",
        name="empagliflozin glucuronidation (EMP2EG)",
        equation="emp_ext -> eg_ext",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""Uptake of the empagliflozin of the plasma of the liver and its
        glucuronidation; the export of EG into the plasma is fast.""",
        formula=(
            "f_oatp * f_ugt * EMP2EG_k * Vli * emp_ext",
            U.mmole_per_min,
        ),
    ),
]

model_liver = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.empagliflozin import MODEL_BASE_PATH

    results = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
