"""Intestine model of empagliflozin: dissolution, absorption and fecal excretion."""

import numpy as np
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    Model,
    Parameter,
    RateRule,
    Reaction,
    Species,
    UnitDefinition,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.empagliflozin.models import annotations, templates
from pkdb_models.models.empagliflozin.models.fits import FITS


class U(templates.U):
    """Units of the intestine model."""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "empagliflozin_intestine",
    name="Model for empagliflozin absorption in the small intestine",
    notes="""
    # Model for empagliflozin absorption

    - Dissolution of the oral dose and absorption of empagliflozin; the
      fraction absorbed `F_emp_abs` of the empagliflozin of the lumen is
      absorbed, the rest is excreted in the feces after the transit through
      the intestine (`METEXC_k`). After an oral solution of 50 mg of
      14C-empagliflozin 41.2 % of the radioactivity is recovered in the feces,
      34.2 % of the dose as unchanged empagliflozin, no glucuronide; 54.4 % in
      the urine [Chen2015a]. The fecal excretion is slow: 10.9, 29.0, 37.9 and
      39.9 % of the dose at 48, 72, 96 and 120 hr.
    - Food (`fed`) slows the absorption (`FED_f_absorption`): a high fat
      breakfast lowers Cmax by 29 % and delays tmax from 1.0 to 2.5 hr, the AUC
      is lowered by 10 % [Seman2013]; a Chinese high fat meal lowers Cmax by
      20 % and the AUC by 5 % [Chen2020].
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "NCIT:C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine
        (BQB.HAS_PROPERTY, "NCIT:C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
    ],
)

_m.compartments = [
    Compartment(
        "Vext",
        1.0,
        name="plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["plasma"],
    ),
    Compartment(
        "Vgu",
        1.2825,  # 0.0171 [l/kg] * 75 kg
        name="intestine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["gu"],
    ),
    Compartment(
        "Vlumen",
        1.2825 * 0.9,  # 0.0171 [l/kg] * 75 kg * 0.9,
        name="intestinal lumen (inner part of intestine)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["apical"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["stomach"],
    ),
]


_m.species = [
    Species(
        "emp_stomach",
        metaId="meta_emp_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name="empagliflozin (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
        boundaryCondition=True,
    ),
    Species(
        "emp_lumen",
        initialConcentration=0.0,
        name="empagliflozin (intestinal lumen)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
        port=True,
    ),
    Species(
        "emp_ext",
        initialConcentration=0.0,
        name="empagliflozin (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
        port=True,
    ),
    Species(
        "emp_feces",
        initialConcentration=0.0,
        name="empagliflozin (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        "F_emp_abs",
        **FITS.parameter("GU__F_emp_abs", unit=U.dimensionless),
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name="fraction absorbed empagliflozin",
        notes="""
        Fraction absorbed: `F_emp_abs` of the empagliflozin of the intestinal
        lumen is absorbed, `(1 - F_emp_abs)` is excreted in the feces. The
        absolute bioavailability is not known; the radioactivity of the urine
        (54.4 %) and the oxidative metabolites of the feces (3.6 %) give a
        fraction absorbed of about 0.58 of the oral solution [Chen2015a].
        """,
    ),
    Parameter(
        "EMPABS_k",
        **FITS.parameter("GU__EMPABS_k", unit=U.per_min),
        name="rate of empagliflozin absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        First-order rate of the empagliflozin leaving the lumen, absorbed
        (`F_emp_abs`) or to the fecal transit: tmax 1.0-1.5 hr fasted
        [Seman2013, Heise2013].
        """,
    ),
    Parameter(
        "METEXC_k",
        **FITS.parameter("GU__METEXC_k", unit=U.per_min),
        name="rate of the fecal transit",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        Transit of the empagliflozin which is not absorbed through the
        intestine into the feces, `n_chain` steps with this rate (mean transit
        time of about 70 hr): 10.9, 29.0, 37.9 and 39.9 % of the dose in the
        feces at 48, 72, 96 and 120 hr [Chen2015a].
        """,
    ),
    Parameter(
        "fed",
        0.0,
        unit=U.dimensionless,
        name="fed state (0 fasted, 1 fed)",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        The dose is taken in the fasted state (0) or with a meal (1). The effect
        of food on the absorption rate is `FED_f_absorption`.
        """,
    ),
    Parameter(
        "FED_f_absorption",
        **FITS.parameter("GU__FED_f_absorption", unit=U.dimensionless),
        name="fed effect on the absorption rate",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        Factor of the absorption rate with food: Cmax 0.71 and the AUC 0.90 of
        the fasted state, tmax 2.5 instead of 1.0 hr [Seman2013]; Cmax 0.80,
        AUC 0.95 [Chen2020]. The fraction absorbed is not changed.
        """,
    ),
    Parameter(
        "f_absorption",
        np.nan,
        unit=U.dimensionless,
        constant=False,
        name="factor of the absorption rate",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
]

_m.rules.extend(
    [
        AssignmentRule(
            "f_absorption",
            "1 dimensionless + fed * (FED_f_absorption - 1 dimensionless)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "absorption",
            value="f_absorption * EMPABS_k * Vlumen * emp_lumen",
            unit=U.mmole_per_min,
            name="empagliflozin leaving the lumen",
        ),
    ]
)

_m.reactions = [
    Reaction(
        "EMPABS",
        name="absorption empagliflozin (EMPABS)",
        equation="emp_lumen -> emp_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("F_emp_abs * absorption", U.mmole_per_min),
    ),
    Reaction(
        sid="EMPEXC",
        name="empagliflozin to the fecal transit (EMPEXC)",
        compartment="Vlumen",
        equation="emp_lumen -> emp_intestine_0",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "(1 dimensionless - F_emp_abs) * absorption",
            U.mmole_per_min,
        ),
    ),
]

#: the steps of the fecal transit
n_chain = 5
for k in range(n_chain):
    _m.compartments.append(
        Compartment(
            f"Vintestine_{k}",
            1.2825 * 0.9 / n_chain,  # 0.0171 [l/kg] * 75 kg * 0.9,
            name=f"intestinal lumen (transit {k})",
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            unit=U.liter,
            constant=True,
            annotations=annotations.compartments["gu_lumen"],
        )
    )
    _m.species.append(
        Species(
            f"emp_intestine_{k}",
            initialConcentration=0.0,
            name=f"empagliflozin (intestinal transit {k})",
            compartment=f"Vintestine_{k}",
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["emp"],
        )
    )
    product = f"emp_intestine_{k + 1}" if k < n_chain - 1 else "emp_feces"
    _m.reactions.append(
        Reaction(
            sid=f"EMPEXC_{k}",
            name=f"fecal transit empagliflozin {k}",
            compartment=f"Vintestine_{k}",
            equation=f"emp_intestine_{k} -> {product}",
            sboTerm=SBO.TRANSPORT_REACTION,
            formula=(
                f"METEXC_k * Vintestine_{k} * emp_intestine_{k}",
                U.mmole_per_min,
            ),
        )
    )


_m.parameters.extend(
    [
        Parameter(
            "PODOSE_emp",
            0,
            U.mg,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="oral dose empagliflozin [mg]",
            port=True,
        ),
        Parameter(
            "Ka_dis_emp",
            2.0,
            U.per_hr,
            constant=True,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="dissolution rate empagliflozin",
            port=True,
        ),
        Parameter(
            "Mr_emp",
            450.909,
            U.g_per_mole,
            constant=True,
            name="Molecular weight empagliflozin [g/mole]",
            sboTerm=SBO.MOLECULAR_MASS,
            port=True,
        ),
    ]
)

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid="dissolution_emp",
            name="dissolution empagliflozin",
            formula=(
                "Ka_dis_emp/60 min_per_hr * PODOSE_emp/Mr_emp",
                U.mmole_per_min,
            ),
            equation="emp_stomach -> emp_lumen",
            compartment="Vgu",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """,
        ),
    ]
)
_m.rate_rules.append(
    RateRule("PODOSE_emp", "-dissolution_emp * Mr_emp", U.mg_per_min),
)

model_intestine = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.empagliflozin import MODEL_BASE_PATH

    results = create_model(
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        model=model_intestine,
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
