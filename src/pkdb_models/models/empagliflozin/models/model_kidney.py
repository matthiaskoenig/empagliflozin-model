"""Kidney model of empagliflozin: renal excretion, renal glucuronidation and glucose excretion."""

from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    Model,
    Parameter,
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
    """Units of the kidney model."""

    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")
    l_per_ml = UnitDefinition("l_per_ml", "l/ml")


mid = "empagliflozin_kidney"

_m = Model(
    sid=mid,
    name="Model for the renal excretion of empagliflozin and EG and the urinary glucose excretion.",
    notes="""
    Model for the renal excretion of empagliflozin and its glucuronide (EG), the
    renal glucuronidation of empagliflozin and the urinary glucose excretion.

    - 16-23 % of an oral dose of empagliflozin is excreted unchanged in the
      urine with a renal clearance of 28-39 ml/min [Seman2013, Sarashina2013,
      Macha2014, Macha2014f]; the glucuronides in the urine are 23 % of the
      dose [Chen2015a].
    - The renal excretion is relative to the glomerular filtration rate (GFR):
      the filtration of the unbound fraction (fu 0.14) and a tubular secretion
      by OAT3 (`f_oat3`); probenecid lowers the renal clearance from 32.7 to
      15.2 ml/min, the filtration [Macha2014].
    - The kidney glucuronidates empagliflozin (UGT1A9, UGT2B7) after its uptake
      by OAT3; the renal glucuronidation scales with the renal function (the
      functional renal mass), which makes the exposure of empagliflozin rise in
      renal impairment (AUC ratio 1.20, 1.23, 1.67 [Macha2014f]).
    - Empagliflozin lowers the renal threshold of glucose (RTG) by the
      inhibition of SGLT2; glucose above the RTG is excreted with the GFR.
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7203"),  # kidney
        (BQB.OCCURS_IN, "bto/BTO:0000671"),  # kidney
        (BQB.OCCURS_IN, "NCIT:C12415"),  # kidney
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
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
        "Vki",
        value=0.3,  # 0.4 % of bodyweight
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True,
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
        port=True,
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
    Species(
        "eg_urine",
        name="empagliflozin-glucuronide (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
        port=True,
    ),
    Species(
        "emptot_urine",
        name="total empagliflozin (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emptot"],
        notes="""Empagliflozin and EG in the urine, compared with the total
        radioactivity in the urine.""",
        port=True,
    ),
]

_m.parameters.extend(
    [
        Parameter(
            "f_renal_function",
            name="renal function",
            value=1.0,
            unit=U.dimensionless,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            >1.0: increased kidney function
            1.0: normal kidney function (healthy control)
            <1.0: decreased kidney function

            Renal function relative to the healthy glomerular filtration rate
            GFR_base, which is 101 ml/min at 75 kg; a measured creatinine
            clearance or GFR of X ml/min is the renal function X/101. It scales
            the filtration (renal excretion, glucose excretion) and the renal
            glucuronidation (the functional renal mass).
            """,
        ),
        Parameter(
            "GFR_base",
            101,
            U.ml_per_min,
            name="glomerular filtration rate (healthy)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
        ),
        Parameter(
            "f_ugt",
            1,
            U.dimensionless,
            name="scaling factor UGT activity",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""Scaling factor of the UGT activity of the kidney (and the
            liver, the parameter of the whole-body model replaces both).
            1.0: normal activity; <1.0 decreased activity; >1.0 increased activity.
            """,
            port=True,
        ),
        Parameter(
            "f_oat3",
            1,
            U.dimensionless,
            name="scaling factor renal uptake (OAT3)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""Scaling factor of the basolateral uptake of empagliflozin
            into the proximal tubule by OAT3: the tubular secretion and the renal
            glucuronidation. Probenecid, an inhibitor of OAT1/OAT3, lowers the
            renal clearance of empagliflozin from 32.7 to 15.2 ml/min, the
            filtration of the unbound fraction [Macha2014].
            1.0: normal activity; <1.0 inhibited.
            """,
        ),
        Parameter(
            "fu_emp",
            0.14,
            U.dimensionless,
            name="fraction unbound in plasma empagliflozin",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""Empagliflozin is 80-86 % bound to plasma proteins
            [Chen2015a, Macha2014f], independent of the renal function.""",
        ),
        Parameter(
            "EMPSEC_f",
            **FITS.parameter("KI__EMPSEC_f", unit=U.dimensionless),
            name="tubular secretion of empagliflozin relative to its filtration",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The renal clearance of empagliflozin is the filtration of its unbound
            fraction, fu_emp * GFR (14 ml/min), and the tubular secretion by
            OAT3, `EMPSEC_f` times the filtration. The renal clearance of
            28-39 ml/min [Seman2013, Sarashina2013, Macha2014] is about twice the
            filtration.
            """,
        ),
        Parameter(
            "EGEX_f",
            **FITS.parameter("KI__EGEX_f", unit=U.dimensionless),
            name="renal clearance of EG relative to the GFR",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The renal clearance of the glucuronides is EGEX_f times the GFR; the
            glucuronides are about 0.2 of empagliflozin in the plasma (17.7 %
            vs 79.6 % of the radioactivity, AUC 0-12 hr) [Chen2015a], their
            kinetics are limited by their formation.
            """,
        ),
        Parameter(
            "EMP2EG_k",
            0.0566,
            U.per_min,
            name="rate of the renal glucuronidation of empagliflozin",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            First-order rate of the renal glucuronidation relative to the volume
            of the kidney; the whole-body model sets it from the hepatic rate
            and the renal fraction of the glucuronidation.
            """,
            port=True,
        ),
    ]
)

_m.rules.extend(
    [
        AssignmentRule(
            "GFR",
            "f_renal_function * GFR_base",
            unit=U.ml_per_min,
            name="glomerular filtration rate",
        ),
        AssignmentRule(
            "emptot_urine",
            "emp_urine + eg_urine",
            U.mmole,
            name="total empagliflozin (urine)",
        ),
    ]
)

_m.reactions = [
    Reaction(
        sid="EMP2EG",
        name="empagliflozin glucuronidation (EMP2EG)",
        equation="emp_ext -> eg_ext",
        compartment="Vki",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""The glucuronidation of the empagliflozin of the renal plasma
        after the uptake by OAT3; the export of EG is fast.""",
        formula=(
            "f_renal_function * f_oat3 * f_ugt * EMP2EG_k * Vki * emp_ext",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="EMPEX",
        name="empagliflozin excretion (EMPEX)",
        equation="emp_ext -> emp_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        # [ml/min] * [l/ml] * [-] * [-] * [mmole/l]
        formula=(
            "GFR * 0.001 l_per_ml * fu_emp * (1 dimensionless + f_oat3 * EMPSEC_f) * emp_ext",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="EGEX",
        name="EG excretion (EGEX)",
        equation="eg_ext -> eg_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "GFR * 0.001 l_per_ml * EGEX_f * eg_ext",
            U.mmole_per_min,
        ),
    ),
]

# ---------------------------------------------------------------------------------------------------------------------
# Pharmacodynamics
# ---------------------------------------------------------------------------------------------------------------------
_m.species.extend(
    [
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
            Fasting plasma glucose (FPG), set by the studies. It shifts the RTG
            (SGLT2 is induced in type 2 diabetes) and, with the difference of
            the mean plasma glucose of the day to the FPG (`PG_delta`), gives
            the glucose which is filtered.
            """,
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
            port=True,
        ),
    ]
)

_m.parameters.extend(
    [
        Parameter(
            "Mr_glc",
            180.16,
            U.g_per_mole,
            name="molecular weight glucose",
            sboTerm=SBO.MOLECULAR_MASS,
        ),
        Parameter(
            "cf_mg_per_g",
            1000,
            U.mg_per_g,
            name="conversion factor mg per g",
        ),
        Parameter(
            "fpg_healthy",
            5.0,
            U.mM,
            name="fasting plasma glucose (healthy)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "PG_delta",
            0.7,
            U.mM,
            name="mean plasma glucose of the day minus the fasting plasma glucose",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The glucose which is filtered over the day is the mean plasma glucose
            of the day, above the fasting plasma glucose by the meals. In healthy
            subjects the 24-hr mean plasma glucose is 5.6-5.7 mM [Seman2013,
            Sarashina2013], 0.7 mM above a fasting plasma glucose of 5 mM. Studies
            which report the 24-hr mean plasma glucose set the difference to the
            FPG.
            """,
        ),
        Parameter(
            "RTG_E50",
            **FITS.parameter("KI__RTG_E50", unit=U.mM),
            name="EC50 reduction in RTG",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            EC50 of the reduction of the RTG by plasma empagliflozin (total).
            The IC50 of empagliflozin for human SGLT2 is 3.1 nM [Grempler2012];
            SGLT2 is inhibited from the tubular lumen, where empagliflozin is
            concentrated by the filtration, the secretion and the reabsorption
            of water. The urinary glucose excretion of 1 mg is a third of the
            maximal one [Sarashina2013]; the fit gives about 4 nM total plasma
            empagliflozin (0.6 nM unbound).
            """,
        ),
        Parameter(
            "RTG_gamma",
            1,
            U.dimensionless,
            name="hill coefficient reduction in RTG",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "RTG_base",
            **FITS.parameter("KI__RTG_base", unit=U.mM),
            name="baseline RTG",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""RTG without SGLT2 inhibition in healthy subjects (FPG 5 mM),
            which corresponds to the SGLT2 capacity; about 10-11 mM (180-200
            mg/dl).""",
        ),
        Parameter(
            "RTG_m_fpg",
            **FITS.parameter("KI__RTG_m_fpg", unit=U.dimensionless),
            name="FPG effect on RTG",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""Increase of the RTG per mM of FPG above the healthy FPG:
            SGLT2 is induced in type 2 diabetes [Rahmoune2005]. Phenomenological:
            it gives the urinary glucose excretion of type 2 diabetes (78-102 g
            with 25-50 mg [Heise2013, Sarashina2014, Macha2014f]) with the
            constant glucose of the model, which has no fall of the glucose
            during the treatment.""",
        ),
        Parameter(
            "RTG_max_inhibition",
            **FITS.parameter("KI__RTG_max_inhibition", unit=U.dimensionless),
            name="maximal relative reduction of the RTG",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""Maximal reduction of the RTG by SGLT2 inhibition relative to
            the RTG without inhibition: SGLT1 reabsorbs the glucose which SGLT2
            does not.""",
        ),
        Parameter(
            "RTG_splay",
            1.0,
            U.mM,
            name="splay of the glucose excretion around the RTG",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""The thresholds of the nephrons differ, so the glucose
            excretion rises smoothly around the RTG (splay) instead of at a
            sharp threshold; the softplus with a width of 1 mM.""",
        ),
    ]
)

_m.rules.extend(
    [
        AssignmentRule(
            "pg",
            "fpg + PG_delta",
            U.mM,
            name="mean plasma glucose of the day",
            notes="""The glucose which is filtered.""",
        ),
        AssignmentRule(
            "RTG_fpg",
            "RTG_base + RTG_m_fpg * (fpg - fpg_healthy)",
            U.mM,
            name="RTG without SGLT2 inhibition",
            notes="""
            FPG dependent RTG without SGLT2 inhibition. SGLT2 is induced in T2DM
            [Rahmoune2005].
            """,
        ),
        AssignmentRule(
            "RTG_delta",
            "RTG_fpg * RTG_max_inhibition",
            U.mM,
            name="maximal reduction of the RTG",
        ),
        AssignmentRule(
            variable="RTG",
            value="RTG_fpg - RTG_delta * power(emp_ext, RTG_gamma)/ (power(RTG_E50, RTG_gamma) + power(emp_ext, RTG_gamma))",
            unit=U.mM,
            name="renal threshold glucose (RTG)",
            notes="""
            Renal threshold glucose.

            The renal threshold for glucose (RTG) is the plasma glucose
            concentration at which the tubular reabsorption of glucose begins to
            saturate; glucose is excreted into the urine in proportion to the
            glucose concentration above this threshold.
            """,
        ),
    ]
)

# Glucose excretion (UGE)
_m.reactions.extend(
    [
        Reaction(
            sid="GLCEX",
            name="glucose excretion (GLCEX)",
            equation="fpg -> glc_urine [emp_ext]",
            compartment="Vki",
            sboTerm=SBO.TRANSPORT_REACTION,
            notes="""The filtered glucose above the RTG, the softplus of the
            difference with the width `RTG_splay`: the excretion approaches
            GFR * (pg - RTG) above the RTG and zero below it.""",
            # [ml/min] * [l/ml] * [mmole/l] = [mmole/min]
            formula=(
                "GFR * 0.001 l_per_ml * RTG_splay * ln(1 dimensionless + exp((pg - RTG)/RTG_splay))",
                U.mmole_per_min,
            ),
        ),
    ]
)

_m.rules.extend(
    [
        AssignmentRule(
            "UGE",
            "glc_urine * Mr_glc/cf_mg_per_g",
            unit=U.gram,
            name="urinary glucose excretion (UGE)",
            notes="""
            Urinary glucose excretion, the cumulative amount of glucose in the
            urine.
            """,
        )
    ]
)


model_kidney = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.empagliflozin import MODEL_BASE_PATH

    # create SBML
    results = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
