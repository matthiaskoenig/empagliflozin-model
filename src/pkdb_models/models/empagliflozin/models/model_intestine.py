"""Empagliflozin intestine model."""
import numpy as np
from sbmlutils.converters import odefac

from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.empagliflozin.models import annotations
from pkdb_models.models.empagliflozin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "empagliflozin_intestine",
    name="Model for empagliflozin absorption in the small intestine",
    notes="""
    # Model for empagliflozin absorption

    - absorption empagliflozin (emp), ~60% fraction absorbed
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=annotations.model + [
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:45615"),  # gut
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "NCIT:C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine

        (BQB.HAS_PROPERTY, "NCIT:C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
    ]
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
        "Ventero",
        1.0,
        name="intestinal lining (enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
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
        "Vbaso",
        np.nan,
        name="basolateral membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["basolateral"],
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
        f"emp_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name=f"empagliflozin (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["emp"],
        boundaryCondition=True,
    ),
    Species(
        "emp_lumen",
        initialConcentration=0.0,
        name="empagliflozin (intestinal volume)",
        compartment="Vlumen",
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
        "eg_lumen",
        initialConcentration=0.0,
        name="empagliflozin glucuronide (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
        port=True,
    ),
    Species(
        "eg_feces",
        initialConcentration=0.0,
        name="empagliflozin glucuronide (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["eg"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        f"F_emp_abs",
        0.65,
        U.dimensionless,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"fraction absorbed empagliflozin",
        notes="""
        Fraction absorbed, i.e., only a fraction of the empagliflozin in the intestinal lumen
        is absorbed. This parameter determines how much of the empagliflozin is excreted.
        
        `F_emp_abs` of dose is absorbed. `(1-F_emp_abs)` is excreted in feces.
        """,
    ),
    Parameter(
        "EMPABS_k",
        0.01,
        unit=U.per_min,
        name="rate of empagliflozin absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
    Parameter(
        "METEXC_k",
        0.001,
        unit=U.per_min,
        name="rate of empagliflozin fecal excretion",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
]

_m.rules.append(
    AssignmentRule(
        "absorption",
        value="EMPABS_k * Vgu * emp_lumen",
        unit=U.mmole_per_min,
        name="absorption empagliflozin",
    ),
)

_m.reactions = [
    Reaction(
        "EMPABS",
        name="absorption empagliflozin",
        equation="emp_lumen -> emp_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("F_emp_abs * absorption", U.mmole_per_min),
    ),
    Reaction(
        sid=f"EMPEXC",
        name=f"excretion empagliflozin (feces)",
        compartment="Vlumen",
        equation=f"emp_lumen -> emp_intestine_0",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"(1 dimensionless - F_emp_abs) * absorption",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid=f"EGEXC",
        name=f"excretion empagliflozin glucuronide (feces)",
        compartment="Vlumen",
        equation=f"eg_lumen -> eg_intestine_0",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"METEXC_k * Vgu * eg_lumen",
            U.mmole_per_min,
        ),
    )
]

substance_info = {
    "emp": "empagliflozin",
    "eg": "empagliflozin glucuronide",
}


n_chain = 5
for k in range(n_chain):
    _m.compartments.append(
        Compartment(
            f"Vintestine_{k}",
            1.2825 * 0.9 / n_chain,  # 0.0171 [l/kg] * 75 kg * 0.9,
            name="intestinal lumen",
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            unit=U.liter,
            constant=True,
            annotations=annotations.compartments["gu_lumen"],
        )
    )
    for sid, name in substance_info.items():
        _m.species.append(
            Species(
                f"{sid}_intestine_{k}",
                initialConcentration=0.0,
                name=f"{name} (intestine)",
                compartment=f"Vintestine_{k}",
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            )
        )

        if k < n_chain - 1:
            product = f"{sid}_intestine_{k + 1}"
        else:
            product = f"{sid}_feces"
        _m.reactions.append(
            Reaction(
                sid=f"{sid.upper()}EXC_{k}",
                name=f"excretion {name} {k}",
                compartment=f"Vintestine_{k}",
                equation=f"{sid}_intestine_{k} -> {product}",
                sboTerm=SBO.TRANSPORT_REACTION,
                formula=(
                    f"METEXC_k * Vintestine_{k} * {sid}_intestine_{k}",
                    U.mmole_per_min,
                )
            )
        )


_m.parameters.extend([
    Parameter(
        f"PODOSE_emp",
        0,
        U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose empagliflozin [mg]",
        port=True,
    ),
    Parameter(
        f"Ka_dis_emp",
        2.0,
        U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"Ka_dis [1/hr] dissolution empagliflozin",
        port=True
    ),
    Parameter(
        f"Mr_emp",
        450.909,
        U.g_per_mole,
        constant=True,
        name=f"Molecular weight empagliflozin [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
        port=True,
    ),
])

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid=f"dissolution_emp",
            name=f"dissolution empagliflozin",
            formula=(
                f"Ka_dis_emp/60 min_per_hr * PODOSE_emp/Mr_emp",
                U.mmole_per_min,
            ),
            equation=f"emp_stomach -> emp_lumen",
            compartment="Vgu",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """
        ),
    ]
)
_m.rate_rules.append(
    RateRule(f"PODOSE_emp", f"-dissolution_emp * Mr_emp", U.mg_per_min),
)

model_intestine = _m


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import MODEL_BASE_PATH

    results = create_model(
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        model=model_intestine, sbml_level=3, sbml_version=2
    )
    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    visualize_sbml(results.sbml_path, delete_session=False)
