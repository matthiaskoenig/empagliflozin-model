"""Empagliflozin body model."""

from typing import NotRequired, TypedDict

import numpy as np
from sbmlutils.factory import (
    PORT_SUFFIX,
    AssignmentRule,
    Compartment,
    ExternalModelDefinition,
    Model,
    Parameter,
    RateRule,
    Reaction,
    ReplacedElement,
    Species,
    Submodel,
    UnitDefinition,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.empagliflozin.models import annotations, templates
from pkdb_models.models.empagliflozin.models.fits import FITS


# -----------------------------------------------------------------------------
# Whole Body Metabolism
# -----------------------------------------------------------------------------
class U(templates.U):
    """Units of the whole-body model."""

    kg = UnitDefinition("kg")
    l_per_kg = UnitDefinition("l_per_kg", "l/kg")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")
    ml_per_s = UnitDefinition("ml_per_s", "ml/s")
    ml_per_s_kg = UnitDefinition("ml_per_s_kg", "ml/s/kg")
    s_per_min = UnitDefinition("s_per_min", "s/min")


_m = Model(
    "empagliflozin_body",
    name="Empagliflozin body model",
    notes="""
    # Empagliflozin body model

    Distribution of empagliflozin and its glucuronide (EG) with the blood flow,
    the tissue distribution of empagliflozin into a fast and a deep pool, the
    hepatic and renal glucuronidation, the absorption and fecal excretion and
    the renal excretion of empagliflozin and EG.
    """
    + templates.terms_of_use,
    units=U,
    model_units=templates.model_units,
    creators=templates.creators,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "bto/BTO:0001489"),  # whole body
        (BQB.OCCURS_IN, "fma/FMA:256135"),  # body
        (BQB.HAS_PROPERTY, "NCIT:C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
        (BQB.HAS_PROPERTY, "NCIT:C79370"),  # Pharmacokinetics: Distribution
    ],
)

_m.ports = []
_m.deletions = []
_m.replaced_elements = []


class SubstanceBody(TypedDict):
    """Definition of a substance transported in the whole-body model."""

    name: str
    unit: UnitDefinition
    cinit: float  # initial concentration [mmole/l]
    Mr: float  # molecular weight [g/mole]
    IVDOSE: NotRequired[float]  # intravenous dose [mg]
    PODOSE: NotRequired[float]  # oral dose [mg]
    #: distribution in the tissues (ftissue, Kp, fitted)
    tissue_distribution: NotRequired[bool]
    #: distribution in the deep tissue of the rest of the body (ftissue_deep,
    #: Kp_deep, fitted)
    deep_tissue: NotRequired[bool]


SUBSTANCES_BODY: dict[str, SubstanceBody] = {
    "emp": {
        "name": "empagliflozin",
        "unit": U.mmole,
        # initial concentration
        "cinit": 0.0,  # [mmole/l]
        # Molecular weight
        "Mr": 450.909,  # [g/mole]
        # iv doses [mg]
        "IVDOSE": 0,  # venous plasma
        # oral absorption
        "PODOSE": 0,  # dose
        # distribution into the tissues: Vss/F 73.8 l, CL/F 10.6 l/hr, t1/2
        # 12.4 hr (label); the plasma falls after the peak with a half-life of
        # 3-4 hr, then 7-10 hr (12-48 hr) and 16-20 hr in the studies sampled to
        # 72-96 hr [Seman2013, Macha2014b, Macha2014f, Sarashina2014]; the
        # exchange flow and the tissue/plasma partition coefficient (total,
        # 86 % bound in plasma) are fitted, see `fits.FITS`
        "tissue_distribution": True,
        # a deep, slowly equilibrating pool of the rest of the body: one tissue
        # pool has one phase after the distribution, the multi-phasic decline
        # of the oral profiles with a terminal half-life of 12-20 hr needs a
        # second, slow pool (fitted)
        "deep_tissue": True,
    },
    # the glucuronide is hydrophilic, it has no tissue distribution: its
    # kinetics are limited by its formation; EG is the sum of the 2-O-, 3-O-
    # and 6-O-glucuronides [Chen2015a]
    "eg": {
        "name": "empagliflozin-glucuronide",
        "unit": U.mmole,
        # initial concentration
        "cinit": 0.0,  # [mmole/l]
        # Molecular weight
        "Mr": 627.0,  # [g/mole]
    },
}

# -----------------------------------------------------------------------------
# Submodels
# -----------------------------------------------------------------------------
COMPARTMENTS_BODY: dict[str, str] = {
    "ki": "kidney",
    "li": "liver",
    "lu": "lung",
    "gu": "gut",
    "re": "rest",
    "ar": "arterial blood",
    "ve": "venous blood",
    "po": "portal vein",
    "hv": "hepatic vein",
}

SUBMODEL_SID_DICT: dict[str, str] = {  # tissue to submodel mapping
    "ki": "KI",  # kidney
    "li": "LI",  # liver
    "gu": "GU",  # gut/intestine
}

kidney_id = "empagliflozin_kidney"
liver_id = "empagliflozin_liver"
gut_id = "empagliflozin_intestine"

_m.external_model_definitions = [
    ExternalModelDefinition(
        sid="kidney",
        name="kidney model definition",
        source=f"{kidney_id}.xml",
        modelRef=kidney_id,
    ),
    ExternalModelDefinition(
        sid="liver",
        name="liver model definition",
        source=f"{liver_id}.xml",
        modelRef=liver_id,
    ),
    ExternalModelDefinition(
        sid="gut", name="gut model definition", source=f"{gut_id}.xml", modelRef=gut_id
    ),
]

_m.submodels = [
    Submodel(sid=SUBMODEL_SID_DICT["ki"], modelRef="kidney", name="kidney submodel"),
    Submodel(sid=SUBMODEL_SID_DICT["li"], modelRef="liver", name="liver submodel"),
    Submodel(sid=SUBMODEL_SID_DICT["gu"], modelRef="gut", name="gut submodel"),
]

# -------------------------------------------------------------------------------------------------
# Compartments
# -------------------------------------------------------------------------------------------------
_m.compartments = [
    Compartment(
        "Vre",
        metaId="meta_Vre",
        value=1,
        unit=U.liter,
        constant=False,
        name="rest of body",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["re"],
    ),
    Compartment(
        "Vgu",
        metaId="meta_Vgu",
        value=1,
        unit=U.liter,
        constant=False,
        name="gut",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["gu"],
    ),
    Compartment(
        "Vki",
        metaId="meta_Vki",
        value=1,
        unit=U.liter,
        constant=False,
        name="kidney",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
    ),
    Compartment(
        "Vli",
        metaId="meta_Vli",
        value=1,
        unit=U.liter,
        constant=False,
        name="liver",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
    ),
    Compartment(
        "Vlu",
        metaId="meta_Vlu",
        value=1,
        unit=U.liter,
        constant=False,
        name="lung",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["lu"],
    ),
    Compartment(
        "Vve",
        metaId="meta_Vve",
        value=1,
        unit=U.liter,
        constant=False,
        name="venous blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ve"],
    ),
    Compartment(
        "Var",
        metaId="meta_Var",
        value=1,
        unit=U.liter,
        constant=False,
        name="arterial blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ar"],
    ),
    Compartment(
        "Vurine",
        metaId="meta_Vurine",
        value=1,
        unit=U.liter,
        constant=True,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["urine"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["stomach"],
    ),
    Compartment(
        "Vlumen",
        metaId="meta_Vlumen",
        value=1,
        unit=U.liter,
        constant=False,
        name="intestinal lumen (inner part of intestine)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["gu_lumen"],
    ),
    Compartment(
        "Vpo",
        metaId="meta_Vpo",
        value=1,
        unit=U.liter,
        constant=False,
        name="portal plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["po"],
    ),
    Compartment(
        "Vhv",
        metaId="meta_Vhv",
        value=1,
        unit=U.liter,
        constant=False,
        name="hepatic venous plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["hv"],
    ),
]

# create plasma and tissue compartments (for correct blood volume)
for cid in COMPARTMENTS_BODY:
    if cid not in ["ar", "ve", "po", "hv"]:
        _m.compartments.extend(
            [
                Compartment(
                    f"V{cid}_tissue",
                    value=1,
                    unit=U.liter,
                    constant=False,
                    name=f"{cid} tissue",
                    metaId=f"meta_V{cid}_tissue",
                    port=True,
                    annotations=annotations.compartments["parenchyma"]
                    + [
                        (BQB.IS_PART_OF, aid[1])
                        for aid in annotations.compartments[cid]
                    ],
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                ),
                Compartment(
                    f"V{cid}_plasma",
                    value=1,
                    unit=U.liter,
                    constant=False,
                    name=f"{cid} plasma",
                    metaId=f"meta_V{cid}_plasma",
                    port=True,
                    annotations=annotations.compartments["plasma"]
                    + [
                        (BQB.IS_PART_OF, aid[1])
                        for aid in annotations.compartments[cid]
                    ],
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                ),
            ]
        )

# replace volumes of submodels
_m.replaced_elements.extend(
    [
        # kidney
        ReplacedElement(
            sid="Vki_tissue_RE",
            metaId="Vki_tissue_RE",
            elementRef="Vki_tissue",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vki{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vki_plasma_RE",
            metaId="Vki_plasma_RE",
            elementRef="Vki_plasma",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vurine_RE",
            metaId="Vurine_RE",
            elementRef="Vurine",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vurine{PORT_SUFFIX}",
        ),
        # liver
        ReplacedElement(
            sid="Vli_tissue_RE",
            metaId="Vli_tissue_RE",
            elementRef="Vli_tissue",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"Vli{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vli_plasma_RE",
            metaId="Vli_plasma_RE",
            elementRef="Vli_plasma",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        # intestine
        ReplacedElement(
            sid="Vgu_plasma_RE",
            metaId="Vgu_plasma_RE",
            elementRef="Vgu_plasma",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vgu_RE",
            metaId="Vgu_RE",
            elementRef="Vgu",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vgu{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vfeces_RE",
            metaId="Vfeces_RE",
            elementRef="Vfeces",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vfeces{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vstomach_RE",
            metaId="Vstomach_RE",
            elementRef="Vstomach",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vstomach{PORT_SUFFIX}",
        ),
        # the lumen of the intestine
        ReplacedElement(
            sid="Vlumen_RE",
            metaId="Vlumen_RE",
            elementRef="Vlumen",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vlumen{PORT_SUFFIX}",
        ),
    ]
)

# -------------------------------------------------------------------------------------------------
# Species
# -------------------------------------------------------------------------------------------------
_m.species = []
for sid, sdict in SUBSTANCES_BODY.items():
    for cid, cname in COMPARTMENTS_BODY.items():
        if cid in ["ve", "ar", "po", "hv"]:
            sid_ex = f"C{cid}_{sid}"
            cid_ex = f"V{cid}"
        else:
            # plasma compartment
            sid_ex = f"C{cid}_plasma_{sid}"
            cid_ex = f"V{cid}_plasma"

        _m.species.append(
            Species(
                sid_ex,
                metaId=f"meta_{sid_ex}",
                initialConcentration=sdict["cinit"],
                compartment=cid_ex,
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} ({cname} plasma)",
                hasOnlySubstanceUnits=False,
                port=True,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            )
        )

    # tissue species
    for cid, cname in COMPARTMENTS_BODY.items():
        if sdict.get("tissue_distribution") and cid not in [
            "ve",
            "ar",
            "li",
            "ki",
            "gu",
            "po",
            "hv",
        ]:
            sid_ex = f"C{cid}_{sid}"
            cid_ex = f"V{cid}_tissue"

            _m.species.extend(
                [
                    Species(
                        sid_ex,
                        metaId=f"meta_{sid_ex}",
                        initialConcentration=sdict["cinit"],
                        compartment=cid_ex,
                        substanceUnit=sdict["unit"],
                        name=f"{sdict['name']} ({cname})",
                        hasOnlySubstanceUnits=False,
                        port=True,
                        sboTerm=SBO.SIMPLE_CHEMICAL,
                        annotations=annotations.species[sid],
                    )
                ]
            )

    # deep tissue of the rest of the body
    if sdict.get("deep_tissue"):
        _m.species.append(
            Species(
                f"Cre_deep_{sid}",
                metaId=f"meta_Cre_deep_{sid}",
                initialConcentration=sdict["cinit"],
                compartment="Vre_tissue",
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} (rest, deep tissue)",
                hasOnlySubstanceUnits=False,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            )
        )

    # urine
    _m.species.append(
        Species(
            f"Aurine_{sid}",
            metaId=f"meta_Aurine_{sid}",
            initialConcentration=0,
            compartment="Vurine",
            substanceUnit=sdict["unit"],
            name=f"{sdict['name']} (urine)",
            hasOnlySubstanceUnits=True,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[sid],
        )
    )

# intestinal lumen and feces: empagliflozin; Chen2015a found no glucuronide in
# the feces, EG has no biliary excretion
_m.species.extend(
    [
        Species(
            "Cgu_emp",
            metaId="meta_Cgu_emp",
            initialConcentration=0,
            compartment="Vlumen",
            substanceUnit=U.mmole,
            name="empagliflozin (intestinal lumen)",
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["emp"],
        ),
        Species(
            "Afeces_emp",
            metaId="meta_Afeces_emp",
            initialConcentration=0,
            compartment="Vfeces",
            substanceUnit=U.mmole,
            name="empagliflozin (feces)",
            hasOnlySubstanceUnits=True,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["emp"],
        ),
        Species(
            "Cve_emptot",
            metaId="meta_Cve_emptot",
            initialConcentration=0,
            compartment="Vve",
            substanceUnit=U.mmole,
            name="total empagliflozin (plasma)",
            notes="""Empagliflozin and EG in the plasma, compared with the total
            radioactivity of the plasma.""",
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["emptot"],
        ),
        Species(
            "Aurine_emptot",
            metaId="meta_Aurine_emptot",
            initialConcentration=0,
            compartment="Vurine",
            substanceUnit=U.mmole,
            name="total empagliflozin (urine)",
            notes="""Empagliflozin and EG in the urine, compared with the total
            radioactivity of the urine.""",
            hasOnlySubstanceUnits=True,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["emptot"],
        ),
        Species(
            "Afeces_emptot",
            metaId="meta_Afeces_emptot",
            initialConcentration=0,
            compartment="Vfeces",
            substanceUnit=U.mmole,
            name="total empagliflozin (feces)",
            notes="""Empagliflozin in the feces, compared with the total
            radioactivity of the feces (the oxidative metabolites of the feces,
            3.6 % of the dose [Chen2015a], are not modelled).""",
            hasOnlySubstanceUnits=True,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["emptot"],
        ),
    ]
)

# replace species
replaced_species = {
    "li": ["emp", "eg"],
    "ki": ["emp", "eg"],
    "gu": ["emp"],
}

# blood to external species
for tkey, skey_list in replaced_species.items():
    for skey in skey_list:
        _m.replaced_elements.extend(
            [
                ReplacedElement(
                    sid=f"C{tkey}_plasma_{skey}_RE",
                    metaId=f"C{tkey}_plasma_{skey}_RE",
                    elementRef=f"C{tkey}_plasma_{skey}",
                    submodelRef=SUBMODEL_SID_DICT[tkey],
                    portRef=f"{skey}_ext{PORT_SUFFIX}",
                ),
            ]
        )

# urine species replacements
for skey in ["emp", "eg", "emptot"]:
    _m.replaced_elements.append(
        ReplacedElement(
            sid=f"Aurine_{skey}_RE",
            metaId=f"Aurine_{skey}_RE",
            elementRef=f"Aurine_{skey}",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"{skey}_urine{PORT_SUFFIX}",
        )
    )

# feces and intestinal lumen species replacements
_m.replaced_elements.append(
    ReplacedElement(
        sid="Afeces_emp_RE",
        metaId="Afeces_emp_RE",
        elementRef="Afeces_emp",
        submodelRef=SUBMODEL_SID_DICT["gu"],
        portRef=f"emp_feces{PORT_SUFFIX}",
    )
)
_m.replaced_elements.append(
    ReplacedElement(
        sid="Cgu_emp_RE",
        metaId="Cgu_emp_RE",
        elementRef="Cgu_emp",
        submodelRef=SUBMODEL_SID_DICT["gu"],
        portRef=f"emp_lumen{PORT_SUFFIX}",
    )
)

# -------------------------------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------------------------------
_m.parameters.extend(
    [
        # whole body data
        Parameter(
            "BW",
            75,
            U.kg,
            constant=True,
            name="body weight",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "BW_ref",
            75,
            U.kg,
            constant=True,
            name="reference body weight [kg]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=(
                "Body weight at which the exchange flows of the tissue "
                "distribution (`ftissue_<sid>`) are given; they scale with the "
                "body weight as the volumes, the blood flows and the renal "
                "clearance do."
            ),
        ),
        Parameter(
            "COBW",
            1.548,
            U.ml_per_s_kg,
            constant=True,
            name="cardiac output per bodyweight",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_cardiac_function",
            1.0,
            U.dimensionless,
            constant=False,
            name="heart function",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            1.0: normal function;
            <1.0: reduced function as reduced cardiac output;
            """,
        ),
        Parameter(
            "CO",
            108.33,
            U.ml_per_s,
            constant=False,
            name="cardiac output [ml/s]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "QC",
            108.33 / 1000 * 60,
            U.l_per_min,
            constant=False,
            name="cardiac output [l/min]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # fractional tissue volumes
        Parameter(
            "Fblood",
            0.02,
            U.dimensionless,
            constant=False,
            name="blood fraction of organ volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "HCT",
            0.51,
            U.dimensionless,
            constant=True,
            metaId="meta_HCT",
            name="hematocrit",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # fractional tissue volumes
        Parameter(
            "FVgu",
            0.0171,
            U.l_per_kg,
            constant=True,
            name="gut fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVki",
            0.0044,
            U.l_per_kg,
            constant=True,
            name="kidney fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVli",
            0.0210,
            U.l_per_kg,
            constant=True,
            name="liver fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVlu",
            0.0076,
            U.l_per_kg,
            constant=True,
            name="lung fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVre",
            0,
            U.l_per_kg,
            constant=False,
            name="rest of body fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # calculated based on other tissues
        Parameter(
            "FVve",
            0.0514,
            U.l_per_kg,
            constant=True,
            name="venous fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVar",
            0.0257,
            U.l_per_kg,
            constant=True,
            name="arterial fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVpo",
            0.001,
            U.l_per_kg,
            constant=True,
            name="portal fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVhv",
            0.001,
            U.l_per_kg,
            constant=True,
            name="hepatic venous fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # fractional tissue blood flows
        Parameter(
            "FQgu",
            0.18,
            U.dimensionless,
            constant=True,
            name="gut fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
        Assumed FQgu is the sum of FQgu + FQsp + FQpa for correct portal blood flow.
        with FQgu=0.146, FQsp=0.017, FQpa=0.017
        """,
        ),
        Parameter(
            "FQki",
            0.190,
            U.dimensionless,
            constant=True,
            name="kidney fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQh",
            0.215,
            U.dimensionless,
            constant=True,
            name="hepatic (venous side) fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQlu",
            1,
            U.dimensionless,
            constant=True,
            name="lung fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQre",
            0,
            U.dimensionless,
            constant=False,
            name="rest of body fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_cirrhosis",
            0,
            U.dimensionless,
            constant=True,
            name="severity of cirrhosis [0, 0.95]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            annotations=[
                (BQB.IS, "NCIT:C2951"),  # Cirrhosis
                (BQB.IS, "efo/0001422"),  # cirrhosis of liver
            ],
        ),
        Parameter(
            "f_shunts",
            0,
            U.dimensionless,
            constant=False,
            name="fraction of portal venous blood shunted by the liver",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            annotations=[
                (
                    BQB.IS,
                    "hp/HP:0002629",
                ),  # Gastrointestinal arteriovenous malformation
            ],
        ),
        Parameter(
            "f_tissue_loss",
            0,
            U.dimensionless,
            constant=False,
            name="fraction of lost parenchymal liver volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    ]
)


# renal function, the glucuronidation of liver and kidney and its activity
_m.parameters.extend(
    [
        Parameter(
            "f_GFR_Qki",
            101 / 1.32354,
            U.ml_per_l,
            constant=True,
            name="relationship GFR and Qki",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            GFR per renal blood flow (76.3 ml/l, a ratio of 0.0763). Calibrated
            so that the healthy GFR is 101 ml/min at the renal blood flow of
            1.32354 l/min of 75 kg, the normal GFR of the classes of renal
            impairment (`pkpdsim.impairment`).
            """,
        ),
        Parameter(
            "GFR_base",
            np.nan,
            U.ml_per_min,
            constant=False,
            name="glomerular filtration rate (healthy)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The GFR of healthy kidneys, which scales with the renal blood flow
            and so with the body weight; the GFR is `KI__f_renal_function`
            times this rate.
            """,
        ),
        Parameter(
            "f_ugt",
            1.0,
            U.dimensionless,
            constant=True,
            name="scaling factor UGT activity",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Scaling factor of the glucuronidation of empagliflozin (UGT2B7,
            UGT1A3, UGT1A8, UGT1A9) of the liver and the kidney; an inducer or
            inhibitor changes both.
            1.0: normal activity; <1.0 decreased activity; >1.0 increased activity.
            """,
        ),
        Parameter(
            "EMP2EG_k",
            **FITS.parameter("EMP2EG_k", unit=U.per_min),
            constant=True,
            name="rate of the hepatic glucuronidation of empagliflozin",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            First-order rate of the hepatic clearance of empagliflozin by uptake
            (OATP1B1/1B3) and glucuronidation, relative to the volume of the
            liver tissue (`LI__EMP2EG_k` of the liver). Empagliflozin is far
            below the Km of the UGTs, its pharmacokinetics are dose proportional
            from 0.5 to 200 mg [Seman2013]. The renal glucuronidation is
            `f_ugt_renal` of the glucuronidation of both organs.
            """,
        ),
        Parameter(
            "f_ugt_renal",
            0.30,
            U.dimensionless,
            constant=True,
            name="renal fraction of the glucuronidation",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Fraction of the glucuronidation of empagliflozin in the kidney at
            normal renal and hepatic function, the rest is hepatic. The split is
            not identified by the data of healthy subjects; it is the one which
            gives the exposure in hepatic and in renal impairment together: the
            AUC ratios of severe hepatic impairment (Child-Pugh C, 1.75
            [Macha2014b]) and of severe renal impairment (eGFR 15-29, 1.67
            [Macha2014f]) are reached with 30 % of the glucuronidation in the
            kidney, which scales with the renal function.
            """,
        ),
        Parameter(
            "EMP2EG_k_renal",
            np.nan,
            U.per_min,
            constant=False,
            name="rate of the renal glucuronidation of empagliflozin",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            First-order rate of the renal glucuronidation relative to the volume
            of the kidney (`KI__EMP2EG_k` of the kidney), from the hepatic rate
            and the renal fraction `f_ugt_renal`.
            """,
        ),
    ]
)
_m.replaced_elements.extend(
    [
        ReplacedElement(
            sid="GFR_base_RE",
            metaId="GFR_base_RE",
            elementRef="GFR_base",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"GFR_base{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="f_ugt_ki_RE",
            metaId="f_ugt_ki_RE",
            elementRef="f_ugt",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"f_ugt{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="f_ugt_li_RE",
            metaId="f_ugt_li_RE",
            elementRef="f_ugt",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"f_ugt{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="EMP2EG_k_li_RE",
            metaId="EMP2EG_k_li_RE",
            elementRef="EMP2EG_k",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"EMP2EG_k{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="EMP2EG_k_ki_RE",
            metaId="EMP2EG_k_ki_RE",
            elementRef="EMP2EG_k_renal",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"EMP2EG_k{PORT_SUFFIX}",
        ),
    ]
)

for sid, sdict in SUBSTANCES_BODY.items():
    if "PODOSE" in sdict:
        # Parameter for intestine model
        _m.parameters.extend(
            [
                Parameter(
                    f"PODOSE_{sid}",
                    sdict["PODOSE"],
                    U.mg,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"oral dose {sid} [mg]",
                ),
            ]
        )
        _m.replaced_elements.extend(
            [
                ReplacedElement(
                    sid=f"PODOSE_{sid}_RE",
                    metaId=f"PODOSE_{sid}_RE",
                    elementRef=f"PODOSE_{sid}",
                    submodelRef=SUBMODEL_SID_DICT["gu"],
                    portRef=f"PODOSE_{sid}{PORT_SUFFIX}",
                ),
            ]
        )

_m.rules = [
    *_m.rules,
    AssignmentRule("f_shunts", "f_cirrhosis", unit=U.dimensionless),
    AssignmentRule("f_tissue_loss", "f_cirrhosis", unit=U.dimensionless),
]

replaced_parameters: dict[str, list[str]] = {
    "li": [],
    "gu": ["Mr_emp"],
    "ki": [],
}

for ckey, pids in replaced_parameters.items():
    for pid in pids:
        _m.replaced_elements.append(
            ReplacedElement(
                sid=f"{pid}_{ckey}_RE",
                metaId=f"{pid}_{ckey}_RE",
                elementRef=f"{pid}",
                submodelRef=SUBMODEL_SID_DICT[ckey],
                portRef=f"{pid}{PORT_SUFFIX}",
            )
        )

# species specific parameters
for sid, sdict in SUBSTANCES_BODY.items():
    _m.parameters.extend(
        [
            # molecular weights
            Parameter(
                f"Mr_{sid}",
                sdict["Mr"],
                U.g_per_mole,
                constant=True,
                name=f"Molecular weight {sid} [g/mole]",
                sboTerm=SBO.MOLECULAR_MASS,
            ),
        ]
    )

    # tissue distribution
    if sdict.get("tissue_distribution"):
        _m.parameters.extend(
            [
                Parameter(
                    f"ftissue_{sid}",
                    **FITS.parameter(f"ftissue_{sid}", unit=U.l_per_min),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"tissue distribution {sid}",
                    notes="Exchange flow at the reference body weight `BW_ref`.",
                ),
                Parameter(
                    f"Kp_{sid}",
                    **FITS.parameter(f"Kp_{sid}", unit=U.dimensionless),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"tissue/plasma partition coefficient {sid}",
                ),
            ]
        )
    if sdict.get("deep_tissue"):
        _m.parameters.extend(
            [
                Parameter(
                    f"ftissue_deep_{sid}",
                    **FITS.parameter(f"ftissue_deep_{sid}", unit=U.l_per_min),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"deep tissue distribution {sid}",
                    notes=(
                        "Exchange flow of the deep tissue pool of the rest of the "
                        "body at the reference body weight `BW_ref`."
                    ),
                ),
                Parameter(
                    f"Kp_deep_{sid}",
                    **FITS.parameter(f"Kp_deep_{sid}", unit=U.dimensionless),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"deep tissue/plasma partition coefficient {sid}",
                ),
            ]
        )

    # dosing
    if "IVDOSE" in sdict:
        _m.parameters.extend(
            [
                Parameter(
                    f"IVDOSE_{sid}",
                    sdict["IVDOSE"],
                    U.mg,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"IV bolus dose {sid} [mg]",
                ),
                # iv kinetics after application
                Parameter(
                    f"ti_{sid}",
                    10,
                    U.second,
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"injection time {sid} [s]",
                ),
                Parameter(
                    f"Ki_{sid}",
                    0.02,
                    U.per_min,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"Ki [1/min] injection {sid}",
                ),
                # continuous infusion
                Parameter(
                    f"Ri_{sid}",
                    0,
                    U.mg_per_min,
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"Ri [mg/min] rate of infusion {sid}",
                ),
                Parameter(
                    f"cum_dose_{sid}",
                    0,
                    U.mg,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"Cumulative dose due to infusion {sid}",
                ),
            ]
        )

# -------------------------------------------------------------------------------------------------
# AssignmentRules
# -------------------------------------------------------------------------------------------------
_m.rules = [
    *_m.rules,
    AssignmentRule(
        "Cve_emptot",
        "Cve_emp + Cve_eg",
        unit=U.mM,
        name="sum of empagliflozin and EG (plasma)",
    ),
    AssignmentRule(
        "Afeces_emptot",
        "Afeces_emp",
        unit=U.mmole,
        name="total empagliflozin (feces)",
    ),
    # the renal glucuronidation is f_ugt_renal of the glucuronidation of liver
    # and kidney (at normal renal and hepatic function)
    AssignmentRule(
        "EMP2EG_k_renal",
        "EMP2EG_k * FVli/FVki * f_ugt_renal/(1 dimensionless - f_ugt_renal)",
        U.per_min,
        name="rate of the renal glucuronidation of empagliflozin",
    ),
    # Rest body volume
    AssignmentRule(
        "FVre",
        "1.0 l_per_kg - (FVgu + FVki + FVli + FVlu + FVve + FVar)",
        U.l_per_kg,
    ),
    # Rest body perfusion
    AssignmentRule("FQre", "1.0 dimensionless - (FQki + FQh)", U.dimensionless),
    # cardiac output (depending on heart rate and bodyweight)
    AssignmentRule("CO", "f_cardiac_function*BW*COBW", U.ml_per_s),
    # cardiac output (depending on bodyweight)
    AssignmentRule("QC", "CO/1000 ml_per_l * 60 s_per_min", U.l_per_min),
    # volumes
    AssignmentRule("Vgu", "BW*FVgu", U.liter),
    AssignmentRule("Vki", "BW*FVki", U.liter),
    AssignmentRule("Vli", "BW*FVli", U.liter),
    AssignmentRule("Vlu", "BW*FVlu", U.liter),
    AssignmentRule("Vre", "BW*FVre", U.liter),
    # venous and arterial blood volume (corrected for tissue blood volumes)
    AssignmentRule(
        "Vve",
        "BW*FVve - FVve/(FVar+FVve) * BW * Fblood * (1 l_per_kg - FVve - FVar)",
        U.liter,
    ),
    AssignmentRule(
        "Var",
        "BW*FVar - FVar/(FVar+FVve) * BW * Fblood * (1 l_per_kg - FVve - FVar)",
        U.liter,
    ),
    AssignmentRule(
        "Vpo",
        "(1 dimensionless - HCT) * (BW*FVpo - FVpo/(FVar+FVve+FVpo+FVhv) "
        "* BW * Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv)))",
        U.liter,
    ),
    AssignmentRule(
        "Vhv",
        "(1 dimensionless - HCT) * (BW*FVhv - FVhv/(FVar+FVve+FVpo+FVhv) "
        "* BW * Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv)))",
        U.liter,
    ),
    # blood flows
    AssignmentRule("Qgu", "QC*FQgu", U.l_per_min, name="gut blood flow"),
    AssignmentRule("Qki", "QC*FQki", U.l_per_min, name="kidney blood flow"),
    AssignmentRule(
        "Qh", "QC*FQh", U.l_per_min, name="hepatic (venous side) blood flow"
    ),
    AssignmentRule("Qha", "Qh - Qgu", U.l_per_min, name="hepatic artery blood flow"),
    AssignmentRule("Qlu", "QC*FQlu", U.l_per_min, name="lung blood flow"),
    AssignmentRule("Qre", "QC*FQre", U.l_per_min, name="rest of body blood flow"),
    AssignmentRule("Qpo", "Qgu", U.l_per_min, name="portal blood flow"),
    # glomerular filtration rate of healthy kidneys
    AssignmentRule(
        "GFR_base",
        "f_GFR_Qki * Qki",
        U.ml_per_min,
        name="glomerular filtration rate (healthy)",
    ),
    # the lumen is the inner part of the intestine
    AssignmentRule(
        "Vlumen",
        "0.9 dimensionless * Vgu",
        U.liter,
        name="volume of the intestinal lumen",
    ),
]

# Volumes for explicit tissue models
for cid, cname in COMPARTMENTS_BODY.items():
    if cid not in ["ve", "ar", "po", "hv"]:
        _m.rules.append(
            # plasma volume associated with tissue
            AssignmentRule(
                f"V{cid}_plasma",
                value=f"V{cid} * Fblood * (1 dimensionless - HCT)",
                unit=U.liter,
                name=f"plasma volume of {cname}",
            ),
        )

        if cid == "li":
            _m.rules.append(
                # Cirrhosis: Adjustment of fractional liver tissue volume
                AssignmentRule(
                    f"V{cid}_tissue",
                    value=f"V{cid}*(1 dimensionless - f_tissue_loss) * (1 dimensionless - Fblood)",
                    unit=U.liter,
                    name=f"tissue volume of {cname}",
                ),
            )
        else:
            _m.rules.append(
                AssignmentRule(
                    f"V{cid}_tissue",
                    value=f"V{cid}*(1 dimensionless - Fblood)",
                    unit=U.liter,
                    name=f"tissue volume of {cname}",
                ),
            )

for sid, sdict in SUBSTANCES_BODY.items():
    # injection
    if "IVDOSE" in sdict:
        _m.rules.extend(
            [
                AssignmentRule(
                    f"Ki_{sid}",
                    f"0.693 dimensionless/ti_{sid} * 60 s_per_min",
                    U.per_min,
                    name="injection rate IV",
                    notes="""
                log(2) ~ 0.693 used for conversion of half-life injection time to
                injection rate.
                """,
                ),
            ]
        )

# --------------------------------------------------------------------------------------------------
# Reactions
# --------------------------------------------------------------------------------------------------
_m.reactions = []

for sid, sdict in SUBSTANCES_BODY.items():
    sname = sdict["name"]

    # --------------------
    # tissue distribution
    # --------------------
    if sdict.get("tissue_distribution"):
        for cid in COMPARTMENTS_BODY:
            if cid not in ["ve", "ar", "li", "gu", "ki", "po", "hv"]:
                _m.reactions.append(
                    Reaction(
                        sid=f"transport_{cid}_{sid}",
                        name=f"transport {sname}",
                        formula=(
                            f"ftissue_{sid} * BW/BW_ref * (C{cid}_plasma_{sid}*Kp_{sid} - C{cid}_{sid})",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_plasma_{sid} <-> C{cid}_{sid}",
                        compartment=f"V{cid}_tissue",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                )

    if sdict.get("deep_tissue"):
        _m.reactions.append(
            Reaction(
                sid=f"transport_re_deep_{sid}",
                name=f"transport {sname} (deep tissue)",
                formula=(
                    f"ftissue_deep_{sid} * BW/BW_ref * "
                    f"(Cre_plasma_{sid}*Kp_deep_{sid} - Cre_deep_{sid})",
                    U.mmole_per_min,
                ),
                equation=f"Cre_plasma_{sid} <-> Cre_deep_{sid}",
                compartment="Vre_tissue",
                sboTerm=SBO.TRANSPORT_REACTION,
            ),
        )

    if "IVDOSE" in sdict:
        _m.reactions.extend(
            [
                # --------------------
                # iv dose application
                # --------------------
                Reaction(
                    sid=f"iv_{sid}",
                    name=f"iv {sname}",
                    formula=(f"Ki_{sid}*IVDOSE_{sid}/Mr_{sid}", U.mmole_per_min),
                    equation=f"-> Cve_{sid}",
                    sboTerm=SBO.TRANSPORT_REACTION,
                    annotations=[(BQB.IS, "NCIT:C38276")],
                    compartment="Vve",
                ),
            ]
        )

    for cid, cname in COMPARTMENTS_BODY.items():
        if cid in ["ve", "ar", "po", "hv"]:
            continue
        # --------------------
        # ve -> lung -> ar
        # --------------------
        if cid == "lu":
            rid_in = f"Flow_ve_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"

            rid_out = f"Flow_{cid}_ar_{sid}"
            name_out = f"outflow {cname} {sname}"

            _m.reactions.extend(
                [
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Cve_{sid}", U.mmole_per_min),
                        equation=f"Cve_{sid} -> C{cid}_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_plasma_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_plasma_{sid} -> Car_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> organ -> ve
        # --------------------
        if cid in ["ki", "re"]:
            rid_in = f"Flow_ar_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"
            rid_out = f"Flow_{cid}_ve_{sid}"
            name_out = f"outflow {cname} {sname}"

            # only distribution in plasma volume
            _m.reactions.extend(
                [
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> C{cid}_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_plasma_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_plasma_{sid} -> Cve_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> organ -> po
        # --------------------
        if cid in ["gu"]:
            rid_in = f"Flow_ar_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"
            rid_out = f"Flow_{cid}_po_{sid}"
            name_out = f"outflow {cname} {sname}"

            _m.reactions.extend(
                [
                    # only distribution in plasma volume
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> C{cid}_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_plasma_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_plasma_{sid} -> Cpo_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> li
        # po -> li
        # li -> hv -> ve
        # --------------------
        if cid == "li":
            # liver
            _m.reactions.extend(
                [
                    Reaction(
                        sid=f"Flow_arli_li_{sid}",
                        name=f"arterial inflow liver {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*Qha*Car_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Car_{sid} -> Cli_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # shunted arterial flow
                    Reaction(
                        sid=f"Flow_arli_hv_{sid}",
                        name="flow arterial shunts",
                        formula=(f"f_shunts*Qha*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # (unshunted) portal flow
                    Reaction(
                        sid=f"Flow_po_li_{sid}",
                        name=f"outflow po {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*Qpo*Cpo_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Cpo_{sid} -> Cli_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # portal shunts
                    Reaction(
                        sid=f"Flow_po_hv_{sid}",
                        name="flow portal shunts",
                        formula=(f"f_shunts*Qpo*Cpo_{sid}", U.mmole_per_min),
                        equation=f"Cpo_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=f"Flow_li_hv_{sid}",
                        name=f"outflow liver {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*(Qpo+Qha)*Cli_plasma_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Cli_plasma_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=f"Flow_hv_ve_{sid}",
                        name=f"outflow hepatic vein {sname}",
                        formula=(f"Qh*Chv_{sid}", U.mmole_per_min),
                        equation=f"Chv_{sid} -> Cve_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )

# --------------------------------------------------------------------------------------------------
# RateRules
# --------------------------------------------------------------------------------------------------
_m.rate_rules = []

for sid, sdict in SUBSTANCES_BODY.items():
    if "IVDOSE" in sdict:
        _m.rate_rules.extend(
            [
                # injection of dose
                RateRule(
                    f"IVDOSE_{sid}", f"-iv_{sid}*Mr_{sid} + Ri_{sid}", U.mg_per_min
                ),
                # cumulative infusion dose
                RateRule(f"cum_dose_{sid}", f"Ri_{sid}", U.mg_per_min),
            ]
        )

model_body: Model = _m


if __name__ == "__main__":
    from pkdb_models.models.empagliflozin import MODEL_BASE_PATH

    result = create_model(
        filepath=MODEL_BASE_PATH / f"{model_body.sid}.xml",
        model=model_body,
        sbml_level=3,
        sbml_version=2,
    )
