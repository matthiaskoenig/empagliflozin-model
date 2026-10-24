"""Annotations of the empagliflozin models: the model and its species.

The compartments are the shared ones of `pkpdsim.annotations`.
"""

from pkpdsim.annotations import compartments
from sbmlutils.metadata import (
    BQB,
)

__all__ = ["compartments", "model", "species"]

model = [
    # taxonomy
    (BQB.HAS_TAXON, "taxonomy/9606"),  # human
    (BQB.HAS_TAXON, "snomedct/337915000"),  # human
    # modelling approach
    (BQB.HAS_PROPERTY, "mamo/MAMO_0000046"),  # ordinary differential equation model
    # biological process explained by model (GO/NCIT)
    (BQB.HAS_PROPERTY, "NCIT:C158136"),  # empagliflozin
    (BQB.HAS_PROPERTY, "NCIT:C98083"),  # SGLT2 inhibitor
    (BQB.HAS_PROPERTY, "NCIT:C15299"),  # pharmacokinetics
    # model relevance to a particular area
    (BQB.HAS_PROPERTY, "NCIT:C2985"),  # Diabetes Mellitus
    (BQB.HAS_PROPERTY, "NCIT:C99532"),  # Diabetes Therapy
    # model reference
    (BQB.IS_DESCRIBED_BY, "doi/10.5281/zenodo.18480677"),
]

species = {
    "glc": [  # glucose
        (BQB.IS, "chebi/CHEBI:4167"),
        (BQB.IS, "ncit/C2831"),
    ],
    "emp": [  # empagliflozin
        (BQB.IS, "chebi/CHEBI:82720"),
        (BQB.IS, "pubchem.compound/11949646"),
        (BQB.IS, "inchikey/OBWASQILIWPZMG-QZMOQZSNSA-N"),
        (BQB.IS, "ncit/C158136"),
        (BQB.IS, "snomedct/703894008"),
    ],
    "eg": [  # empagliflozin-glucuronide
    ],
    "emptot": [],
}
