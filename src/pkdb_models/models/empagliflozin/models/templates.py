"""Template definition."""

from pkpdsim.templates import terms_of_use
from sbmlutils.factory import (
    Creator,
    ModelUnits,
    UnitDefinition,
    Units,
)

__all__ = ["U", "creators", "model_units", "terms_of_use"]


class U(Units):
    """UnitDefinitions."""

    mmole = UnitDefinition("mmole")
    min = UnitDefinition("min")
    mg = UnitDefinition("mg")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    mmole_per_min_l = UnitDefinition("mmole_per_min_l", "mmole/min/l")
    l_per_min = UnitDefinition("l_per_min", "l/min")
    per_min = UnitDefinition("per_min", "1/min")
    per_mmole = UnitDefinition("per_mmole", "1/mmole")
    mg_per_s_l = UnitDefinition("mg_per_s_l", "mg/s/l")
    g_per_mole = UnitDefinition("g_per_mole", "g/mole")
    min_per_hr = UnitDefinition("min_per_hr", "min/hour")


model_units = ModelUnits(
    time=U.min,
    extent=U.mmole,
    substance=U.mmole,
    length=U.meter,
    area=U.m2,
    volume=U.liter,
)

creators = [
    Creator(
        familyName="Alejandro",
        givenName="Jennesse",
        email="jennesse_alejandro@brown.edu",
        organization="Brown University, Providence, Rohde Island",
    ),
    Creator(
        familyName="Elias",
        givenName="Michelle",
        email="michelleelias66@gmail.com",
        organization="Humboldt-University Berlin, Institute for Biology",
    ),
    Creator(
        familyName="König",
        givenName="Matthias",
        email="koenigmx@hu-berlin.de",
        organization="Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver; University of Stuttgart, Institute of Structural Mechanics and Dynamics in Aerospace Engineering",
        site="https://livermetabolism.com",
    ),
]
