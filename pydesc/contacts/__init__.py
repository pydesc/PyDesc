"""Module storing all sort of objects responsible for contact maps calculation.

Submodules:
    base: base for contact criteria and auxiliary classes like combined criteria.
    geometrical: contact criteria based on some geometrical features of sets of atoms.
    maps: contact and frequency maps and classes related with their calculations.

    criteria: ready-to-use criteria.

"""

from pydesc.contacts.base import *
from pydesc.contacts.maps import ContactMapCalculator
