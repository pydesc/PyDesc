# Copyright 2017 Tymoteusz Oleniecki
#
# This file is part of PyDesc.
#
# PyDesc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyDesc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.
"""Set of pre-defined ready-to-use contact criteria."""

from pydesc.contacts.base import ContactsAlternative
from pydesc.contacts.base import ContactsConjunction
from pydesc.contacts.geometrical import DistancesDifferenceCriterion
from pydesc.contacts.geometrical import PointsDistanceCriterion
from pydesc.mers.full_atom import Residue
from pydesc.selection import MerSubclasses

RESIDUE_SELECTION = MerSubclasses(Residue)

CA_DISTANCE = PointsDistanceCriterion("ca", 6.0, 0.5)
CA_DISTANCE.set_selections(RESIDUE_SELECTION, RESIDUE_SELECTION)
CBX_DISTANCE = PointsDistanceCriterion("cbx", 6.5, 0.5)
CBX_DISTANCE.set_selections(RESIDUE_SELECTION, RESIDUE_SELECTION)
RC_DISTANCE = PointsDistanceCriterion("rc", 7.5, 0.5)

CA_CBX_DISTANCE_DIFFERENCE = DistancesDifferenceCriterion("ca", "cbx", 0.75, 0.05)
CA_CBX_DISTANCE_DIFFERENCE.set_selections(RESIDUE_SELECTION, RESIDUE_SELECTION)

DEFAULT_PROTEIN = ContactsAlternative(
    CA_DISTANCE, ContactsConjunction(CBX_DISTANCE, CA_CBX_DISTANCE_DIFFERENCE)
)
