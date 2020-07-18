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


def _get_residue_selection():
    return MerSubclasses(Residue)


def get_ca_distance_criterion():
    ca_distance = PointsDistanceCriterion("ca", 6.0, 0.5)
    ca_distance.set_selection(_get_residue_selection())
    return ca_distance


def get_cbx_distance_criterion():
    cbx_distance = PointsDistanceCriterion("cbx", 6.5, 0.5)
    cbx_distance.set_selection(_get_residue_selection())
    return cbx_distance


def get_rc_distance_criterion():
    return PointsDistanceCriterion("rc", 7.5, 0.5)


def get_ca_cbx_vectors_difference_criterion():
    ca_cbx_vectors_difference = DistancesDifferenceCriterion("ca", "cbx", 0.75, 0.05)
    ca_cbx_vectors_difference.set_selection(_get_residue_selection())
    return ca_cbx_vectors_difference


def get_default_protein_criterion():
    ca_distance = get_ca_distance_criterion()
    cbx_distance = get_cbx_distance_criterion()
    ca_cbx_vector_criterion = get_ca_cbx_vectors_difference_criterion()
    criterion = ContactsAlternative(
        ca_distance, ContactsConjunction(cbx_distance, ca_cbx_vector_criterion)
    )
    return criterion
