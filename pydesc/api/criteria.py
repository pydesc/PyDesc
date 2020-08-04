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

from pydesc.chemistry.full_atom import FullAtomMer
from pydesc.chemistry.full_atom import Residue
from pydesc.contacts.base import ContactsAlternative
from pydesc.contacts.base import ContactsConjunction
from pydesc.contacts.geometrical import DistancesDifferenceCriterion
from pydesc.contacts.geometrical import PointsDistanceCriterion
from pydesc.selection import AtomSetSubclass


def _get_residue_selection():
    return AtomSetSubclass(Residue)


def get_ca_distance_criterion(threshold=6.0, margin=0.5):
    """Get alpha carbon distance criterion.

    Args:
        threshold: approximate distance at which residues CAs are no longer
            considered in contact.
        margin: margin of tolerance.

    Returns:
        : contact criterion.

    """
    ca_distance = PointsDistanceCriterion("ca", threshold, margin)
    ca_distance.set_selection(_get_residue_selection())
    return ca_distance


def get_cbx_distance_criterion(threshold=6.5, margin=0.5):
    """Get default CBX distance criterion.

    Args:
        threshold: approximate distance at which residues CBXs are no longer considered
            in contact. 6.5A by default.
        margin: margin of tolerance. 0.5A by default.

    Returns:
        : contact criterion.

    """
    cbx_distance = PointsDistanceCriterion("cbx", threshold, margin)
    cbx_distance.set_selection(_get_residue_selection())
    return cbx_distance


def get_rc_distance_criterion(threshold=7.5, margin=0.5):
    """Get geometrical center of side chain distance criterion.

    Applies to Mers in full-atom representation only. That can be changed with
    'set_selection' method. Requires AtomSet subclasses having 'rc' pseudoatom.

    Args:
        threshold: approximate distance at which mers geometrical centers are no longer
            considered in contact. 7.5 by default.
        margin: margin of tolerance. 0.5 by default.

    Returns:
        : contact criterion.

    """
    rc_distance = PointsDistanceCriterion("rc", threshold, margin)
    rc_distance.set_selection(AtomSetSubclass(FullAtomMer))
    return rc_distance


def get_gc_distance_criterion(threshold=8.5, margin=0.5):
    """Get default geometrical center distance criterion.

    Args:
        threshold: approximate distance at which geometrical centers of sets of
            atoms are no longer considered in contact. 8.5 by default.
        margin: margin of tolerance. 0.5 by default.

    Returns:
        : contact criterion.

    """
    return PointsDistanceCriterion("gc", threshold, margin)


def get_ca_cbx_vectors_difference_criterion(threshold=0.75, margin=0.05):
    """Get criterion satisfied when two mers cbx vectors point at each other.

    Args:
        threshold: approximate difference in distances between CAs and CBXs at which
        contact criterion is satisfied. 0.75A by default.
        margin: uncertain margin. 0.05A by default.

    Returns:
        : contact criterion.

    """
    ca_cbx_vectors_difference = DistancesDifferenceCriterion(
        "ca", "cbx", threshold, margin
    )
    ca_cbx_vectors_difference.set_selection(_get_residue_selection())
    return ca_cbx_vectors_difference


def get_default_protein_criterion(
    ca_threshold=6.0,
    ca_margin=0.5,
    cbx_threshold=6.5,
    cbx_margin=0.5,
    ca_cbx_threshold=0.75,
    ca_cbx_margin=0.05,
):
    """Get default contact criterion for proteins.

    Its defined as alternative of CA distance OR (conjunction of CBX and proper
    orientation of CBXs in space).
    Second criterion is satisfied if cbx vectors points at each other to some degree.
    It is calculated as difference between two distances: between two mers CAs and CBXs.

    Args:
        ca_threshold: threshold for CA criterion. 6.0A by default.
        ca_margin: margin for CA criterion. 0.5A by default.
        cbx_threshold: threshold for CBX criterion. 6.5A by default.
        cbx_margin: margin for CBX criterion. 0.5A by default.
        ca_cbx_threshold: threshold for CA-CBX difference criterion. 0.75A by default.
        ca_cbx_margin: margin for CA-CBX difference criterion. 0.05A by default.

    Returns:
        : contact criterion.

    """
    ca_distance = get_ca_distance_criterion(ca_threshold, ca_margin)
    cbx_distance = get_cbx_distance_criterion(cbx_threshold, cbx_margin)
    ca_cbx_vector_criterion = get_ca_cbx_vectors_difference_criterion(
        ca_cbx_threshold, ca_cbx_margin
    )
    criterion = ContactsAlternative(
        ca_distance, ContactsConjunction(cbx_distance, ca_cbx_vector_criterion)
    )
    return criterion
