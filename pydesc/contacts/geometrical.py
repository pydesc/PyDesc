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
"""Contact criteria based on geometrical features of atom sets."""

import numpy
from scipy.spatial.distance import cdist

from pydesc.contacts.base import ContactCriterion
from pydesc.warnexcept import warn, Info


def get_inds(atom_sets):
    """Get list of ind of given AtomSet instances.

    Args:
        atom_sets: sequence of AtomSet instances.

    Returns:
        numpy.array: inds in order corresponding to argument.

    """
    return numpy.array([atom_set.ind for atom_set in atom_sets])


def get_distance_matrix(atom_sets1, atom_sets2, point):
    """Calculate distances between given (pseudo)atoms of two sequences of AtomSet
    instances.

    If any instance from any of given sequences does not have attribute *point*,
    this function will raise AttributeError.

    Args:
        atom_sets1: sequence of AtomSet instances.
        atom_sets2: second sequence of AtomSet instances (can be the same as
            atom_sets1).
        point(str): name of point to measure distance from (e.g. 'rc'), atom or
            pseudoatom.

    Raises:
        TypeError: when given atom sets are empty.

    Returns:
        : matrix of shape (len(atom_sets1), len(atom_sets2)) storing distances between
        appropriate (pseudo)atoms.

    """
    points1 = numpy.array([getattr(atom_set, point).vector for atom_set in atom_sets1])
    points2 = numpy.array([getattr(atom_set, point).vector for atom_set in atom_sets2])
    try:
        dist_mtx = cdist(points1, points2)
    except ValueError:
        # meaning arrays are empty
        msg = "Empty atom sets given."
        raise TypeError(msg)
    return dist_mtx


class PointsDistanceCriterion(ContactCriterion):
    """Criterion based on distance between atoms or pseudoatoms.

    This criterion measures distance between two atoms or pseudoatoms of two atom sets
    (both have to have such (pseudo)atom) than compares it to threshold.
    Everything below threshold - margin counts as sure contact (2) while distances
    between threshold - margin and threshold + margin are considered uncertain (1).

    Args:
        atom_name(str): name of (pseudo)atom.
        threshold(float): distance threshold.
        margin(float): uncertain range from threshold.

    """

    def __str__(self):
        return f"distance between {self.atom_name}s"

    def __repr__(self):
        return f"<PointsDistanceCriterion based on {self.atom_name} distances>"

    def __init__(self, atom_name, threshold, margin):
        super().__init__()
        self.atom_name = atom_name
        self.threshold = threshold
        self.margin = margin

    def _fill_contact_matrix(self, atom_sets1, atom_sets2, matrix):
        structure_ids1 = get_inds(atom_sets1)
        structure_ids2 = get_inds(atom_sets2)

        try:
            dist_mtx = get_distance_matrix(atom_sets1, atom_sets2, self.atom_name)
        except TypeError:
            msg = (
                f"Criterion '{str(self)}' does not apply to any part of given "
                f"structure. Maybe this criterion is inappropriate for the kind of "
                f"representation of given structure?"
            )
            warn(Info(msg))
            return matrix

        possible_contacts = dist_mtx <= (self.threshold + self.margin)
        points1_indexes, points2_indexes = numpy.where(possible_contacts)
        inds1, inds2 = structure_ids1[points1_indexes], structure_ids2[points2_indexes]
        matrix[inds1, inds2] = 1

        sure_contacts = dist_mtx <= (self.threshold - self.margin)
        points1_indexes, points2_indexes = numpy.where(sure_contacts)
        inds1, inds2 = structure_ids1[points1_indexes], structure_ids2[points2_indexes]
        matrix[inds1, inds2] = 2

        return matrix


class DistancesDifferenceCriterion(ContactCriterion):
    """Criterion based on difference between distances of two different points.

    For each pair it calculates distances between two (pseudo)atoms of the same kind,
    e.g. 'ca' and 'cbx' for residues. Then it subtracts those distances and compares
    difference with threshold (taking uncertain margin into account). That implements
    more abstract idea of checking if two vectors present in atom set (in example above:
    ca->cbx) are pointing more or less at each other.
    Note that during comparison actual thresholds for values 2 and 1 are calculated
    as, respectively, threshold + margin and threshold - margin.

    Args:
        point1(str): name of 1st (pseudo)atom.
        point2(str): name of 2nd (pseudo)atom.
        threshold(float): distance difference threshold.
        margin(float): distance difference uncertain margin.

    """

    def __str__(self):
        return f"distance difference between {self.point1} and {self.point2}"

    def __repr__(self):
        return (
            f"<DistanceDifferenceCriterion based on difference between "
            f"distances of {self.point1}s and {self.point2}s>"
        )

    def __init__(self, point1, point2, threshold, margin):
        super().__init__()
        self.point1 = point1
        self.point2 = point2
        self.threshold = threshold
        self.margin = margin

    def _fill_contact_matrix(self, atom_sets1, atom_sets2, matrix):
        structure_ids1 = get_inds(atom_sets1)
        structure_ids2 = get_inds(atom_sets2)

        try:
            dist1_mtx = get_distance_matrix(atom_sets1, atom_sets2, self.point1)
            dist2_mtx = get_distance_matrix(atom_sets1, atom_sets2, self.point2)
        except TypeError:
            msg = (
                f"Criterion '{str(self)}' does not apply to any part of given "
                f"structure. Maybe this criterion is inappropriate for the kind of "
                f"representation of given structure?"
            )
            warn(Info(msg))
            return matrix

        difference_mtx = dist1_mtx - dist2_mtx

        possible_contacts = difference_mtx >= (self.threshold - self.margin)
        points1_indexes, points2_indexes = numpy.where(possible_contacts)
        inds1, inds2 = structure_ids1[points1_indexes], structure_ids2[points2_indexes]
        matrix[inds1, inds2] = 1

        sure_contacts = difference_mtx >= (self.threshold + self.margin)
        points1_indexes, points2_indexes = numpy.where(sure_contacts)
        inds1, inds2 = structure_ids1[points1_indexes], structure_ids2[points2_indexes]
        matrix[inds1, inds2] = 2

        return matrix
