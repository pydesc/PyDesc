"""Contact criteria based on geometrical features of mers."""
import numpy
from scipy.spatial.distance import cdist

from pydesc.contacts.base import ContactCriterion


def get_mer_inds(mers):
    """Get list of ind of given mers.

    Args:
        mers: sequence of mers.

    Returns:
        numpy.array: of mers inds in the same order as mers in given sequence.

    """
    return numpy.array([mer.ind for mer in mers])


def get_distance_matrix(mers1, mers2, point):
    """Calculate distances between given (pseudo)atoms of two sets of mers.

    If any mer from any of given sequence does not have attribute *point*,
    this function will raise AttributeError.

    Args:
        mers1: sequence of mer instances.
        mers2: second sequence of mer instances (can be the same as mers1).
        point(str): name of point to measure distance from (e.g. 'rc'), atom or
        pseudoatom.

    Returns:
        : matrix of shape (len(mers1), len(mers2)) storing distances between
        appropriate (pseudo)atoms.

    """
    points1 = numpy.array([getattr(mer, point).vector for mer in mers1])
    points2 = numpy.array([getattr(mer, point).vector for mer in mers2])
    dist_mtx = cdist(points1, points2)
    return dist_mtx


class PointsDistanceCriterion(ContactCriterion):
    """Criterion based on distance between atoms or pseudoatoms.

    This criterion measures distance between two atoms or pseudoatoms of two mers (
    both mers have to have such (pseudo)atom) than compares it to threshold.
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

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        mer_inds1 = get_mer_inds(mers1)
        mer_inds2 = get_mer_inds(mers2)

        dist_mtx = get_distance_matrix(mers1, mers2, self.atom_name)

        possible_contacts = dist_mtx <= (self.threshold + self.margin)
        points1_indexes, points2_indexes = numpy.where(possible_contacts)
        inds1, inds2 = mer_inds1[points1_indexes], mer_inds2[points2_indexes]
        matrix[inds1, inds2] = 1

        sure_contacts = dist_mtx <= (self.threshold - self.margin)
        points1_indexes, points2_indexes = numpy.where(sure_contacts)
        inds1, inds2 = mer_inds1[points1_indexes], mer_inds2[points2_indexes]
        matrix[inds1, inds2] = 2

        return matrix


class DistancesDifferenceCriterion(ContactCriterion):
    """Criterion based on difference between distances of two different points.

    For each pair it calculates distances between two (pseudo)atoms of the same kind,
    e.g. 'ca' and 'cbx' for residues. Then it subtracts those distances and compares
    difference with threshold (taking uncertain margin into account). That implements
    more abstract idea of checking if two vectors present in mer (in example above:
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

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        mer_inds1 = get_mer_inds(mers1)
        mer_inds2 = get_mer_inds(mers2)

        dist1_mtx = get_distance_matrix(mers1, mers2, self.point1)
        dist2_mtx = get_distance_matrix(mers1, mers2, self.point2)

        difference_mtx = dist1_mtx - dist2_mtx

        possible_contacts = difference_mtx >= (self.threshold - self.margin)
        points1_indexes, points2_indexes = numpy.where(possible_contacts)
        inds1, inds2 = mer_inds1[points1_indexes], mer_inds2[points2_indexes]
        matrix[inds1, inds2] = 1

        sure_contacts = difference_mtx >= (self.threshold + self.margin)
        points1_indexes, points2_indexes = numpy.where(sure_contacts)
        inds1, inds2 = mer_inds1[points1_indexes], mer_inds2[points2_indexes]
        matrix[inds1, inds2] = 2

        return matrix
