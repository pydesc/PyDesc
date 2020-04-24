import numpy
from scipy.spatial.distance import cdist

from pydesc.contacts.base import ContactCriterion


def get_mer_inds(mers):
    return numpy.array([mer.ind for mer in mers])


def get_distance_matrix(mers1, mers2, point):
    points1 = numpy.array([getattr(mer, point).vector for mer in mers1])
    points2 = numpy.array([getattr(mer, point).vector for mer in mers2])
    try:
        dist_mtx = cdist(points1, points2)
    except:
        import pdb; pdb.set_trace()
    return dist_mtx


class PointsDistanceCriterion(ContactCriterion):

    def __str__(self):
        return f"distance between {self.atom_name}"

    def __init__(self, atom_name, threshold, margin):
        super().__init__()
        self.atom_name = atom_name
        self.threshold = threshold
        self.margin = margin

    def _get_vector(self, mer):
        return getattr(mer, self.atom_name).vector

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

    def __str__(self):
        return f"distance difference between {self.point1} and {self.point2}"

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

        possible_contacts = difference_mtx <= (self.threshold + self.margin)
        points1_indexes, points2_indexes = numpy.where(possible_contacts)
        inds1, inds2 = mer_inds1[points1_indexes], mer_inds2[points2_indexes]
        matrix[inds1, inds2] = 1

        sure_contacts = difference_mtx <= (self.threshold - self.margin)
        points1_indexes, points2_indexes = numpy.where(sure_contacts)
        inds1, inds2 = mer_inds1[points1_indexes], mer_inds2[points2_indexes]
        matrix[inds1, inds2] = 2

        return matrix
