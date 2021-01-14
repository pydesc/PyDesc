import numpy
from scipy.linalg import get_blas_funcs

norm = get_blas_funcs("nrm2")


def get_distance_vector(mers1, mers2, point):
    points1 = get_points_array(mers1, point)
    points2 = get_points_array(mers2, point)
    vector = calculate_distances(points1, points2)
    return vector


def get_points_array(mers, point):
    return numpy.array([mer.get_point(point).vector for mer in mers])


def calculate_distances(points1, points2):
    diffs = points1 - points2
    distances = numpy.apply_along_axis(norm, 1, diffs)
    return distances


class RMSCLDCalculator:
    def __init__(
        self,
        structure1,
        structure2,
        cmap1,
        cmap2,
        uncertain_contact_weight=0.5,
        point="gc",
    ):
        self.conformation1 = structure1
        self.conformation2 = structure2
        self.certain_contacts = []
        self.uncertain_contacts = []
        self.set_contacts_union(cmap1, cmap2)
        self.uncertain_contacts_weight = uncertain_contact_weight
        self.point_name = point

    def set_contacts_union(self, cmap1, cmap2):
        merged_map = cmap1.combine(cmap2)
        certain_contacts = set()
        uncertain_contacts = set()
        for pair, value in merged_map:
            if value == 2:
                certain_contacts.add(frozenset(pair))
            else:
                uncertain_contacts.add(frozenset(pair))
        if not certain_contacts and not uncertain_contacts:
            raise ValueError("Combination of given contact maps cannot be empty.")
        self.certain_contacts = sorted(certain_contacts)
        self.uncertain_contacts = sorted(uncertain_contacts)

    def calculate(self):
        certain_contacts_sum = self._calculate_cld_sum(self.certain_contacts)
        uncertain_contacts_sum = self._calculate_cld_sum(self.uncertain_contacts)
        uncertain_contacts_sum *= self.uncertain_contacts_weight
        all_sum = certain_contacts_sum + uncertain_contacts_sum
        no_contacts = len(self.certain_contacts) + len(self.uncertain_contacts)
        score = numpy.sqrt(all_sum / no_contacts)
        return score

    def _calculate_cld_sum(self, contacts):
        if not contacts:
            return 0.0
        difference_matrix = self._calculate_cld(contacts)
        sqrt_sum = numpy.sum(difference_matrix * difference_matrix)
        return sqrt_sum

    def _calculate_cld(self, contacts_set):
        inds1, inds2 = zip(*contacts_set)
        mers1 = [self.conformation1[ind] for ind in inds1]
        mers2 = [self.conformation1[ind] for ind in inds2]
        distance_vector = get_distance_vector(mers1, mers2, self.point_name)
        mers1 = [self.conformation2[ind] for ind in inds1]
        mers2 = [self.conformation2[ind] for ind in inds2]
        distance_matrix2 = get_distance_vector(mers1, mers2, self.point_name)
        return distance_matrix2 - distance_vector
