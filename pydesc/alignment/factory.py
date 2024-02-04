import numpy

from pydesc.alignment.base import DASH
from pydesc.alignment.base import Alignment


def get_column_alignment_class(array):
    """Given array of aligned inds, return class that best suits it or raise
    ValueError for invalid array."""
    _, n_structures = array.shape
    if n_structures == 2:
        return PairAlignment
    elif n_structures > 2:
        return MultipleAlignment
    raise ValueError("Not enough columns to create alignment.")


class AlignmentFactory:
    def create_from_list_of_inds(self, structures, inds_lists):
        array = numpy.array(tuple(zip(*inds_lists)), dtype=object)
        array[array == None] = DASH
        return self.create_from_array(structures, array)

    def create_from_array(self, structures, array):
        return Alignment(structures, array)
