import numpy

from pydesc.alignment.base import DASH
from pydesc.alignment.loaders import get_column_alignment_class


class AlignmentFactory:
    def create_from_list_of_inds(self, structures, inds_lists):
        array = numpy.array(tuple(zip(*inds_lists)), dtype=object)
        array[array == None] = DASH
        klass = get_column_alignment_class(array)
        return klass(structures, array)
