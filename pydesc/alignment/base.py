"""Basic classes storing alignments."""

from abc import ABCMeta

import numpy

DASH = object()


class Alignment(metaclass=ABCMeta):
    pass


class AbstractColumnAlignment(Alignment):
    def __init__(self, structures, inds_rows):
        self.structures = tuple(structures)
        self.inds = inds_rows

    def simplify(self):
        _, n_structures = self.inds.shape
        n_nans = numpy.count_nonzero(self.inds == DASH, axis=1)
        self.inds = self.inds[n_nans < (n_structures - 1)]


class PairAlignment(AbstractColumnAlignment):
    def __init__(self, structures, inds_rows):
        if len(structures) != 2:
            raise ValueError("Pair alignment requires exactly two structures.")
        super().__init__(structures, inds_rows)

    def transit(self, alignment):
        pass


class MultipleColumnsAlignment(AbstractColumnAlignment):
    pass


class JoinedPairAlignments(Alignment):
    def __init__(self, pair_alignments):
        self.pair_alignments = pair_alignments
