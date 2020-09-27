"""Basic classes storing alignments."""

from abc import ABCMeta


class Alignment(metaclass=ABCMeta):
    pass


class AbstractColumnAlignment(Alignment):

    def __init__(self, structures, inds_rows):
        self.structures = structures
        self.inds = inds_rows


class PairAlignment(AbstractColumnAlignment):

    def __init__(self, structures, inds_rows):
        if len(structures) != 2:
            raise ValueError("Pair alignment requires exactly two structures.")
        super().__init__(structures, inds_rows)


class MultipleColumnsAlignment(AbstractColumnAlignment):
    pass


class JoinedPairAlignments(Alignment):

    def __init__(self, pair_alignments):
        self.pair_alignments = pair_alignments
