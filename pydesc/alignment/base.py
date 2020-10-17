"""Basic classes storing alignments."""

from abc import ABC
from abc import abstractmethod

import numpy

DASH = object()


class AbstractAlignment(ABC):

    @abstractmethod
    def limit_to_structures(self, *structures):
        pass


class AbstractJoinedPairAlignments(AbstractAlignment):

    @property
    @abstractmethod
    def pair_alignments(self):
        pass

    @abstractmethod
    def to_columns(self):
        pass

    def join(self, other):
        pair_alignments = set(other.pair_alignments)
        pair_alignments |= set(self.pair_alignments)
        return JoinedPairAlignments(tuple(pair_alignments))


class AbstractColumnAlignment(AbstractAlignment):
    def __init__(self, structures, inds_rows):
        self.structures = tuple(structures)
        self.inds = inds_rows
        self._mer_map = {structure: {} for structure in structures}

    @abstractmethod
    def to_joined_pairs(self):
        pass

    def _fill_mer_map(self):
        for no, row in enumerate(self.inds):
            for structure, ind in zip(self.structures, row):
                self._mer_map[structure][ind] = no

    def get_mers_aligned_with(self, mer, structure):
        array_index = self._mer_map[structure][mer.ind]
        inds = self.inds[array_index]
        generator = zip(self.structures, inds)
        mers = {structure: structure[ind] for structure, ind in generator}
        return mers

    def prune(self):
        _, n_structures = self.inds.shape
        n_nans = numpy.count_nonzero(self.inds == DASH, axis=1)
        new_array = self.inds[n_nans < (n_structures - 1)]
        klass = type(self)
        new_alignment = klass(self.structures, new_array)
        return new_alignment

    def concatenate(self, other):
        # adding rows (more aligned mers)
        other_indices = other.get_structure_indices()
        column_inds = [other_indices[structure] for structure in self.structures]
        other_array = other.inds[:, column_inds]
        all_inds = numpy.concatenate((self.inds, other_array))
        inds = numpy.unique(all_inds, axis=0)
        return MultipleColumnsAlignment(self.structures, inds)

    def get_structure_indices(self):
        indices = {structure: i for i, structure in enumerate(self.structures)}
        return indices


class PairAlignment(AbstractColumnAlignment, AbstractJoinedPairAlignments):
    def __init__(self, structures, inds_rows):
        if len(structures) != 2:
            raise ValueError("Pair alignment requires exactly two structures.")
        super().__init__(structures, inds_rows)

    def __eq__(self, other):
        if self is other:
            return True
        if set(self.structures) != set(other.structures):
            return False
        order = [0, 1]
        if self.structures != other.structures:
            order = [1, 0]
        comparison = self.inds == other.inds[:, order]
        return numpy.all(comparison)

    def __hash__(self):
        return hash(self.structures)

    @property
    def pair_alignments(self):
        return [self]

    def limit_to_structures(self, *structures):
        if set(structures) != set(self.structures):
            msg = "Pair alignment limit can only be trivial." \
                  "Make sure passed structures are already in it."
            raise ValueError(msg)
        return self

    def transit(self, alignment):
        pass

    def to_joined_pairs(self):
        return self

    def to_columns(self):
        return self


class MultipleColumnsAlignment(AbstractColumnAlignment):

    def limit_to_structures(self, *structures):
        if len(structures) < 2:
            msg = "At least two structure are necessary to perform this operation."
            raise ValueError(msg)
        structures_inds = self.get_structure_indices()
        column_inds = [structures_inds[structure] for structure in structures]
        inds = self.inds[:, column_inds]
        if len(column_inds) == 2:
            return PairAlignment(structures, inds)
        return MultipleColumnsAlignment(structures, inds)

    def to_joined_pairs(self):
        structure_indices = self.get_structure_indices()
        alignments = []
        for i, stc1 in enumerate(self.structures):
            for stc2 in self.structures[i + 1:]:
                structures = stc1, stc2
                columns_inds = [structure_indices[stc] for stc in structures]
                array = self.inds[:, columns_inds]
                alignment = PairAlignment((stc1, stc2), array)
                alignment = alignment.prune()
                alignments.append(alignment)
        alignment = JoinedPairAlignments(alignments)
        return alignment


class JoinedPairAlignments(AbstractJoinedPairAlignments):
    def __init__(self, pair_alignments):
        self._pair_alignments = pair_alignments

    @property
    def pair_alignments(self):
        return self._pair_alignments

    def limit_to_structures(self, *structures):
        alignments = []
        for alignment in self.pair_alignments:
            if set(alignment.structures).issubset(structures):
                alignments.append(alignment)
        if len(alignments) == 1:
            return max(alignments)
        return JoinedPairAlignments(alignments)

    def to_columns(self):
        structures = set()
        length = 0
        for alignment in self.pair_alignments:
            structures |= set(alignment.structures)
            length += alignment.inds.shape[0]
        structures = tuple(structures)

        structure_indices = {structure: i for i, structure in enumerate(structures)}
        array = numpy.full((length, len(structures)), DASH)

        start = 0
        for alignment in self.pair_alignments:
            columns_indices = [
                structure_indices[structure] for structure in alignment.structures
            ]
            end = start + alignment.inds.shape[0]
            array[start:end, columns_indices] = alignment.inds
            start = end

        alignment = MultipleColumnsAlignment(structures, array)
        return alignment
