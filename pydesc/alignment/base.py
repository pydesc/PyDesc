"""Basic classes storing alignments."""

from abc import ABC
from abc import abstractmethod
from collections import defaultdict

import numpy


class _Dash:
    def __repr__(self):
        return "<->"

    def __gt__(self, other):
        return False

    def __lt__(self, other):
        return True


DASH = _Dash()


def not_dash_if_possible(ind1, ind2):
    if ind1 is DASH:
        return ind2
    if ind2 is DASH:
        return ind1
    if ind1 == ind2:
        return ind1
    raise ValueError("Two different integers given.")


not_dash_if_possible = numpy.vectorize(not_dash_if_possible, otypes=[object])


def drop_single_mer_rows(array):
    _, n_structures = array.shape
    n_nans = numpy.count_nonzero(array == DASH, axis=1)
    new_array = array[n_nans < (n_structures - 1)]
    return new_array


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
        self.mer_map = {structure: {} for structure in structures}
        self._fill_mer_map()

    @abstractmethod
    def to_joined_pairs(self):
        pass

    def _fill_mer_map(self):
        mer_map = {structure: defaultdict(list) for structure in self.structures}
        for no, row in enumerate(self.inds):
            for structure, ind in zip(self.structures, row):
                if ind is DASH:
                    continue
                mer_map[structure][ind].append(no)
        self.mer_map = {structure: {} for structure in self.structures}
        for structure in self.structures:
            stc_mer_map = mer_map[structure]
            for ind, occurs in stc_mer_map.items():
                self.mer_map[structure][ind] = numpy.array(occurs, dtype=numpy.uint32)

    def iter_rows(self):
        return iter(self.inds)

    def iter_columns(self):
        return iter(self.inds.T)

    def get_inds_aligned_with(self, structure, ind):
        row_indices = self.mer_map[structure][ind]
        aligned_rows = self.inds[row_indices]
        generator = zip(self.structures, aligned_rows.T)
        aligned_map = {stc: inds[inds != DASH].tolist() for stc, inds in generator}
        aligned_map = {k: v for k, v in aligned_map.items() if v and (k != structure)}
        return aligned_map

    def prune(self):
        new_array = drop_single_mer_rows(self.inds)
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

    def get_common_structures(self, other):
        self_structures = set(self.structures)
        other_structures = set(other.structures)
        common = self_structures & other_structures
        return common

    def sort(self):
        longest_structure = max(self.mer_map, key=lambda stc: len(self.mer_map[stc]))
        stc_col_dct = self.get_structure_indices()
        column = stc_col_dct[longest_structure]
        new_row_order = self.inds[:, column].argsort()
        partially_sorted_rows = self.inds[new_row_order, :]
        rows_to_move = numpy.where(partially_sorted_rows[:, column] == DASH)
        rows_to_push = partially_sorted_rows[rows_to_move].tolist()
        no_dash = partially_sorted_rows[:, column] != DASH
        accepted_rows = partially_sorted_rows[no_dash].tolist()
        for new_row in rows_to_push:
            for ind, row in enumerate(accepted_rows):
                comparison = new_row < row
                if numpy.any(comparison):
                    accepted_rows.insert(ind, new_row)
                    break
            else:
                accepted_rows.append(new_row)
        accepted_rows = numpy.array(accepted_rows)
        alignment = type(self)(self.structures, accepted_rows)
        return alignment


class PairAlignment(AbstractColumnAlignment, AbstractJoinedPairAlignments):
    def __init__(self, structures, inds_rows):
        if len(structures) != 2:
            raise ValueError("Pair alignment requires exactly two structures.")
        AbstractColumnAlignment.__init__(self, structures, inds_rows)

    def __eq__(self, other):
        if self is other:
            return True
        if set(self.structures) != set(other.structures):
            return False
        order = [0, 1]
        if self.structures != other.structures:
            order = [1, 0]
        if self.inds.shape != other.inds.shape:
            return False
        comparison = self.inds == other.inds[:, order]
        return numpy.all(comparison)

    def __hash__(self):
        return hash(self.structures)

    @property
    def pair_alignments(self):
        return [self]

    def limit_to_structures(self, *structures):
        if set(structures) != set(self.structures):
            msg = (
                "Pair alignment limit can only be trivial."
                "Make sure passed structures are already in it."
            )
            raise ValueError(msg)
        return self

    def transit(self, other):
        # note that if some mers are aligned multiple times,
        # arbitrary one will be picked
        try:
            (common_stc,) = self.get_common_structures(other)
        except ValueError:
            msg = (
                "To perform transition, both alignments have to have exactly one "
                "common structure."
            )
            raise ValueError(msg)
        own_map = self.mer_map[common_stc]
        own_mer_inds = set(own_map)
        other_map = other.mer_map[common_stc]
        other_mer_inds = set(other_map)
        common_mer_inds = sorted(own_mer_inds & other_mer_inds)
        own_rows = [own_map[i].max() for i in common_mer_inds]
        other_rows = [other_map[i].max() for i in common_mer_inds]
        # tricky way to get index of structure other than common
        # it works for pair alignments. column index is simply
        # int of bool value of predicate "is common structure at column 0?"
        own_2nd_structure_col = int(self.structures[0] == common_stc)
        other_2nd_structure_col = int(other.structures[0] == common_stc)

        inds_array = numpy.empty((len(common_mer_inds), 2), dtype=object)
        inds_array[:, 0] = self.inds[own_rows, own_2nd_structure_col]
        inds_array[:, 1] = other.inds[other_rows, other_2nd_structure_col]
        inds_array = drop_single_mer_rows(inds_array)

        structures = (
            self.structures[own_2nd_structure_col],
            other.structures[other_2nd_structure_col],
        )
        transit_alignment = PairAlignment(structures, inds_array)
        return transit_alignment

    def is_consistent_with(self, other):
        common_structures = self.get_common_structures(other)
        if common_structures != set(self.structures):
            return True
        for structure in common_structures:
            self_inds = set(self.mer_map[structure])
            other_inds = set(other.mer_map[structure])
            common_inds = self_inds & other_inds
            for ind in common_inds:
                self_rows = self.mer_map[structure][ind]
                self_aligned = self.inds[self_rows]
                other_rows = other.mer_map[structure][ind]
                other_aligned = other.inds[other_rows]
                try:
                    if numpy.any(self_aligned != other_aligned, axis=1):
                        return False
                except ValueError:
                    msg = (
                        "One of pair alignments is internally inconsistent."
                        "Make sure each mer is aligned with only one other mer."
                    )
                    raise ValueError(msg)
        return True

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
            for stc2 in self.structures[i + 1 :]:
                structures = stc1, stc2
                columns_inds = [structure_indices[stc] for stc in structures]
                array = self.inds[:, columns_inds]
                array = drop_single_mer_rows(array)
                alignment = PairAlignment((stc1, stc2), array)
                alignments.append(alignment)
        alignment = JoinedPairAlignments(alignments)
        return alignment

    def close(self):
        # unexpected behaviour for inconsistent alignments
        seen_mers = {}
        new_rows = []
        current_ind = -1
        for row in self.inds:
            current_ind += 1
            new_row = row.copy()
            for i, ind in enumerate(row):
                if ind is DASH:
                    continue
                seen_index = seen_mers.get((i, ind), None)
                if seen_index is None:
                    continue
                seen_row = new_rows[seen_index]
                try:
                    new_row = not_dash_if_possible(new_row, seen_row)
                except ValueError:
                    msg = "Inconsistent alignment cannot be closed."
                    raise ValueError(msg)
            new_rows.append(new_row)
            for i, ind in enumerate(new_row):
                if ind is DASH:
                    continue
                seen_mers[(i, ind)] = current_ind
        valid_new_rows_indices = sorted(set(seen_mers.values()))
        n_rows = len(valid_new_rows_indices)
        n_cols = len(self.structures)
        new_array = numpy.empty((n_rows, n_cols), dtype=object)
        for no, index in enumerate(valid_new_rows_indices):
            new_array[no] = new_rows[index]
        alignment = type(self)(self.structures, new_array)
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
