# Copyright 2020 Tymoteusz 'vdhert' Oleniecki
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
"""Basic classes storing alignments."""

from abc import ABC
from abc import abstractmethod
from collections import defaultdict

import numpy


def not_dash_if_possible(ind1, ind2):
    """Pick int if int and DASH given, int if two equal ints or raise ValueError
    otherwise."""
    if ind1 is DASH:
        return ind2
    if ind2 is DASH:
        return ind1
    if ind1 == ind2:
        return ind1
    raise ValueError("Two different integers given.")


def drop_single_mer_rows(array):
    """Return new array without rows having only single mer."""
    _, n_structures = array.shape
    n_nans = numpy.count_nonzero(array == DASH, axis=1)
    new_array = array[n_nans < (n_structures - 1)]
    return new_array


not_dash_if_possible = numpy.vectorize(not_dash_if_possible, otypes=[object])


class _Dash:
    """Singleton representing missing value in alignment."""

    def __repr__(self):
        return "<->"

    def __gt__(self, other):
        return False

    def __lt__(self, other):
        return True


DASH = _Dash()


class AbstractAlignment(ABC):
    """Abstract superclass for aligned pairs and larger sets of structures."""

    def __init__(self, structures, inds_rows):
        self.structures = tuple(structures)
        self.inds = inds_rows
        self.mer_map = {structure: {} for structure in structures}
        self._fill_mer_map()

    def __getitem__(self, item):
        inds_array = self.inds[item, :]
        alignment = type(self)(self.structures, inds_array)
        return alignment

    def __len__(self):
        return self.inds.shape[0]

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
        """Return iterator that runs over rows (returning arrays of pydesc inds)"""
        return iter(self.inds)

    def iter_columns(self):
        """Return iterator that runs over columns (returning arrays of inds)
        corresponding with structures."""
        return iter(self.inds.T)

    def get_inds_aligned_with(self, structure, inds):
        """Return inds aligned with mers of given inds from given structure.

        Args:
            structure: one of aligned structures.
            inds: sequence of inds present in both: given structure and alignment.

        Returns:
            dict: other than given structures as keys and list of inds (present in
                key structure) as values. For structures not aligned with given mer
                value will be an empty list. For inconsistent alignments this method
                will list all inds aligned with given inconsistent (aligned with more
                than one ind from other structure) mers.

        """
        aligned_rows = self._get_rows_aligned_with(structure, inds)
        generator = zip(self.structures, aligned_rows.T)
        aligned_map = {stc: inds[inds != DASH].tolist() for stc, inds in generator}
        aligned_map = {k: v for k, v in aligned_map.items() if v and (k != structure)}
        return aligned_map

    def _get_rows_aligned_with(self, structure, inds):
        row_indices = set()
        for ind in inds:
            row_indices = row_indices.union(self.mer_map[structure][ind])
        row_indices = sorted(row_indices)
        aligned_rows = self.inds[row_indices]
        return aligned_rows

    def extract_aligned_with(self, structure, inds):
        """Create new alignment narrowed down to rows aligning mers with given inds from
        given structure.

        Args:
            structure: structure aligned in this alignment.
            inds: sequence of mer indices.

        Returns:
            pair or multiple alignment.

        """
        aligned_rows = self._get_rows_aligned_with(structure, inds)
        alignment = type(self)(self.structures, aligned_rows)
        return alignment

    def prune(self):
        """Return new alignment without single mer rows."""
        new_array = drop_single_mer_rows(self.inds)
        klass = type(self)
        new_alignment = klass(self.structures, new_array)
        return new_alignment

    def sum_rows(self, other):
        """Append rows from second alignment at the end of first one and return new
        alignment.

        Second alignment must align at least the same structures as self. That means
        the same objects, not just structures loaded from the same file.

        Args:
            other: alignment aligning the same structures (or more).

        Returns:
            AbstractAlignment: instance of the same type as self with extra rows (
            PairAlignment or MultipleAlignment).

        """
        other_indices = other.get_structure_indices()
        column_inds = [other_indices[structure] for structure in self.structures]
        other_array = other.inds[:, column_inds]
        all_inds = numpy.concatenate((self.inds, other_array))
        inds = numpy.vstack({tuple(row) for row in all_inds})
        klass = type(self)
        return klass(self.structures, inds)

    def get_structure_indices(self):
        """Return dict with structures as keys and their column indices as values."""
        indices = {structure: i for i, structure in enumerate(self.structures)}
        return indices

    def get_common_structures(self, other):
        """Return set of structures common for this and given alignment."""
        self_structures = set(self.structures)
        other_structures = set(other.structures)
        common = self_structures & other_structures
        return common

    def sort(self):
        """Return new alignment with sorted rows.

        Rows are sorted first for longest aligned structure, then every other row is
        inserted, if possible, to be form continuous bridge for leftmost structure.
        In case of disjoint alignments -- shorter aligned fragments will be placed at
        at the end (in random order, possibly mixed; in such case it is best to
        separate disjoint fragments and sort them separately).

        """
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
            new_row = numpy.array(new_row)
            for ind, row in enumerate(accepted_rows):
                row = numpy.array(row)
                mask = new_row != DASH
                mask &= row != DASH
                comparison = new_row[mask] < row[mask]
                if numpy.any(comparison):
                    accepted_rows.insert(ind, new_row)
                    break
            else:
                accepted_rows.append(new_row)
        accepted_rows = numpy.array(accepted_rows)
        alignment = type(self)(self.structures, accepted_rows)
        return alignment


class PairAlignment(AbstractAlignment):
    """Alignment of two structures."""

    def __init__(self, structures, inds_rows):
        if len(structures) != 2:
            raise ValueError("Pair alignment requires exactly two structures.")
        AbstractAlignment.__init__(self, structures, inds_rows)

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
        """List of pair alignments (in this case only containing this alignment)."""
        return [self]

    def transit(self, other):
        """Calculate new alignment assuming transition based on two alignments such
        that both have single common structure.

        Assuming that this alignment aligns structures A and B and other alignment
        aligns structures B and C, create new alignment aligning A and C.

        Note that if any of given alignment was inconsistent, i.e. aligned single mer
        with two mers from other structure, single arbitrary mer will be picked.

        """
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
        """Return True if two pair alignments are not contradictory.

        Only internally consistent alignments are allowed: those aligning single mer
        with more than one mer from other structure will raise ValueError (lazy).

        """
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
        """Return pruned version of this alignment (new object)."""
        return self.prune()

    def to_columns(self):
        """Return this alignment."""
        return self


class MultipleAlignment(AbstractAlignment):
    """Alignment of more than three structures."""

    def __repr__(self):
        rows, cols = self.inds.shape
        return f"<MultipleAlignment of {rows} mers ({cols} structures)>"

    def limit_to_structures(self, *structures):
        """Return new alignment cropped to given structures.

        Args:
            *structures: any number of structures present in this alignment (objects
            equality).

        """
        if len(structures) < 2:
            msg = "At least two structure are necessary to perform this operation."
            raise ValueError(msg)
        structures_inds = self.get_structure_indices()
        column_inds = [structures_inds[structure] for structure in structures]
        inds = self.inds[:, column_inds]
        if len(column_inds) == 2:
            return PairAlignment(structures, inds)
        return MultipleAlignment(structures, inds)

    def to_joined_pairs(self):
        """Cast to series of PairAlignments and pack then to JoinedPairAlignments."""
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
        """Calculate closure and return new, closed alignment.

        Closure is possible for consistent alignments only, for others behaviours is
        not predicted.

        Examples:
            Given alignment:
            A   B   C
            1   1   -
            -   1   1
            1   -   1
            2   2   -
            -   2   2
            4   4   4
            after closure it should be:
            A   B   C
            1   1   1
            2   2   2
            4   4   4
            so this method deals with fully connected graphs as well as performs its
            own transition to fully close circles

        """
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

    def drop_empty_structures(self):
        """Return new alignment without structures that had no aligned mers in this
        one."""
        non_empty_cols = []
        for col_ind, column in enumerate(self.iter_columns()):
            if numpy.any(column != DASH):
                non_empty_cols.append(col_ind)
        structures = [self.structures[i] for i in non_empty_cols]
        inds_array = self.inds[:, non_empty_cols]
        _, n_cols = inds_array.shape
        if n_cols > 2:
            klass = type(self)
        elif n_cols == 2:
            klass = PairAlignment
        else:
            msg = "Not enough non-empty structures."
            raise ValueError(msg)
        alignment = klass(structures, inds_array)
        return alignment


class JoinedPairAlignments:
    """Container for pair alignments."""

    def __init__(self, pair_alignments):
        self._pair_alignments = pair_alignments

    @property
    def pair_alignments(self):
        """List of pair alignments."""
        return self._pair_alignments

    def limit_to_structures(self, *structures):
        """Return new container storing only alignments aligning given structures.

        Args:
            *structures: any number of structures present in this alignment (objects
            equality).

        """
        alignments = []
        for alignment in self.pair_alignments:
            if set(alignment.structures).issubset(structures):
                alignments.append(alignment)
        if len(alignments) == 1:
            return max(alignments)
        return JoinedPairAlignments(alignments)

    def to_columns(self):
        """Cast to MultipleAlignment.

        Note that this method concatenates rows from each pair alignments at the end
        of final MultipleAlignment and it might be necessary to close it.

        """
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

        alignment = MultipleAlignment(structures, array)
        return alignment

    def join(self, other):
        """Merge two sets of pair alignments (possibly changing their order)."""
        pair_alignments = set(other.pair_alignments)
        pair_alignments |= set(self.pair_alignments)
        return JoinedPairAlignments(tuple(pair_alignments))
