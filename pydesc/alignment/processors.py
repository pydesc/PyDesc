import operator
from abc import ABC
from abc import abstractmethod
from collections import Counter
from collections import defaultdict
from functools import reduce
from itertools import chain

import numpy

from pydesc.alignment.base import DASH
from pydesc.alignment.factory import AlignmentFactory


def drop_single_mer_rows(array):
    """Return new array without rows having only single mer."""
    _, n_structures = array.shape
    n_nans = numpy.count_nonzero(array == DASH, axis=1)
    new_array = array[n_nans < (n_structures - 1)]
    return new_array


class AbstractProcessor(ABC):
    def __init__(self, factory=None):
        if factory is None:
            factory = AlignmentFactory()
        self._factory = factory

    @abstractmethod
    def process(self, *args, **kwargs):
        pass


class SingleAlignmentAbstractProcessor(AbstractProcessor):
    def __init__(self, alignment):
        self._alignment = alignment
        super().__init__()

    @abstractmethod
    def process(self):
        pass


class PickStructuresMixIn:
    def _get_structures_from_inds(self, indices):
        structures = numpy.fromiter(
            iter(self._alignment.get_structures()), dtype=object
        )
        return structures[indices].tolist()


class Close(SingleAlignmentAbstractProcessor):
    def __init__(self, alignment):
        super().__init__(alignment)

    def process(self):
        self._get_inds_data()
        self._set_empty_row_set()
        self._find_mer_cliques()
        new_structures = self._get_output_structures()
        array = self._get_output_array(new_structures)
        return self._factory.create_from_array(new_structures, array)

    def _set_empty_row_set(self):
        self._row_set = set()

    def _get_inds_data(self):
        self._inds_map = self._alignment.get_inds_map()
        self._inds_table = self._alignment.get_inds_table()
        self._structures = self._alignment.get_structures()
        self._structures_counter = Counter()

    def _find_mer_cliques(self):
        self._rows_to_skip = set()
        for row in self._inds_table:
            row_as_set = set()
            for structure, ind in zip(self._structures, row):
                if ind is DASH:
                    continue
                row_as_set |= self._dfs(structure, ind)
            if not row_as_set:
                continue
            row_structure_count = Counter([stc for stc, _ in row_as_set])
            self._structures_counter |= row_structure_count
            self._row_set.add(frozenset(row_as_set))

    def _prepare_row_as_set(self, row):
        return {
            (structure, ind)
            for structure, ind in zip(self._structures, row)
            if ind is not DASH
        }

    def _dfs(self, structure, ind):
        row_as_set = set()
        for row_index in self._inds_map[structure][ind]:
            if row_index in self._rows_to_skip:
                continue
            self._rows_to_skip.add(row_index)
            row_as_set |= self._prepare_row_as_set(self._inds_table[row_index])
        for structure, ind in frozenset(row_as_set):
            row_as_set |= self._dfs(structure, ind)
        return row_as_set

    def _get_output_structures(self):
        return reduce(
            operator.add,
            [[structure] * n for structure, n in self._structures_counter.items()],
        )

    def _get_output_array(self, structures):
        output = []
        stc_column_indices = {stc: structures.index(stc) for stc in set(structures)}
        for row_as_set in self._row_set:
            row_as_dict = defaultdict(list)
            for structure, ind in sorted(row_as_set, key=lambda x: x[1]):
                row_as_dict[structure].append(ind)
            row = [
                row_as_dict[structure].pop() if row_as_dict[structure] else DASH
                for structure in structures
            ]
            output.append(row)
        return numpy.array(output)


class DropSingleMerRows(SingleAlignmentAbstractProcessor):
    def process(self):
        structures = self._alignment.get_structures()
        table = self._alignment.get_inds_table()
        clean_table = drop_single_mer_rows(table)
        return self._factory.create_from_array(structures, clean_table)


class SelectStructures(SingleAlignmentAbstractProcessor):
    """Processor returning new alignment cropped to given structures.

    Args:
        alignment: alignment to process.
        structures: sequence of structures. Only structures present in that sequence will
            be present in resulting alignment.

    """

    def __init__(self, alignment, structures):
        self._structures = structures
        super().__init__(alignment)

    def process(self):
        current_structures = self._alignment.get_structures()
        columns_to_pick = [
            structure in self._structures for structure in current_structures
        ]
        new_table = self._alignment.get_inds_table()[:, columns_to_pick]
        return self._factory.create_from_array(self._structures, new_table)


class ConvertToPairwiseAlignments(SingleAlignmentAbstractProcessor):
    def __init__(self, alignment):
        super().__init__(alignment)

    def process(self):
        """Cast to series of pairwise Alignments."""
        if len(self._alignment.get_structures()) == 2:
            return [self._alignment]
        alignments = []
        structures = self._alignment.get_structures()
        for i, stc1 in enumerate(structures):
            for stc2 in structures[i + 1 :]:
                alignment = SelectStructures(self._alignment, [stc1, stc2]).process()
                alignment = DropSingleMerRows(alignment).process()
                alignments.append(alignment)
        return alignments


class Sort(SingleAlignmentAbstractProcessor, PickStructuresMixIn):
    def process(self):
        table = self._alignment.get_inds_table()
        col_order_by_len = numpy.flip(numpy.argsort((table != DASH).sum(axis=0)))
        table = table[:, col_order_by_len]

        class Row(tuple):
            def __lt__(self, other):
                for ind1, ind2 in zip(self, other):
                    if DASH in (ind1, ind2) or ind1 == ind2:
                        continue
                    return ind1 < ind2
                return False

        table = numpy.array(sorted(table, key=Row))
        structures = self._get_structures_from_inds(col_order_by_len)
        return self._factory.create_from_array(structures, table)


class DropUnalignedStructures(SingleAlignmentAbstractProcessor, PickStructuresMixIn):
    def process(self):
        table = self._alignment.get_inds_table()
        non_empty_columns_inds = (table != DASH).sum(axis=0) > 0
        table = table[:, non_empty_columns_inds]
        structures = self._get_structures_from_inds(non_empty_columns_inds)
        return self._factory.create_from_array(structures, table)


class Reorder(SingleAlignmentAbstractProcessor, PickStructuresMixIn):
    def __init__(self, alignment, structures):
        super().__init__(alignment)
        self._structures_order = structures
        alignment_structures = alignment.get_structures()
        if any(structure not in alignment_structures for structure in structures):
            msg = "Structures passed to Reorder not present in alignment to reorder."
            raise ValueError(msg)

    def process(self):
        table = self._alignment.get_inds_table()
        structures = self._alignment.get_structures()
        structure_to_indices_map = defaultdict(list)
        for n, structure in enumerate(structures):
            structure_to_indices_map[structure].append(n)
        indices = list(
            chain(
                *(
                    structure_to_indices_map[structure]
                    for structure in self._structures_order
                )
            )
        )
        table = table[:, indices]
        structures = self._get_structures_from_inds(indices)
        return self._factory.create_from_array(structures, table)


class MultipleAlignmentsAbstractProcessor(AbstractProcessor):
    def __init__(self, alignments):
        self._alignments = alignments
        super().__init__()

    @abstractmethod
    def process(self):
        pass


class Merge(MultipleAlignmentsAbstractProcessor):
    def __init__(self, alignments):
        self._structures_counter = Counter()
        self._structure_index_map = {}
        self._structures = ()
        self._total_length = 0
        self._array = None
        super().__init__(alignments)

    def process(self):
        self._count_structures_occurs()
        self._prepare_structures()
        self._count_total_lenght()
        self._prepare_empty_array()
        self._fill_array()
        return self._factory.create_from_array(self._structures, self._array)

    def _count_structures_occurs(self):
        for alignment in self._alignments:
            self._structures_counter |= Counter(alignment.get_structures())

    def _prepare_structures(self):
        self._structures = reduce(
            operator.add,
            [[structure] * n for structure, n in self._structures_counter.items()],
        )
        self._structure_index_map = {
            structure: self._structures.index(structure)
            for structure in self._structures
        }

    def _count_total_lenght(self):
        self._total_length = reduce(
            operator.add, [len(alignment) for alignment in self._alignments]
        )

    def _prepare_empty_array(self):
        self._array = numpy.full((self._total_length, len(self._structures)), DASH)

    def _fill_array(self):
        start_row_index = 0
        for alignment in self._alignments:
            column_inds = self._get_column_indices(alignment)
            self._paste_alignment_into_array(column_inds, start_row_index, alignment)
            start_row_index += len(alignment)

    def _get_column_indices(self, alignment):
        local_counter = Counter()
        indices = []
        for structure in alignment.get_structures():
            indices.append(
                self._structure_index_map[structure] + local_counter[structure]
            )
            local_counter[structure] += 1
        return indices

    def _paste_alignment_into_array(self, column_indices, start_row_index, alignment):
        self._array[
            start_row_index : start_row_index + len(alignment), column_indices
        ] = alignment.get_inds_table()
