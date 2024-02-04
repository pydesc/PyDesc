from collections import Counter

import numpy
import pytest
from common import get_arange_array
from common import get_n_mocked_structures
from common import get_trivial_array

from pydesc.alignment.base import DASH
from pydesc.alignment.base import Alignment
from pydesc.alignment.base import IncorrectAlignmentError
from pydesc.alignment.processors import Close
from pydesc.alignment.processors import ConvertToPairwiseAlignments
from pydesc.alignment.processors import DropSingleMerRows
from pydesc.alignment.processors import DropUnalignedStructures
from pydesc.alignment.processors import Merge
from pydesc.alignment.processors import Reorder
from pydesc.alignment.processors import SelectStructures
from pydesc.alignment.processors import Sort


class TestMerge:
    def test_simple_case(self):
        # GIVEN
        stcs = get_n_mocked_structures(5)
        al1 = Alignment(stcs[:3], get_trivial_array(10, 3))
        al2 = Alignment(stcs[2:], get_trivial_array(10, 3))

        result = Merge([al1, al2]).process()

        assert set(result.get_structures()) == set(stcs)
        assert result.get_inds_table().shape == (20, 5)

    def test_merge_two_with_self_align(self):
        stc1, stc2, stc3 = get_n_mocked_structures(3)
        self_al = Alignment((stc1, stc1, stc2), get_arange_array(10, 3))
        al = Alignment((stc1, stc3), get_trivial_array(10, 2))

        result = Merge((self_al, al)).process()

        assert result.get_inds_table().shape == (20, 4)
        assert Counter(result.get_structures())[stc1] == 2

        assert len(result[(stc1, 0)]) == 2
        assert len(result[(stc1, 10)]) == 1


class TestClose:
    def test_close_trivial(self):
        arr = numpy.array(
            [
                [1, 1, DASH, DASH],
                [DASH, 1, 1, DASH],
                [DASH, DASH, 1, 1],
                [1, DASH, DASH, 1],
                [2, 2, DASH, DASH],
                [DASH, 2, 2, 2],
                [3, 3, DASH, DASH],
                [DASH, DASH, 4, 4],
            ]
        )
        structures = get_n_mocked_structures(4)
        alignment = Alignment(structures, arr)
        closed_alignment = Close(alignment).process()

        inds_table = closed_alignment.get_inds_table()
        assert inds_table.shape == (4, 4)
        rows = [set(zip(closed_alignment.get_structures(), i)) for i in inds_table]
        expected_rows = [
            set(zip(structures, [1, 1, 1, 1])),
            set(zip(structures, [2, 2, 2, 2])),
            set(zip(structures, [3, 3, DASH, DASH])),
        ]
        for expected_row in expected_rows:
            assert expected_row in rows

    def test_close_inconsistent(self, mocked_structures3):
        arr = numpy.array(
            [
                [1, 1, DASH],
                [DASH, 1, 1],
                [1, DASH, 2],
            ]
        )
        alignment = Alignment(mocked_structures3, arr)
        res = Close(alignment).process()
        assert len(res.get_structures()) == 4
        # duplicates col with one of the structures


class TestDropSingleMerRows:
    def test_basic(self):
        payload = numpy.array(
            [
                [0, DASH, 0],
                [1, 0, 1],
                [DASH, DASH, 4],
            ]
        )
        alignment = Alignment(get_n_mocked_structures(3), payload)
        alignment = DropSingleMerRows(alignment).process()

        assert len(alignment) == 2
        assert len(alignment.get_structures()) == 3


class TestSelectStructures:
    def test_limit_structures(self):
        structures = get_n_mocked_structures(4)
        array = numpy.array([[i] * 4 for i in range(5)])
        alignment = Alignment(structures, array)
        stc1, stc2, stc3, stc4 = structures

        # WHEN
        with pytest.raises(IncorrectAlignmentError):
            SelectStructures(alignment, (stc1,)).process()
        new_pa = SelectStructures(alignment, [stc1, stc2]).process()
        new_ma = SelectStructures(alignment, [stc2, stc1, stc4]).process()

        # THEN
        assert new_pa.get_inds_table().shape == (5, 2)
        assert new_ma.get_inds_table().shape == (5, 3)


class TestConvertToPairwiseAlignments:
    def test_split_to_joined_pairs(self, triple_nontrivial_alignment):
        triple_nontrivial_alignment._table[0, 0] = DASH
        pairwise_alignments = ConvertToPairwiseAlignments(
            triple_nontrivial_alignment
        ).process()

        assert len(pairwise_alignments) == 3
        shapes = set()
        for pa in pairwise_alignments:
            shapes.add(pa.get_inds_table().shape)

        assert shapes == {(9, 2), (10, 2)}

    def test_split_to_joined_pairs_on_pairwise(self, pairwise_nontrivial_alignment):
        result = ConvertToPairwiseAlignments(pairwise_nontrivial_alignment).process()
        assert pairwise_nontrivial_alignment == result[0]

    def test_columns_to_joined_pairs(self, mocked_structures3):
        payload = numpy.array([[0, 0, 0], [1, 1, 1], [DASH, 2, 2]])
        alignment = Alignment(mocked_structures3, payload)
        jp_alignment = ConvertToPairwiseAlignments(alignment).process()

        assert len(jp_alignment) == 3
        shapes = set()
        for pa in jp_alignment:
            shapes.add(pa.get_inds_table().shape)

        assert shapes == {(3, 2), (2, 2), (2, 2)}


class TestSort:
    def test_sort(self):
        col0 = numpy.array([DASH] * 5 + list(range(8)))
        col1 = numpy.array(list(range(10)) + [DASH] * 3)
        col2 = numpy.array(list(range(5)) + [DASH] * 8)
        col3 = numpy.array([1, 2, DASH, DASH, 5, 6, 7, 8, 9, 10, 11, 12, 13])
        arr = numpy.empty((13, 4), dtype=object)
        arr[:, 0] = col0
        arr[:, 1] = col1
        arr[:, 2] = col2
        arr[:, 3] = col3
        messy_arr = arr[[2, 6, 1, 12, 5, 7, 0, 11, 9, 4, 10, 3, 8], :]
        structures = get_n_mocked_structures(4)

        alignment = Alignment(structures, messy_arr)

        sorted_al = Sort(alignment).process()

        expected_col_order = [3, 1, 0, 2]
        numpy.testing.assert_array_equal(
            sorted_al.get_inds_table(), arr[:, expected_col_order]
        )

    def test_sort_disjoint(self):
        arr = numpy.array(
            [
                [1, 1, DASH, DASH],
                [DASH, DASH, 43, 43],
                [3, 3, DASH, DASH],
                [2, 2, DASH, DASH],
                [DASH, DASH, 44, 44],
                [4, 4, DASH, DASH],
            ]
        )

        expected_left = numpy.array(
            [
                [1, 1],
                [2, 2],
                [3, 3],
                [4, 4],
            ]
        )
        expected_right = numpy.array(
            [
                [43, 43],
                [44, 44],
            ]
        )

        structures = get_n_mocked_structures(4)
        alignment = Alignment(structures, arr)
        sorted_al = Sort(alignment).process()
        assert sorted_al.get_inds_table().shape == (6, 4)
        clean_left_al = DropSingleMerRows(
            SelectStructures(sorted_al, structures[:2]).process()
        ).process()
        numpy.testing.assert_array_equal(clean_left_al.get_inds_table(), expected_left)
        clean_right_al = DropSingleMerRows(
            SelectStructures(sorted_al, structures[2:]).process()
        ).process()
        numpy.testing.assert_array_equal(
            clean_right_al.get_inds_table(), expected_right
        )


class TestDropUnalignedStructures:
    def test_process(self, mocked_structures3):
        arr = numpy.array(
            [
                [1, 1, DASH],
                [2, 2, DASH],
                [3, 3, DASH],
                [4, 4, DASH],
                [5, 5, DASH],
            ]
        )
        alignment = Alignment(mocked_structures3, arr)
        result = DropUnalignedStructures(alignment).process()

        stc1, stc2 = result.get_structures()
        assert stc1 == mocked_structures3[0]
        assert stc2 == mocked_structures3[1]
        assert DASH not in result.get_inds_table()


class TestReorder:
    def test_reorder(self, mocked_structures3):
        stc1, stc2, stc3 = mocked_structures3
        arr = numpy.array(
            [
                [1, 1, 1],
                [2, DASH, 2],
                [DASH, 2, 2],
            ]
        )
        alignment = Alignment([stc1, stc2, stc3], arr)
        reordered_al = Reorder(alignment, [stc2, stc1, stc3]).process()
        new_arr = reordered_al.get_inds_table()
        col1 = new_arr[:, 0]
        col2 = new_arr[:, 1]
        numpy.testing.assert_array_equal(col1, [1, DASH, 2])
        numpy.testing.assert_array_equal(col2, [1, 2, DASH])

    def test_reorder_negative(self, triple_alignment):
        (stc4,) = get_n_mocked_structures(1)
        stc1, stc2, stc3 = triple_alignment.get_structures()
        with pytest.raises(IncorrectAlignmentError):
            Reorder(triple_alignment, [stc1]).process()
        with pytest.raises(ValueError):
            Reorder(triple_alignment, [stc4]).process()

