from unittest.mock import MagicMock

import numpy
import pytest

from pydesc.alignment.base import JoinedPairAlignments
from pydesc.alignment.base import MultipleColumnsAlignment
from pydesc.alignment.base import PairAlignment
from pydesc.alignment.loaders import DASH


def get_n_mocked_structures(n, *, start=0):
    return [MagicMock(name=f"stc{i}") for i in range(start, start + n)]


def get_trivial_pair_alignment(stc1, stc2):
    arr = numpy.array([[i, i] for i in range(10)])
    return PairAlignment((stc1, stc2), arr)


def get_3_pair_alignments():
    stc1, stc2, stc3 = get_n_mocked_structures(3)
    pa1 = get_trivial_pair_alignment(stc1, stc2)
    pa2 = get_trivial_pair_alignment(stc1, stc3)
    pa3 = get_trivial_pair_alignment(stc2, stc3)
    return pa1, pa2, pa3


@pytest.fixture
def mocked_structures3():
    return get_n_mocked_structures(3)


class TestColumnAlignment:
    def test_prune(self):
        payload = numpy.array([[0, DASH, 0], [1, 0, 1], [DASH, DASH, 4],])
        alignment = MultipleColumnsAlignment([None, None, None], payload)
        alignment = alignment.prune()

        assert alignment.inds.shape == (2, 3)

    @pytest.fixture
    def trivial_array3(self):
        array = numpy.array([[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]])
        return array

    def test_columns_to_joined_pairs(self, mocked_structures3):
        payload = numpy.array([[0, 0, 0], [1, 1, 1], [DASH, 2, 2]])
        alignment = MultipleColumnsAlignment(mocked_structures3, payload)
        jp_alignment = alignment.to_joined_pairs()

        assert len(jp_alignment.pair_alignments) == 3
        shapes = set()
        for pa in jp_alignment.pair_alignments:
            shapes.add(pa.inds.shape)

        assert shapes == {(3, 2), (2, 2), (2, 2)}

    def test_concatenate(self, mocked_structures3):
        # GIVEN
        arr1 = numpy.array([[1, 1, 1], [2, 2, 3]])
        arr2 = numpy.array([[2, 3, 2], [3, 4, 5], [6, 7, 8]])
        stc1, stc2, stc3 = mocked_structures3
        structure_new_order = stc1, stc3, stc2
        al1 = MultipleColumnsAlignment(mocked_structures3, arr1)
        al2 = MultipleColumnsAlignment(structure_new_order, arr2)

        # WHEN
        al3 = al1.concatenate(al2)

        # THEN
        assert al3.inds.shape == (4, 3)

    def test_limit_structures(self):
        structures = get_n_mocked_structures(4)
        array = numpy.array([[i] * 4 for i in range(5)])
        alignment = MultipleColumnsAlignment(structures, array)
        stc1, stc2, stc3, stc4 = structures

        # WHEN
        with pytest.raises(ValueError):
            alignment.limit_to_structures(stc1)
        new_pa = alignment.limit_to_structures(stc1, stc2)
        new_ma = alignment.limit_to_structures(stc2, stc1, stc4)

        # THEN
        assert new_pa.inds.shape == (5, 2)
        assert new_ma.inds.shape == (5, 3)

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
        alignment = MultipleColumnsAlignment(structures, arr)

        closed_alignment = alignment.close()

        assert closed_alignment.inds.shape == (4, 4)
        expected_row0 = [1, 1, 1, 1]
        numpy.testing.assert_array_equal(closed_alignment.inds[0], expected_row0)
        expected_row1 = [2, 2, 2, 2]
        numpy.testing.assert_array_equal(closed_alignment.inds[1], expected_row1)
        expected_row2 = [3, 3, DASH, DASH]
        numpy.testing.assert_array_equal(closed_alignment.inds[2], expected_row2)

    def test_close_inconsistent(self):
        arr = numpy.array([[1, 1, DASH], [DASH, 1, 1], [1, DASH, 2],])
        structures = get_n_mocked_structures(3)
        alignment = MultipleColumnsAlignment(structures, arr)

        with pytest.raises(ValueError):
            alignment.close()

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

        alignment = MultipleColumnsAlignment(structures, messy_arr)

        sorted_al = alignment.sort()

        numpy.testing.assert_array_equal(sorted_al.inds, arr)

    def test_get_aligned_inds(self):
        arr = numpy.array([[1, 1, 1, 1], [1, 2, DASH, DASH], [3, 3, 3, 3],])
        structures = get_n_mocked_structures(4)
        stc0, stc1, stc2, stc3 = structures
        alignment = MultipleColumnsAlignment(structures, arr)

        expected_stc0_mer1 = {
            stc1: [1, 2],
            stc2: [1],
            stc3: [1],
        }
        stc0_mer1 = alignment.get_inds_aligned_with(stc0, 1)
        assert stc0_mer1 == expected_stc0_mer1

        expected_stc1_mer2 = {stc0: [1]}
        stc1_mer2 = alignment.get_inds_aligned_with(stc1, 2)
        assert stc1_mer2 == expected_stc1_mer2

        expected_stc0_mer3 = {stc: [3] for stc in structures[1:]}
        stc0_mer3 = alignment.get_inds_aligned_with(stc0, 3)
        assert stc0_mer3 == expected_stc0_mer3

        with pytest.raises(KeyError):
            alignment.get_inds_aligned_with(stc0, 42)


class TestJoinedAlignments:
    def test_join(self):
        pair_alignments1 = get_3_pair_alignments()
        alignment1 = JoinedPairAlignments(pair_alignments1)
        different_pair_alignments = get_3_pair_alignments()
        alignment2 = JoinedPairAlignments(different_pair_alignments)

        single_pa = different_pair_alignments[0]
        already_included = pair_alignments1[0]

        # WHEN
        alignment3 = alignment1.join(alignment2)
        alignment4 = alignment1.join(single_pa)
        alignment5 = alignment1.join(already_included)

        # THEN
        assert len(alignment3.pair_alignments) == 6
        assert len(alignment4.pair_alignments) == 4
        assert len(alignment5.pair_alignments) == 3

    def test_joined_pair_to_column(self, mocked_structures3):
        stc1, stc2, stc3 = mocked_structures3

        pair_array1 = numpy.array([[0, 0], [1, 1], [2, 2]])
        pair_array2 = numpy.array([[0, 0], [1, 1], [3, 3]])
        pair_al1 = PairAlignment([stc1, stc2], pair_array1)
        pair_al2 = PairAlignment([stc1, stc3], pair_array2)

        alignment = JoinedPairAlignments((pair_al1, pair_al2))

        column_alignment = alignment.to_columns()

        assert column_alignment.inds.shape == (6, 3)
        for row in column_alignment.inds:
            assert DASH in row

    def test_limit_to_2_structures(self):
        pair_alignments = get_3_pair_alignments()
        alignment = JoinedPairAlignments(pair_alignments)
        stc1, stc2 = pair_alignments[0].structures

        limited_pa = alignment.limit_to_structures(stc1, stc2)

        assert isinstance(limited_pa, PairAlignment)

    def test_limit_to_multiple_structures(self):
        structures = get_n_mocked_structures(4)
        pas = []
        for i, stc1 in enumerate(structures):
            for stc2 in structures[i + 1 :]:
                arr = numpy.array([[l, l] for l in range(6)])
                pa = PairAlignment((stc1, stc2), arr)
                pas.append(pa)
        alignment = JoinedPairAlignments(pas)

        limited = alignment.limit_to_structures(*structures[:-1])

        assert len(limited.pair_alignments) == 3


class TestPairAlignment:
    @pytest.fixture
    def pair_alignment(self):
        stc1, stc2 = get_n_mocked_structures(2)
        pa = get_trivial_pair_alignment(stc1, stc2)
        return pa

    def test_prune(self):
        payload = numpy.array([[0, DASH], [1, 0]])
        alignment = PairAlignment([None, None], payload)
        alignment = alignment.prune()

        assert alignment.inds.shape == (1, 2)

    def test_init(self, mocked_structures3):
        arr = numpy.array([[1, 2], [3, 4], [5, 6],])
        with pytest.raises(ValueError):
            PairAlignment(mocked_structures3, arr)
        pa = PairAlignment(mocked_structures3[:-1], arr)

        assert mocked_structures3[0] in pa.mer_map
        assert mocked_structures3[1] in pa.mer_map
        assert pa.inds.shape == (3, 2)

    def test_conversion(self, pair_alignment):
        same_pa = pair_alignment.to_columns()
        also_this_pa = pair_alignment.to_joined_pairs()

        assert pair_alignment is same_pa
        assert pair_alignment is also_this_pa

    def test_limit(self, pair_alignment):
        (stc3,) = get_n_mocked_structures(1, start=2)
        with pytest.raises(ValueError):
            pair_alignment.limit_to_structures(stc3)

        new_pa = pair_alignment.limit_to_structures(*pair_alignment.structures)

        assert new_pa is pair_alignment

    def test_equality(self, pair_alignment):
        very_similar_pa = get_trivial_pair_alignment(*pair_alignment.structures)

        assert pair_alignment == very_similar_pa

        stc3, stc4 = get_n_mocked_structures(2, start=2)
        pa_w_different_structures = get_trivial_pair_alignment(stc3, stc4)
        assert not pa_w_different_structures == pair_alignment
        assert pa_w_different_structures != pair_alignment

        pa_different_mers = get_trivial_pair_alignment(*pair_alignment.structures)
        pa_different_mers.inds = pa_different_mers.inds[:3]

        assert not pa_different_mers == pair_alignment

    def test_transit_trivial(self):
        stc1, stc2, stc3 = get_n_mocked_structures(3)
        pa1_2 = get_trivial_pair_alignment(stc1, stc2)
        pa2_3 = get_trivial_pair_alignment(stc2, stc3)
        expected1_3 = get_trivial_pair_alignment(stc1, stc3)

        # WHEN
        pa1_3 = pa1_2.transit(pa2_3)

        # THEN
        numpy.testing.assert_array_equal(pa1_3.inds, expected1_3.inds)

    def test_transit_with_dash(self):
        stc1, stc2, stc3 = get_n_mocked_structures(3)
        arr1 = numpy.array([[1, 1], [DASH, 2]])
        arr2 = numpy.array([[1, 0], [2, 2]])
        pa1_2 = PairAlignment((stc1, stc2), arr1)
        pa2_3 = PairAlignment((stc2, stc3), arr2)
        pa1_3 = pa1_2.transit(pa2_3)

        assert pa1_3.inds.shape == (1, 2)
        assert DASH not in pa1_3.inds.ravel()

    def test_consistency_positive(self):
        stc1, stc2 = get_n_mocked_structures(2)
        pa1 = get_trivial_pair_alignment(stc1, stc2)
        pa2 = get_trivial_pair_alignment(stc1, stc2)
        pa3 = get_trivial_pair_alignment(stc2, stc1)
        arr = numpy.array([[1, 1], [42, 42],])
        pa4 = PairAlignment([stc1, stc2], arr)

        assert pa1.is_consistent_with(pa2)
        assert pa2.is_consistent_with(pa3)
        assert pa1.is_consistent_with(pa3)
        assert pa2.is_consistent_with(pa1)
        assert pa3.is_consistent_with(pa1)
        assert pa1.is_consistent_with(pa4)

    def test_consistency_double_mer_alignment(self):
        stc1, stc2 = get_n_mocked_structures(2)
        pa1 = get_trivial_pair_alignment(stc1, stc2)
        arr = numpy.array([[1, 1], [1, 42],])
        pa2 = PairAlignment([stc1, stc2], arr)

        with pytest.raises(ValueError):
            pa1.is_consistent_with(pa2)

    def test_consistency_negative(self):
        stc1, stc2 = get_n_mocked_structures(2)
        pa1 = get_trivial_pair_alignment(stc1, stc2)
        arr = numpy.array([[1, 1], [2, 42],])
        pa2 = PairAlignment([stc1, stc2], arr)
        assert not pa1.is_consistent_with(pa2)

    def test_consistency_different_structures(self):
        stc1, stc2, stc3, stc4 = get_n_mocked_structures(4)
        pa1 = get_trivial_pair_alignment(stc1, stc2)
        pa2 = get_trivial_pair_alignment(stc2, stc3)
        pa3 = get_trivial_pair_alignment(stc3, stc4)

        assert pa1.is_consistent_with(pa2)
        assert pa1.is_consistent_with(pa3)
