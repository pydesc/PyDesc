import os.path
from unittest.mock import MagicMock

import numpy
import pytest

from pydesc.alignment.base import JoinedPairAlignments
from pydesc.alignment.base import MultipleColumnsAlignment
from pydesc.alignment.base import PairAlignment
from pydesc.alignment.loaders import CSVLoader
from pydesc.alignment.loaders import DASH
from pydesc.alignment.loaders import FASTALoader
from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file


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


class TestLoaders:
    def test_csv_loader_multi(self, alignments_dir):
        path = os.path.join(alignments_dir, "csv", "artificial_multi.csv")
        loader = CSVLoader(path)
        structures = [MagicMock() for _ in range(4)]
        for mocked_structure in structures:
            mocked_structure.converter.get_ind.side_effect = [i for i in range(1, 7)]
        alignment = loader.load_alignment(structures)

        assert alignment.inds.shape == (6, 4)
        expected_not_nans = [5, 4, 3, 5]
        real_not_nans = numpy.count_nonzero(alignment.inds != DASH, axis=0)
        numpy.testing.assert_equal(real_not_nans, expected_not_nans)

        assert len(alignment.structures) == 4

        stc_w_icodes = structures[2]
        calls = stc_w_icodes.converter.get_ind.call_args_list
        for call in calls[1:]:
            arg = call.args[0]
            assert arg.icode

        stc_wo_icodes = structures[1]
        calls = stc_wo_icodes.converter.get_ind.call_args_list
        for call in calls[1:]:
            arg = call.args[0]
            assert arg.icode is None

        stc_w_2letter_chain = structures[3]
        calls = stc_w_2letter_chain.converter.get_ind.call_args_list
        for call in calls[1:]:
            arg = call.args[0]
            assert len(arg.chain) == 2

        assert isinstance(alignment, MultipleColumnsAlignment)

    def test_csv_loader_pair(self, alignments_dir):
        path = os.path.join(alignments_dir, "csv", "artificial_pair.csv")
        loader = CSVLoader(path)
        structures = [MagicMock() for _ in range(2)]
        alignment = loader.load_alignment(structures)

        assert isinstance(alignment, PairAlignment)

    def test_csv_wrong_structures(self, alignments_dir):
        path = os.path.join(alignments_dir, "csv", "artificial_pair.csv")
        loader = CSVLoader(path)
        structures = [MagicMock() for _ in range(4)]

        with pytest.raises(ValueError):
            loader.load_alignment(structures)

    def test_pal_artificial_multi(self, alignments_dir):
        path = os.path.join(alignments_dir, "pal", "artificial_multi.pal")
        structures = [MagicMock() for _ in range(3)]
        mocked_chain = MagicMock(chain_name="K")
        structures[0].chains = [mocked_chain]

        loader = PALLoader(path)
        metadata = loader.read_metadata()
        for mol in "MOL1", "MOL2", "MOL3":
            assert mol in metadata["labels"]
        alignment = loader.load_alignment(structures)

        assert len(alignment.pair_alignments) == 3

        def get_stc_calls(index):
            return structures[index].converter.get_ind.call_args_list

        stc1_call_args_set = {call.args[0] for call in get_stc_calls(0)}
        for i in stc1_call_args_set:
            assert i.chain == "K"
        assert ("K", 31, "i") in stc1_call_args_set

        stc2_call_args_set = {call.args[0] for call in get_stc_calls(1)}
        for i in stc2_call_args_set:
            assert i.chain == "AB"

        stc3_call_args_set = {call.args[0] for call in get_stc_calls(2)}
        stc3_chains = {i.chain for i in stc3_call_args_set}
        assert stc3_chains == {"A", "B"}

    @pytest.mark.system
    def test_pal_pair(self, alignments_dir):
        stc_path = os.path.join(alignments_dir, "structures", "sars_pair.pdb")
        structures = get_structures_from_file(stc_path)
        path = os.path.join(alignments_dir, "pal", "sars_pair.pal")

        loader = PALLoader(path)
        alignment = loader.load_alignment(structures)

        assert isinstance(alignment, PairAlignment)
        assert alignment.inds.shape == (148, 2)

        expected_first = numpy.array([0, 2])
        numpy.testing.assert_equal(alignment.inds[0], expected_first)

        expected_7_8 = numpy.array(
            [[8, 10], [9, 12]]
        )  # there is 1 residue gap in right structure
        numpy.testing.assert_equal(alignment.inds[(8, 9), :], expected_7_8)

    def test_fasta_artificial_multi(self, alignments_dir):
        path = os.path.join(alignments_dir, "fasta", "artificial_multi.fasta")
        structures = [MagicMock(name=f"mol{i}") for i in range(5)]
        for stc in structures[:4]:
            stc.__getitem__.return_value = [MagicMock(ind=i) for i in range(10, 15)]
        structures[-1].__getitem__.return_value = [MagicMock(ind=1), MagicMock(ind=2)]

        loader = FASTALoader(path)
        alignment = loader.load_alignment(structures)

        for label in ("mol1", "mol2", "mol3_chainAB", "mol4", "mol5"):
            assert label in loader.structure_labels

        assert alignment.inds.shape == (4, 5)

        expected = numpy.array(
            [
                [DASH, DASH, DASH, 10.0, 1.0],
                [10.0, 11.0, DASH, 11.0, 2.0],
                [11.0, 12.0, 10.0, 12.0, 1.0],
                [12.0, 13.0, 11.0, 13.0, DASH],
            ]
        )
        numpy.testing.assert_equal(alignment.inds, expected)

    @pytest.mark.system
    def test_fasta_pal_from_dama(self, alignments_dir):
        file_name = "mers_sars2.%s"
        fasta_path = os.path.join(alignments_dir, "fasta", file_name % "fasta")
        pal_path = os.path.join(alignments_dir, "pal", file_name % "pal")
        mers_path = os.path.join(alignments_dir, "structures", "mers.pdb")
        (mers_stc,) = get_structures_from_file(mers_path)
        sars_path = os.path.join(alignments_dir, "structures", "sars2.pdb")
        (sars_stc,) = get_structures_from_file(sars_path)
        structures_map = {
            "MOL1": mers_stc,
            "MOL2": sars_stc,
        }

        fasta_loader = FASTALoader(fasta_path)
        fasta = fasta_loader.load_partial_alignment(structures_map)
        pal_loader = PALLoader(pal_path)
        pal = pal_loader.load_partial_alignment(structures_map)

        numpy.testing.assert_equal(fasta.inds, pal.inds)
        assert fasta.structures == pal.structures


class TestColumnAlignment:
    def test_simplify_pair(self):
        payload = numpy.array([[0, DASH], [1, 0]])
        alignment = PairAlignment([None, None], payload)
        alignment.prune()

        assert alignment.inds.shape == (1, 2)

    def test_simplify_multi(self):
        payload = numpy.array([[0, DASH, 0], [1, 0, 1], [DASH, DASH, 4], ])
        alignment = MultipleColumnsAlignment([None, None, None], payload)
        alignment.prune()

        assert alignment.inds.shape == (2, 3)

    @pytest.fixture
    def trivial_array3(self):
        array = numpy.array([
            [1, 1, 1],
            [2, 2, 2],
            [3, 3, 3],
            [4, 4, 4]
        ])
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
        arr1 = numpy.array([
            [1, 1, 1],
            [2, 2, 3]
        ])
        arr2 = numpy.array([
            [2, 3, 2],
            [3, 4, 5],
            [6, 7, 8]
        ])
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
            for stc2 in structures[i + 1:]:
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

    def test_init(self, mocked_structures3):
        with pytest.raises(ValueError):
            PairAlignment(mocked_structures3, None)
        pa = PairAlignment(mocked_structures3[:-1], None)

        assert pa.inds is None

    def test_conversion(self, pair_alignment):
        same_pa = pair_alignment.to_columns()
        also_this_pa = pair_alignment.to_joined_pairs()

        assert pair_alignment is same_pa
        assert pair_alignment is also_this_pa

    def test_limit(self, pair_alignment):
        stc3, = get_n_mocked_structures(1, start=2)
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
