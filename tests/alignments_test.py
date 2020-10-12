import os.path
from unittest.mock import MagicMock

import numpy
import pytest

from pydesc.alignment.base import MultipleColumnsAlignment
from pydesc.alignment.base import PairAlignment
from pydesc.alignment.loaders import CSVLoader
from pydesc.alignment.loaders import FASTALoader
from pydesc.alignment.loaders import PALLoader, DASH
from pydesc.api.structure import get_structures_from_file


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
        alignment.simplify()

        assert alignment.inds.shape == (1, 2)

    def test_simplify_multi(self):
        payload = numpy.array([[0, DASH, 0], [1, 0, 1], [DASH, DASH, 4],])
        alignment = MultipleColumnsAlignment([None, None, None], payload)
        alignment.simplify()

        assert alignment.inds.shape == (2, 3)
