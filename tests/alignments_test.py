import os.path
from unittest.mock import MagicMock

import numpy
import pytest

from pydesc.alignment.base import MultipleColumnsAlignment
from pydesc.alignment.base import PairAlignment
from pydesc.alignment.loaders import CSVLoader
from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file


class TestLoaders:
    def test_csv_loader_multi(self, alignments_dir):
        path = os.path.join(alignments_dir, "csv", "artificial_multi.csv")
        loader = CSVLoader(path)
        structures = [MagicMock() for _ in range(4)]
        for mocked_structure in structures:
            mocked_structure.converter.get_ind.side_effect = [i for i in range(1, 7)]
        alignment = loader.create_alignment(structures)

        assert alignment.inds.shape == (6, 4)
        expected_not_nans = [5, 4, 3, 5]
        real_not_nans = numpy.count_nonzero(~numpy.isnan(alignment.inds), axis=0)
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
        alignment = loader.create_alignment(structures)

        assert isinstance(alignment, PairAlignment)

    def test_csv_wrong_structures(self, alignments_dir):
        path = os.path.join(alignments_dir, "csv", "artificial_pair.csv")
        loader = CSVLoader(path)
        structures = [MagicMock() for _ in range(4)]

        with pytest.raises(Exception):
            loader.create_alignment(structures)

    def test_pal_artificial_multi(self, alignments_dir):
        path = os.path.join(alignments_dir, "pal", "artificial_multi.pal")
        structures = [MagicMock() for _ in range(3)]
        mocked_chain = MagicMock(chain_name="K")
        structures[0].chains = [mocked_chain]

        loader = PALLoader(path)
        metadata = loader.read_metadata()
        for mol in "MOL1", "MOL2", "MOL3":
            assert mol in metadata['labels']
        alignment = loader.create_alignment(structures)

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
        alignment = loader.create_alignment(structures)

        assert isinstance(alignment, PairAlignment)
        assert alignment.inds.shape == (148, 2)

        expected_first = numpy.array([0, 2])
        numpy.testing.assert_equal(alignment.inds[0], expected_first)

        expected_7_8 = numpy.array(
            [[8, 10], [9, 12]]
        )  # there is 1 residue gap in right structure
        numpy.testing.assert_equal(alignment.inds[(8, 9), :], expected_7_8)
