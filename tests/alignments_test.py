import pytest
import os.path
from pydesc.alignment.base import PairAlignment, MultipleColumnsAlignment
from unittest.mock import MagicMock

import numpy

from pydesc.alignment.loaders import CSVLoader


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
