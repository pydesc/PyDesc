from unittest.mock import MagicMock

import numpy
import pytest
from common import get_n_mocked_structures
from common import get_trivial_array

from pydesc.alignment.base import DASH
from pydesc.alignment.base import Alignment
from pydesc.alignment.base import IncorrectAlignmentError


def test_dash_repr():
    assert repr(DASH) == "<->"


class TestAlignment:
    def test_init(self, triple_alignment):
        assert triple_alignment

    def test_init_one_structure(self):
        (stc,) = get_n_mocked_structures(1)
        arr = numpy.zeros(10).reshape(10, 1)
        with pytest.raises(IncorrectAlignmentError):
            Alignment([stc], arr)

    def test_init_wrong_arr(self):
        stcs = get_n_mocked_structures(2)
        arr = get_trivial_array(30, 3)
        with pytest.raises(IncorrectAlignmentError):
            Alignment(stcs, arr)

    def test_len(self, triple_alignment):
        assert len(triple_alignment) == 10

    def test_get_structures(self, triple_alignment):
        recovered_structures = triple_alignment.get_structures()
        assert len(recovered_structures) == 3
        for structure in recovered_structures:
            assert structure.name.startswith("stc")

    def test_get_table(self, triple_alignment):
        table = triple_alignment.get_inds_table()
        assert table.shape == (10, 3)
        table[0, 0] = 42
        new_table = triple_alignment.get_inds_table()
        assert 42 not in new_table

    def test_slicing(self, triple_alignment):
        stc1, stc2, stc3 = triple_alignment.get_structures()
        subalignment = triple_alignment[(stc1, 3):(stc1, 4)]
        assert len(subalignment) == 2

    def test_relative_slicing(self, triple_alignment):
        stc1, stc2, stc3 = triple_alignment.get_structures()
        subalignment = triple_alignment[None:(stc1, 4)]
        assert len(subalignment) == 5

    def test_get_single_item(self, triple_alignment):
        stc1, stc2, stc3 = triple_alignment.get_structures()
        subalignment = triple_alignment[(stc1, 4)]
        assert len(subalignment) == 1

    def test_get_list_of_inds(self, triple_alignment):
        stc1, stc2, stc3 = triple_alignment.get_structures()
        subalignment = triple_alignment[(stc1, [1, 2])]
        assert len(subalignment) == 2

    def test_pdb_slicer_slice(self, triple_alignment):
        stc0, stc1, stc2 = triple_alignment.get_structures()
        stc1.pdb_ids.__getitem__.side_effect = [MagicMock(ind=1), MagicMock(ind=3)]
        subalignment = triple_alignment.pdb_ids["stc1:A:1":"stc1:A:3"]
        assert len(subalignment) == 3

    def test_pdb_slicer_relative_slice(self, triple_alignment):
        stc0, stc1, stc2 = triple_alignment.get_structures()
        stc1.pdb_ids.__getitem__.side_effect = [MagicMock(ind=3)]
        subalignment = triple_alignment.pdb_ids[None:"stc1:A:3"]
        assert len(subalignment) == 4

    def test_pdb_slicer_single(self, triple_alignment):
        stc0, stc1, stc2 = triple_alignment.get_structures()
        stc1.pdb_ids.__getitem__.side_effect = [MagicMock(ind=3)]
        subalignment = triple_alignment.pdb_ids["stc1:A:3"]
        assert len(subalignment) == 1
        assert set(*subalignment._table) == {3}

    def test_pdb_wrong_structure(self, triple_alignment):
        with pytest.raises(KeyError):
            triple_alignment.pdb_ids["wrong-structure:A:3"]

    def test_repr(self, triple_alignment):
        repr_str = repr(triple_alignment)
        assert repr_str.startswith("<")
        assert repr_str.endswith(">")

    def test_iter(self, triple_alignment):
        for ind_s1, ind_s2, ind_s3 in triple_alignment:
            assert isinstance(ind_s1, int)
            assert isinstance(ind_s2, int)
            assert isinstance(ind_s3, int)
