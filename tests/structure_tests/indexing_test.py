from pydesc.api.structure import get_structures_from_file

from pydesc.structure.topology import Segment, PartialStructure
from pathlib import Path
import pytest


@pytest.fixture(scope="session")
def structure(structures_dir):
    pth = Path(structures_dir) / "rna_only" / "1KIS.pdb"
    structure = get_structures_from_file(str(pth))
    return structure[0]


class TestInds:
    def test_mers_positive(self, structure):
        mer = structure[12]
        assert mer.ind == 12

    def test_mers_negative(self, structure):
        with pytest.raises(IndexError):
            structure[42]

    def test_slice_full(self, structure):
        segment = structure[3:14]
        assert isinstance(segment, Segment)
        assert len(segment) == 12

        partial = structure[14:18]
        assert isinstance(partial, PartialStructure)
        assert len(partial) == 5

    def test_slice_step(self, structure):
        with pytest.raises(ValueError):
            structure[1:2:2]

    def test_slice_start(self, structure):
        sub = structure[29:]
        inds = [i.ind for i in sub]
        expected = [29, 30, 31]
        assert inds == expected

    def test_slice_stop(self, structure):
        sub = structure[:4]
        inds = [i.ind for i in sub]
        expected = [0, 1, 2, 3, 4]
        assert inds == expected

    def test_slice_all(self, structure):
        sub = structure[:]
        assert tuple(structure) == tuple(sub)

    def test_iterable(self, structure):
        expected = [1, 2, 3]
        sub = structure[expected]
        inds = [i.ind for i in sub]
        assert inds == expected

    def test_iterable_w_repeats(self, structure):
        expected = [1, 2, 3]
        sub = structure[expected * 2]
        inds = [i.ind for i in sub]
        assert inds == expected


class TestPDBids:
    def test_mer(self, structure):
        expected = "B:22"
        mer = structure.pdb_ids[expected]
        pdb_id = structure.converter.get_pdb_id(mer.ind)
        assert pdb_id.format(chain=True) == expected

    @pytest.mark.skip(reason="TBD")
    def test_mer_wildcard(self, structure):
        sub = structure.pdb_ids["22"]
        assert len(sub) == 2
        chains = set()
        for mer in sub:
            pdb_id = structure.converter.get_pdb_id(mer.ind)
            assert pdb_id.format(chain=False) == "22"
            chains.add(pdb_id.chain)
        assert chains == {"A", "B"}

    @pytest.mark.skip(reason="TBD")
    def test_slice(self, structure):
        sub = structure.pdb_ids["B:22":"B:30"]
        assert len(sub) == 9
        assert isinstance(sub, Segment)
        mers = tuple(sub)
        pdb_id = structure.converter.get_pdb_id(mers[0].ind)
        assert pdb_id.format(chain=True) == "B:22"
        pdb_id = structure.converter.get_pdb_id(mers[1].ind)
        assert pdb_id.format(chain=True) == "B:23"
        pdb_id = structure.converter.get_pdb_id(mers[-1].ind)
        assert pdb_id.format(chain=True) == "B:30"

    @pytest.mark.skip(reason="TBD")
    def test_slice_partial(self, structure):
        sub1 = structure.pdb_ids["B:22":]
        sub2 = structure.pdb_ids[:"B:22"]
        len1 = len(sub1)
        len2 = len(sub2)
        assert len1 + len2 == len(structure) + 1

    @pytest.mark.skip(reason="TBD")
    def test_slice_wildcard_positive(self, structure):
        sub = structure.pdb_ids["22":"24"]
        assert len(sub) == 3

    @pytest.mark.skip(reason="TBD")
    def test_slice_wildcard_negative(self, structure):
        with pytest.raises(ValueError):
            structure.pdb_ids["22":"B:24"]

        with pytest.raises(ValueError):
            structure.pdb_ids["B:22":"24"]

    def test_iterable(self, structure):
        sub = structure.pdb_ids[["B:22", "B:23"]]
        assert len(sub) == 2

    @pytest.mark.skip(reason="TBD")
    def test_iterable_wildcard(self, structure):
        sub = structure.pdb_ids[["B:22", "23"]]
        assert len(sub) == 3
