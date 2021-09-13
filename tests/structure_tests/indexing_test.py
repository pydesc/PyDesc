from pathlib import Path

import pytest

from pydesc.api.structure import get_structures_from_file
from pydesc.structure.topology import PartialStructure
from pydesc.structure.topology import Segment


@pytest.fixture(scope="session")
def rna(structures_dir):
    pth = Path(structures_dir) / "rna_only" / "1KIS.pdb"
    structure = get_structures_from_file(str(pth))
    return structure[0]


@pytest.fixture(scope="session")
def protein(structures_dir):
    pth = Path(structures_dir) / "prots_only" / "4ONK.pdb"
    structure = get_structures_from_file(str(pth))
    return structure[0]


class TestInds:
    def test_mers_positive(self, rna):
        mer = rna[12]
        assert mer.ind == 12

    def test_mers_negative(self, rna):
        with pytest.raises(IndexError):
            rna[42]

    def test_slice_full(self, rna):
        segment = rna[3:14]
        assert isinstance(segment, Segment)
        assert len(segment) == 12

        partial = rna[14:18]
        assert isinstance(partial, PartialStructure)
        assert len(partial) == 5

    def test_slice_step(self, rna):
        with pytest.raises(ValueError):
            rna[1:2:2]

    def test_slice_start(self, rna):
        sub = rna[29:]
        inds = [i.ind for i in sub]
        expected = [29, 30, 31]
        assert inds == expected

    def test_slice_stop(self, rna):
        sub = rna[:4]
        inds = [i.ind for i in sub]
        expected = [0, 1, 2, 3, 4]
        assert inds == expected

    def test_slice_all(self, rna):
        sub = rna[:]
        assert tuple(rna) == tuple(sub)

    def test_iterable(self, rna):
        expected = [1, 2, 3]
        sub = rna[expected]
        inds = [i.ind for i in sub]
        assert inds == expected

    def test_iterable_w_repeats(self, rna):
        expected = [1, 2, 3]
        sub = rna[expected * 2]
        inds = [i.ind for i in sub]
        assert inds == expected


class TestPDBids:
    def test_mer(self, rna):
        expected = "B:22"
        mer = rna.pdb_ids[expected]
        pdb_id = rna.converter.get_pdb_id(mer.ind)
        assert pdb_id.format(chain=True) == expected

    def test_chain_wildcard(self, rna):
        expected = "B:"
        chain = rna.pdb_ids[expected]
        assert chain.chain_name == "B"

    def test_mer_wildcard(self, protein):
        sub = protein.pdb_ids["3"]
        assert len(sub) == 2
        chains = set()
        for mer in sub:
            pdb_id = protein.converter.get_pdb_id(mer.ind)
            assert pdb_id.format(chain=False) == "3"
            chains.add(pdb_id.chain)
        assert chains == {"A", "B"}

    def test_slice(self, rna):
        sub = rna.pdb_ids["B:22":"B:30"]
        assert len(sub) == 9
        assert isinstance(sub, Segment)
        mers = tuple(sub)
        pdb_id = rna.converter.get_pdb_id(mers[0].ind)
        assert pdb_id.format(chain=True) == "B:22"
        pdb_id = rna.converter.get_pdb_id(mers[1].ind)
        assert pdb_id.format(chain=True) == "B:23"
        pdb_id = rna.converter.get_pdb_id(mers[-1].ind)
        assert pdb_id.format(chain=True) == "B:30"

    def test_slice_partial(self, rna):
        sub1 = rna.pdb_ids["B:22":]
        sub2 = rna.pdb_ids[:"B:22"]
        len1 = len(sub1)
        len2 = len(sub2)
        assert len1 + len2 == len(rna) + 1

    def test_slice_wildcard_negative(self, rna):
        with pytest.raises(ValueError):
            rna.pdb_ids["22":"B:24"]

        with pytest.raises(ValueError):
            rna.pdb_ids["B:22":"24"]

    def test_iterable(self, rna):
        sub = rna.pdb_ids[["B:22", "B:23"]]
        assert len(sub) == 2

    def test_iterable_wildcard_negative(self, rna):
        with pytest.raises(ValueError):
            sub = rna.pdb_ids[["B:22", "23"]]
