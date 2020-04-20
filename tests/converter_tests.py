import os.path

import pytest
from Bio.PDB import PDBParser

from pydesc import config
from pydesc import numberconverter
from tests.conftest import TEST_STRUCTURES_DIR


class TestSmithWaterman:
    def test_matrix_creation(self):
        ids = (("b", 1, None), ("b", 1, "a"), ("b", 2, None), ("a", 8, None))
        ids2 = (
            ("b", 1, "a"),
            ("b", 2, None),
            ("b", 3, None),
            ("c", 6, None),
            ("a", 7, None),
        )
        res_mtx = numberconverter.build_smith_waterman_matrix(ids, ids2)

        assert res_mtx.shape == (6, 5, 6)
        # assert chains are stored
        for i in res_mtx[1:, 1:, 3].ravel():
            assert chr(i) in ("a", "b", "c")
        # assert i_codes are 0
        for i in res_mtx[3:, 3:, 5].ravel():
            assert i == 0
        # asert other icodes are not 0
        assert res_mtx[2, 2, 5] != 0

    def test_back_trace(self):
        ids = (("b", 1, None), ("b", 1, "a"), ("b", 2, None), ("a", 8, None))
        ids2 = (
            ("b", 1, "a"),
            ("b", 2, None),
            ("b", 3, None),
            ("c", 6, None),
            ("a", 7, None),
        )
        res_mtx = numberconverter.build_smith_waterman_matrix(ids, ids2)
        ids = numberconverter.go_backwards(res_mtx)

        assert ids == (ids[0],) + ids2[1:] + (ids[-1],)


class TestConverter:
    @pytest.mark.system
    def test_simple_singe_structure(self):
        pdb_structure = PDBParser(QUIET=True).get_structure(
            "5MPV.pdb", os.path.join(TEST_STRUCTURES_DIR, "prots_only", "5MPV.pdb")
        )
        nc = numberconverter.NumberConverter([mod for mod in pdb_structure])

        assert len(nc.pdb2ind) == 78
        assert len(nc.ind2pdb) == 78

        for mod in pdb_structure:
            for ch in mod:
                for res in ch:
                    id_ = res.get_full_id()
                    if "W" in id_[3][0]:
                        continue
                        # skip solvent (according to Bio.PDB)
                    ic = id_[3][2]
                    pdb_id = (id_[2], id_[3][1], ic if ic != " " else None)
                    assert nc.pdb2ind[pdb_id] == int(id_[3][1]) - 10
                    # in case of this protein (5mpv) that works

    @pytest.mark.system
    def test_solvent(self):
        config.ConfigManager.mers.solvent = []
        pdb_structure = PDBParser(QUIET=True).get_structure(
            "5MPV.pdb", os.path.join(TEST_STRUCTURES_DIR, "prots_only", "5MPV.pdb")
        )
        nc = numberconverter.NumberConverter([mod for mod in pdb_structure])

        assert len(nc.pdb2ind) == 149
        assert len(nc.ind2pdb) == 149

        for mod in pdb_structure:
            for ch in mod:
                for res in ch:
                    id_ = res.get_full_id()
                    ic = id_[3][2]
                    pdb_id = (id_[2], id_[3][1], ic if ic != " " else None)
                    assert pdb_id in nc.pdb2ind
        config.ConfigManager.mers.solvent = ["HOH"]

    @pytest.mark.system
    def test_nmr_20_structures(self):
        pdb_structure = PDBParser(QUIET=True).get_structure(
            "1A24.pdb", os.path.join(TEST_STRUCTURES_DIR, "prots_only", "1A24.pdb")
        )
        nc = numberconverter.NumberConverter([mod for mod in pdb_structure])

        assert len(nc.pdb2ind) == 189
        assert len(nc.ind2pdb) == 189
