import os.path

import pytest
from Bio.PDB import PDBParser

from pydesc import config
from pydesc import numberconverter
from pydesc.numberconverter import PDBid

TWISTED_EXAMPLE = (
    (("b", 1, None), ("b", 1, "a"), ("b", 2, None), ("a", 8, None)),
    (("b", 1, "a"), ("b", 2, None), ("b", 3, None), ("c", 6, None), ("a", 7, None),),
)
SIMPLE_EXAMPLE = (
    (("b", 1, None), ("b", 2, None), ("a", 8, None)),
    (("b", 1, None), ("b", 2, None), ("a", 8, None)),
)


@pytest.mark.parametrize("ids,ids2", (TWISTED_EXAMPLE, SIMPLE_EXAMPLE))
class TestSmithWaterman:
    def test_matrix_creation(self, ids, ids2):
        res_mtx = numberconverter.build_smith_waterman_matrix(ids, ids2)

        values_stored = res_mtx[1:, 1:, 3:].reshape((-1, 3))
        values_stored_b = set()
        for chain, ind, i_code in values_stored:
            chain = chr(chain)
            i_code = chr(i_code) if i_code else None
            values_stored_b.add((chain, ind, i_code))

        expected_values = set(ids) | set(ids2)
        assert values_stored_b == expected_values
        expected_shape = len(ids2) + 1, len(ids) + 1, 6
        assert res_mtx.shape == expected_shape

    def test_back_trace(self, ids, ids2):
        res_mtx = numberconverter.build_smith_waterman_matrix(ids, ids2)
        result_ids = numberconverter.go_backwards(res_mtx)

        expected_ids = set(ids) | set(ids2)
        assert set(result_ids) == expected_ids


class TestConverter:
    @pytest.mark.system
    def test_simple_single_structure(self, structures_dir):
        structure_name = "5MPV.pdb"
        file_path = os.path.join(structures_dir, "prots_only", structure_name)
        pdb_structure = PDBParser(QUIET=True).get_structure(structure_name, file_path)

        nc_factory = numberconverter.NumberConverterFactory()
        nc = nc_factory.from_pdb_models([mod for mod in pdb_structure])

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
                    assert nc.pdb2ind[pdb_id] == int(id_[3][1]) - 11
                    # in case of this protein (5mpv) that works

    @pytest.mark.system
    def test_solvent(self, structures_dir):
        old_solvent_setting = list(config.ConfigManager.chemistry.solvent)
        config.ConfigManager.chemistry.solvent = []

        structure_name = "5MPV.pdb"
        file_path = os.path.join(structures_dir, "prots_only", structure_name)
        pdb_structure = PDBParser(QUIET=True).get_structure(structure_name, file_path)

        nc_factory = numberconverter.NumberConverterFactory()
        nc = nc_factory.from_pdb_models([mod for mod in pdb_structure])

        assert len(nc.pdb2ind) == 149
        assert len(nc.ind2pdb) == 149

        for mod in pdb_structure:
            for ch in mod:
                for res in ch:
                    id_ = res.get_full_id()
                    ic = id_[3][2]
                    pdb_id = (id_[2], id_[3][1], ic if ic != " " else None)
                    assert pdb_id in nc.pdb2ind
        config.ConfigManager.chemistry.solvent = old_solvent_setting
        # TODO: this is exactly why implicit config is a bad idea.
        #   can we get rid of this?

    @pytest.mark.system
    def test_nmr_20_structures(self, structures_dir):
        structure_name = "2JRM.pdb"
        file_path = os.path.join(structures_dir, "prots_only_nmr", structure_name)
        pdb_structure = PDBParser(QUIET=True).get_structure(structure_name, file_path)

        nc_factory = numberconverter.NumberConverterFactory()
        nc = nc_factory.from_pdb_models([mod for mod in pdb_structure])

        assert len(nc.pdb2ind) == 60
        assert len(nc.ind2pdb) == 60


class TestPDBid:
    def test_from_string_simple(self):
        pdb_id = PDBid.create_from_string("A:12")
        assert pdb_id.icode is None
        assert pdb_id.chain == "A"
        assert pdb_id.ind == 12

    def test_from_string_icode(self):
        pdb_id = PDBid.create_from_string("A:12i")
        assert pdb_id.icode is "i"
        assert pdb_id.chain == "A"
        assert pdb_id.ind == 12

    def test_from_string_long_chain(self):
        pdb_id = PDBid.create_from_string("AB:12i")
        assert pdb_id.icode is "i"
        assert pdb_id.chain == "AB"
        assert pdb_id.ind == 12
