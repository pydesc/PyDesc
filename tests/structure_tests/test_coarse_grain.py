from os.path import join as path_join

import pytest

from pydesc.chemistry.base import Atom
from pydesc.chemistry.factories import BioPythonAtomSetFactory
from pydesc.chemistry.martini import MartiniResidue
from pydesc.structure import StructureLoader


@pytest.fixture
def martini_structure_path(structures_dir):
    return path_join(structures_dir, "martini")


def test_martini(martini_structure_path):
    factory = BioPythonAtomSetFactory(classes=[MartiniResidue])
    loader = StructureLoader(atom_set_factory=factory)
    path = path_join(martini_structure_path, "gpcr_d.pdb")
    stc = loader.load_structures(path=path)[0]

    assert len(stc) > 0
    for residue in stc:
        assert isinstance(residue, MartiniResidue)
        assert len(tuple(residue.iter_bb_atoms())) == 1
        assert len(tuple(residue.iter_nbb_atoms())) == (len(residue) - 1)
        assert isinstance(residue.last_sc, Atom)
        if residue.name in ("GLY", "ALA"):
            assert residue.last_sc == residue.atoms["BB"]
            continue
        last_cs = max(
            [k for k in residue.atoms if "SC" in k], key=lambda name: int(name[-1])
        )
        assert residue.last_sc == residue.atoms[last_cs]

    assert len(stc.chains) == 2
    ch_A = stc.chains[0]
    prev_res = ch_A[0]
    assert prev_res.prev_mer is None
    for residue in tuple(ch_A)[1:]:
        assert prev_res.next_mer == residue
        assert residue.prev_mer == prev_res
        prev_res = residue
    assert prev_res.next_mer is None