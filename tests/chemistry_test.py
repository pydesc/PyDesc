import os.path

import Bio.PDB
import numpy
import pytest

from pydesc.chemistry.base import Atom
from pydesc.chemistry.base import AtomSet
from pydesc.chemistry.base import Ligand
from pydesc.chemistry.base import Mer
from pydesc.chemistry.base import Pseudoatom
from pydesc.chemistry.base import register_dynamic_feature
from pydesc.chemistry.base import register_pseudoatom
from pydesc.chemistry.factories import BioPythonAtomSetFactory
from pydesc.chemistry.full_atom import MonoatomicIon
from pydesc.chemistry.full_atom import Nucleotide
from pydesc.chemistry.full_atom import Residue
from pydesc.chemistry.full_atom import calculate_residues_angles_vectorized

TYPE_THRESHOLDS = {Nucleotide: 0.25, Residue: 0.01, MonoatomicIon: 0.0}


@pytest.mark.parametrize(
    "superclass, expected_is_chainable",
    ((Mer, True), (AtomSet, False), (Ligand, False),),
)
def test_mer_subclass(superclass, expected_is_chainable):
    class Subclass(superclass):
        def __init__(self, ind, name, chain, atoms):
            super().__init__(ind, name, chain, atoms)

        def is_new(self):
            return True

        @register_pseudoatom
        def my_pa(self):
            return Pseudoatom(0, 0, 0)

        @register_dynamic_feature
        def my_df(self):
            return 42, 3.14

    atoms = {
        "A": Atom(numpy.array([1.0, 1.0, 1.0]), "A"),
        "B": Atom(numpy.array([0.0, 0.0, 0.0]), "B"),
    }
    mer = Subclass(42, "42", "C", atoms)

    numpy.testing.assert_array_equal(mer.gc.vector, (0.5, 0.5, 0.5))
    assert mer.is_new()
    assert mer.is_chainable() == expected_is_chainable
    assert not mer.dynamic_features
    assert mer.my_df == (42, 3.14)
    assert "my_df" in mer.dynamic_features
    assert isinstance(mer.my_pa, Pseudoatom)


class TestAtomSetFactory:
    def test_default_mer_factory_create_from_pdb_res(
        self, structure_file_w_pure_type, structures_dir, pure_types_2_mers
    ):
        root, fname = os.path.split(structure_file_w_pure_type)
        dummy, type_ = os.path.split(root)
        factory = BioPythonAtomSetFactory()
        pth = os.path.join(structures_dir, structure_file_w_pure_type)
        pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(fname, pth)

        expected_type = pure_types_2_mers[type_]
        type_threshold = TYPE_THRESHOLDS[expected_type]

        points = {True: 0, False: 1}
        length = 0

        for model in pdb_structure:
            for chain in model:
                for residue in chain:
                    result, warns = factory.create(residue, warn_in_place=False)

                    if result is None:
                        assert warns is None
                        assert residue.get_id()[0] == "W"
                    else:
                        length += 1
                        points[expected_type in result] += 1

        success_rate = (length - points[True]) / float(length)
        msg = "%i out of %i mers are of expected type" % (points[True], length)
        assert success_rate <= type_threshold, msg

    def test_create_residue_from_pdb_res(self, structures_dir):
        factory = BioPythonAtomSetFactory()
        pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(
            "5MPV.pdb", os.path.join(structures_dir, "prots_only", "5MPV.pdb")
        )

        r1 = pdb_structure[0]["D"][15]
        r2 = pdb_structure[0]["D"][16]
        r3 = pdb_structure[0]["D"][17]
        prev = factory.create(r1, warn_in_place=False)[0][Residue]
        res = factory.create(r2, warn_in_place=False)[0][Residue]
        next = factory.create(r3, warn_in_place=False)[0][Residue]

        res.prev_mer = prev
        res.next_mer = next

        assert len(tuple(res.iter_bb_atoms())) == 3
        assert len(tuple(res.iter_nbb_atoms())) > 0

        assert res.angles
        assert res.cbx
        assert res.rc
        assert res.backbone_average
        assert res.ca

        assert res.angles == pytest.approx((-0.62, -1.30), abs=0.01)
        assert res.next_mer is next
        assert res.prev_mer is prev


class AtomSetTest:
    cache = {}

    @classmethod
    def iter_structure(cls, file_, class_, test_structures_dir):
        factory = BioPythonAtomSetFactory()
        try:
            pdb_structure = cls.cache[(file_, test_structures_dir)]
        except KeyError:
            pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(
                "test", os.path.join(test_structures_dir, file_)
            )
            cls.cache[(file_, test_structures_dir)] = pdb_structure
        to_return = []
        for model in pdb_structure:
            for chain in model:
                for residue in chain:
                    result, warns = factory.create(residue, warn_in_place=False)
                    if result is None or class_ not in result:
                        continue
                    result = result[class_]
                    to_return.append(result)
        return to_return


class TestResidue(AtomSetTest):
    def test_calculate_cbx(self, protein_file, structures_dir):
        checked = 0
        for result in self.iter_structure(protein_file, Residue, structures_dir):
            if result.name == "GLY":
                distance = (result.cbx - result.atoms["CA"]).calculate_length()
                assert distance < 5
                continue
            assert result.cbx
            del result.pseudoatoms["cbx"]
            assert result.cbx
            length = (result.atoms["CB"] - result.cbx).calculate_length()
            assert round(length, 2) == 1.0
            unit_a_bx = (result.cbx - result.atoms["CA"]).get_unit_vector()
            unit_a_b = (result.atoms["CB"] - result.atoms["CA"]).get_unit_vector()
            dot_prod = unit_a_b.dot(unit_a_bx)
            assert pytest.approx(dot_prod) == 1.0
            checked += 1
        assert checked != 0

    def test_rc(self, protein_file, structures_dir):
        for result in self.iter_structure(protein_file, Residue, structures_dir):
            rc_atom = result.rc
            del result.pseudoatoms["rc"]
            assert tuple(rc_atom) == tuple(result.rc)
            for atom in result:
                assert (result.rc - atom).calculate_length() < 8

    def test_angles(self, protein_file, structures_dir):
        mers = self.iter_structure(protein_file, Residue, structures_dir)
        for m1, m2 in zip(mers, mers[1:]):
            if m1.has_bond(m2):
                m1.next_mer = m2
                m2.prev_mer = m1

        calculate_residues_angles_vectorized(mers)
        for mer in mers:
            assert "angles" in mer.dynamic_features
            angs = mer.dynamic_features["angles"]
            del mer.dynamic_features["angles"]
            psi1, phi1 = angs
            psi2, phi2 = mer.angles
            assert pytest.approx(psi1, abs=1.0e-4) == psi2
            assert pytest.approx(phi1) == phi2

    def test_bb_average(self, protein_file, structures_dir):
        mers = self.iter_structure(protein_file, Residue, structures_dir)
        for m1, m2 in zip(mers, mers[1:]):
            if m1.has_bond(m2):
                m1.next_mer = m2
                m2.prev_mer = m1

        checked = 0
        for mer in mers:
            assert (mer.backbone_average - mer.ca).calculate_length() < 8
            a_l = mer.get_adjusted_length()
            if a_l is None:
                continue
            checked += 1
            assert a_l < 5
        assert checked > 0

    def test_basic(self, protein_file, structures_dir):
        mers = self.iter_structure(protein_file, Residue, structures_dir)
        for mer in mers:
            assert mer.is_chainable()
            assert isinstance(mer.rc, Pseudoatom)
            all_atoms = set(mer)
            bb_atoms = set(mer.iter_bb_atoms())
            nbb_atoms = set(mer.iter_nbb_atoms())
            assert all_atoms == (nbb_atoms | bb_atoms)


class TestNucleotide(AtomSetTest):
    def test_calculate_features(self, nuclei_file, structures_dir):
        checked = 0
        for result in self.iter_structure(nuclei_file, Nucleotide, structures_dir):
            assert result.nx
            assert result.ring_center
            assert result.prc
            for pseudoatom in result.nx, result.ring_center, result.prc:
                dist_from_rc = (result.rc - pseudoatom).calculate_length()
                assert dist_from_rc < 10
            checked += 1
        assert checked != 0


class TestIon(AtomSetTest):
    def test_basic(self, ligands_structures_path):
        mers = self.iter_structure("2SRT.pdb", MonoatomicIon, ligands_structures_path)
        for mer in mers:
            assert not mer.is_chainable()
            assert isinstance(mer.rc, Pseudoatom)
            assert (tuple(mer)[0].vector == mer.rc.vector).all()
            assert mer.get_radius() < 10


class TestLigand(AtomSetTest):
    def test_basic(self, ligands_structures_path):
        mers = self.iter_structure("2SRT.pdb", MonoatomicIon, ligands_structures_path)
        for mer in mers:
            assert not mer.is_chainable()
            assert isinstance(mer.rc, Pseudoatom)
            for at in mer:
                assert (at - mer.rc).calculate_length() < 10
