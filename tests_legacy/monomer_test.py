# Copyright 2017 Pawel Daniluk, Tymoteusz Oleniecki
#
# This file is part of PyDesc.
#
# PyDesc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PyDesc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.

"""
Unit tests_legacy for mers.py.

Usage:
    python mers_test.py [-v] [--fast]

    or

    python -m unittest [-v] monomer_test

Pawel, Tymoteusz
"""

import unittest
import mock

from tests_legacy.syntax_check import module_syntax
from tests_legacy.syntax_check import parse_args

import sys
import types
import itertools
import numpy
import random
import os

import Bio.PDB

from argparse import ArgumentParser as AP

import pydesc.structure as structure
import pydesc.mers as monomer
import pydesc.config as config
import tests_legacy
from pydesc.warnexcept import IncompleteChainableParticle
from pydesc.warnexcept import set_filters

config.ConfigManager.warnings_and_exceptions.class_filters.set("UnknownParticleName", "ignore")
config.ConfigManager.warnings_and_exceptions.class_filters.set("IncompleteChainableParticle", "always")
set_filters()

data_dir = os.path.join(os.path.abspath(os.path.dirname(tests_legacy.__file__)), 'data/test_structures/')

# pylint: disable=C0111


def make_monomertestinits(strname, test_structure_dir, cls, names, crucial_atoms):
    """ Create and return a MonomerTestBasic test case for a given structure and monomer class. """

    class MonomerTestInit(unittest.TestCase):

        @classmethod
        def setUpClass(clss):
            clss.mf = monomer.MonomerFactory()

        def setUp(self):
            self.pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(strname, os.path.join(data_dir, test_structure_dir, strname))
            self.stc = mock.MagicMock()
            self.stc._monomers = []

        def test_init_random_name(self):
            errors = {}
            for i in range(100):
                nm = random.choice(names)
                atoms = {i: monomer.Atom(numpy.random.randn(3), 'E') for i in crucial_atoms}
                try:
                    it = cls(self.stc, i, nm, 'A', atoms)
                except Exception as e:
                    errors[i] = e
            self.assertEqual(len(errors), 0, "Failed to create %i mers out of 100" % len(errors))

        def test_init_from_MF_unpacked_Bio(self):
            errors = {}
            for n, res in enumerate(self.pdb_structure.get_residues()):
                data = self.mf.unpack_pdb_residue(res)
                if 'HOH' == data[0]:
                    continue
                try:
                    it = cls(self.stc, n, *data)
                except Exception as e:
                    errors[res] = e
            self.assertEqual(len(errors), 0, "Failed to create mers %s" % ",".join([str(i.id[1]) for i in errors]))

    return MonomerTestInit

def make_monomertestbasic(strname, test_structure_dir, cls, names):
    """ Create and return a MonomerTestBasic test case for a given structure and monomer class. """

    class MonomerTestBasic(unittest.TestCase):

        def test_get_config(self):
            config.ConfigManager.monomer.residue.set_default('test', 1)
            config.ConfigManager.monomer.nucleotide.set_default('test', 2)
            config.ConfigManager.monomer.ligand.set_default('test', 3)
            config.ConfigManager.monomer.ion.set_default('test', 4)
            test_dict = {monomer.Residue: 1, monomer.Nucleotide: 2, monomer.Ligand: 3, monomer.Ion: 4}
            for n, res in enumerate(self.pdb_structure.get_residues()):
                if cls == monomer.MonomerOther:
                    cls = res.__class__
                self.assertEqual(test_dict[cls], res.get_config('test'), "")

                class TestMonomer(cls):
                    pass

                artif = TestMonomer(res.pdb_residue)
                self.assertEqual(test_dict[cls], artif.get_config('test'))

        def test_select(self):
            nc = numcon.NumberConverter(self.pdb_structure)
            class S:
                converter = nc
            for i in self.ress:
                for res in i.values():
                    res.structure = S
                    pid = res.pdb_residue.get_full_id()
                    pids = numcon.PDB_id.create_from_string(pid[2] + str(pid[3][1]) + pid[3][2])
                    res.ind = nc.get_ind(pids)
                    s = res.select()
                    self.assertEquals(len(s.ids), 1)
                    self.assertEquals(pids, s.ids[0])

        def test_adjusted_length(self):
            ress = filter(lambda x: any(isinstance(key, monomer.Residue) for key in x), self.ress)
            for i in zip(ress[:-1], ress[1:]):
                m1 = i[0][monomer.Residue]
                m2 = i[1][monomer.Residue]
                m1.next_mer = m2
                m2.previous_mer = m1
            for i in self.ress:
                if i in ress:
                    self.assertTrue(i[monomer.Residue].adjusted_lenght() <= 3)
                else:
                    if monomer.Residue in i:
                        self.assertIsNone(i[monomer.Residue].adjusted_length())

    return MonomerTestBasic

def make_monomertest(strname, cls):
    """ Create and return a MonomerTest test case for a given structure and monomer class. """

    class MonomerTest(unittest.TestCase):
        cont_excl = {
            '1no5': [52], '1pxq': [], '2lp2': [], '3ftk': [], '3lgb': [],
            '3tk0': [], '1asz': [], '1gax': [], '3g88': [], '3m6x': [], '3npn': [38]}

        def setUp(self):
            self.load_mers(strname, data_dir + '%s.pdb' % strname, cls)

        def load_mers(self, name, filename, cls):
            self.pdb_structure = Bio.PDB.PDBParser(
                QUIET=True).get_structure(name, filename)

            self.pdbl = []
            self.l = []

            for res in self.pdb_structure.get_residues():
                res_id = res.get_id()[0]
                if not res_id.startswith('H_') and not res_id.startswith('W'):
                    try:
                        self.pdbl.append(res)
                        self.l.append(cls(res))
                    except KeyError as e:
                        self.assertTrue(str(e).strip() not in res, "Monomer %s doesn't lack %s, but programe failed to create it (struc: %s)" % (res.get_id(), str(e), self.pdb_structure.id))

        def test_indicators(self):
            for n in self.l:
                for name in n.indicators:
                    self.assertIsInstance(getattr(n, name), geometry.Coord)

                for c in n.representation:
                    self.assertIsInstance(c, geometry.Coord)

        def test_seq(self):
            for n in self.l:
                c = n.seq
                s3 = n.seq_1to3(c)
                self.assertIsInstance(c, types.StringType)
                self.assertEqual(len(c), 1)
                self.assertEqual(len(s3), 3)

        def test_seq3(self):
            for n in self.l:
                c = n.seq3
                s1 = n.seq_3to1(c)
                self.assertIsInstance(c, types.StringType)
                self.assertEqual(len(c), 3)
                self.assertEqual(len(s1), 1)

        def test_backbone(self):
            for n in self.l:
                for c in n.backbone:
                    self.assertIsInstance(c, geometry.Coord)

        def test_calculate_rc(self):
            for m in self.l:
                try:
                    m.calculate_rc()
                except Exception as e:
                    self.fail("Test for %s borks (%s)." % (str(m), str(e)))

        def test_iterators(self):
            for (pdbr, r) in zip(self.pdbl, self.l):
                self.assertEqual(len(pdbr), len(list(r.iter_atoms())))

                for a in r.iter_atoms():
                    self.assertIn(a, r)

                for a in r.iter_atoms():
                    self.assertIn(a, r)

                for a in r.pseudoatoms.values():
                    self.assertIn(a, r)

                ll = list(r.pseudoatoms.values()) + list(r.atoms.values())
                for a in r:
                    self.assertIn(a, ll)

        @unittest.expectedFailure
        def test_getsetattr(self):
            for m in self.l:
                for aname in itertools.chain(iter(m.atoms), iter(m.pseudoatoms)):
                    self.assertEqual(aname, aname.strip())
                    #TO -- commented to avoid deprecation warning
                    #~ self.assertEqual(m[aname], getattr(m, aname.lower()))
                    #~ self.assertEqual(m[aname], getattr(m, aname.upper()))

        if issubclass(cls, monomer.Residue):
            def test_calculate_cbx(self):
                for n in self.l:
                    n.calculate_cbx()

                    if n.name == 'GLY':
                        self.assertAlmostEqual((n.cbx - n.atoms['CA ']).calculate_length(), 0)
                    else:
                        self.assertAlmostEqual((n.cbx - n.atoms['CB ']).calculate_length(), 1)

            @testing(monomer.Residue.ca)
            def test_ca(self):
                for n in self.l:
                    self.assertEqual(n.ca, n.atoms['CA '])

        if issubclass(cls, monomer.Nucleotide):
            def test_ring_geom(self):
                for n in self.l:
                    n.calculate_ring_plane()
                    n.calculate_ring_center()

                    pl = n.ring_plane
                    rc = n.ring_center

                    rcp = pl.ort_projection(rc)

                    self.assertLess((rcp - rc).calculate_length(), 0.01)

    return MonomerTest


def make_monomercreatetest(strname):
    """ Create and return a MonomerTestBasic test case for a given structure and monomer class. """

    class MonomerCreateTest(unittest.TestCase):

        def setUp(self):
            self.pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(strname, data_dir + '%s.pdb' % strname)

        def test_init(self):
            for res in self.pdb_structure.get_residues():
                res_id = res.get_id()[0]
                m = monomer.Monomer.create_monomers(res)
                if res_id.startswith('W'):
                    self.assertIsNone(m, "%s is water. No monomer should be created." % str(res.get_id()))
                    continue
                elif res_id.startswith('H_'):
                    #~ self.assertIsInstance(m, monomer.MonomerOther, "%s is hetero. MonomerOther should be created." % str(res.get_id()))
                    # since create_monomers has been changed and returns dictionary, I need to change the test:
                    self.assertTrue(any(isinstance(i, monomer.MonomerOther) for i in m.values()), "%s is hetero. MonomerOther should be created (struc %s)." % (str(res.get_id()), self.pdb_structure.id))
                else:

                    def hasitem(obj, name):
                        try:
                            obj[name]
                            return True
                        except KeyError:
                            return False
                    mch = [i for i in m.values() if isinstance(i, monomer.MonomerChainable)]
                    if len(mch) == 0:
                        for cls, additionally in [[monomer.Residue, ('CB',)], [monomer.Nucleotide, ('',)]]:
                            needed_atoms = cls.get_config('backbone_atoms') + additionally
                            if not all(i.strip() in res for i in needed_atoms):
                                continue
                            #~ self.assertIsInstance(m, monomer.MonomerChainable, "%s is regular. MonomerChainable should be created." % str(res.get_id()))
                            #TO -- same here:
                            self.assertTrue(False, "%s is regular. MonomerChainable should be created (struc %s)." % (str(res.get_id()), self.pdb_structure.id))

                    #~ self.assertEqual(len(m.atoms.values()), len(res), "%s %s" % (str(list(m)), str(list(res))))
                    #TO -- while m is not one monomer...
                    self.assertTrue(any(len(i.atoms) == len(res) for i in m.values()))
                    self.assertEqual(len(m), len(list(m)))

                for i in m.values():
                    tests_performed["monomer.%s" % i.__class__.__name__] = True

    return MonomerCreateTest

def make_monomerfactorytest(strname, test_structure_dir, cls_):

    class TestMonomerFactory(unittest.TestCase):

        @classmethod
        def setUpClass(cls):
            f = os.path.join(data_dir, test_structure_dir, strname)
            cls.pdbS = Bio.PDB.PDBParser(QUIET=True).get_structure(f, f)
            cls.cls = cls_

        def setUp(self):
            self.mf = monomer.MonomerFactory()

        def test_init(self):
            pass

        def test_create_residues_from_Bio(self):
            fails = []
            done = 0
            for ch in list(self.pdbS.get_models())[0]:
                for m in ch:
                    res = self.mf.create_from_biopdb(m)
                    mn = m.get_resname().strip()
                    if mn == 'HOH':
                        self.assertIs(res, None)
                        continue
                    if self.cls not in res:
                        fails.append((m, res))
                        continue
                    done += 1
                    for i in res.values():
                        self.assertEqual(i.name, mn)
            self.assertLess(len(fails), .1 * done)

    return TestMonomerFactory

def load_tests(loader, standard_tests, pattern):
    """ Add tests_legacy created by make_* functions for all structures. Return a complete TestSuite. """

    # selecting structures
    bio_stc = {}
    bio_stc_sh = {}
    for i in (('prots_only', monomer.Residue), ('rna_only', monomer.Nucleotide), ('dna_only', monomer.Nucleotide)):
        for f in os.listdir(os.path.join(data_dir, i[0])):
            bio_stc.setdefault(i, []).append(f)
    for i in bio_stc:
        bio_stc_sh[i] = bio_stc[i][:3]

    if cla.fast:
        bio_stc = bio_stc_sh

    factory = unittest.TestSuite()

    # mers names
    names = {monomer.Residue: ['ALA', 'GLY', 'TRP', 'TYR', 'LYS'],
             monomer.Nucleotide: ['U', 'A', 'C', 'G'],
             monomer.Nucleotide: ['dT', 'dA', 'dG', 'dC'],
             }

    nuc_ats = ("C1'", "P", "O5'", "C5'", "C4'", "C3'", "O3'", "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9")
    atoms = {monomer.Residue: ['CA', 'C', 'N', 'CB', 'CO'],
             monomer.Nucleotide: nuc_ats,
             monomer.Nucleotide: nuc_ats,
             }

    # making factory tests_legacy
    for pth, tp in bio_stc:
        for stc in bio_stc[(pth, tp)]:
            factory.addTests(loader.loadTestsFromTestCase(make_monomerfactorytest(stc, pth, tp)))

    # making init tests_legacy
    for pth, tp in bio_stc_sh:
        for stc in bio_stc_sh[(pth, tp)]:
            nms = names[tp]
            ats = atoms[tp]
            standard_tests.addTests(loader.loadTestsFromTestCase(make_monomertestinits(stc, pth, tp, nms, ats)))

    standard_tests.addTests(loader.loadTestsFromTestCase(module_syntax(monomer)))

    standard_tests.addTests(factory)

    return standard_tests


#~ class AtomEnhancements(unittest.TestCase):

    #~ """ TestCase for enhancements requested in monomer.Atom. """

    #~ @unittest.expectedFailure
    #~ def test_deprec(self):
        #~ """ Deprecation of __getitem__ interface to Atom. """
        #~ a = monomer.Atom(*randomcoord())
        #~ syntax_check.test_deprec(self, lambda: a['name'], "__getitem__")
        #~ syntax_check.test_deprec(self, lambda: a.__setitem__('name', 'aa'), "__setitem__")

    #~ @unittest.expectedFailure
    #~ def test_hashable(self):
        #~ """ Check whether Atom is hashable. """
        #~ a = monomer.Atom(*randomcoord())
        #~ try:
            #~ set([a])
        #~ except:
            #~ self.fail()

    #~ def test_pseudoatom_init(self):
        #~ """Checks if Pseudoatoms are created both ways."""
        #~ for d in range(100):
            #~ random_coords = randomcoord()
            #~ pa = monomer.Pseudoatom(*random_coords)
            #~ pa2 = monomer.Pseudoatom(numpy_vec=numpy.array(random_coords))
            #~ for i in zip(pa.vector, pa2.vector):
                #~ self.assertEqual(*i)


class MonomerClassMethods(unittest.TestCase):

    def test_seqs(self):
        cls = monomer.Residue
        self.assertEqual(cls.seq_3to1('ALA'), 'A')
        self.assertEqual(cls.seq_1to3('A'), 'ALA')
        self.assertEqual(cls.seq_3to1('CYS'), 'C')
        self.assertEqual(cls.seq_1to3('C'), 'CYS')
        self.assertEqual(cls.seq_3to1('ASP'), 'D')
        self.assertEqual(cls.seq_1to3('D'), 'ASP')
        self.assertEqual(cls.seq_3to1('GLU'), 'E')
        self.assertEqual(cls.seq_1to3('E'), 'GLU')
        self.assertEqual(cls.seq_3to1('PHE'), 'F')
        self.assertEqual(cls.seq_1to3('F'), 'PHE')
        self.assertEqual(cls.seq_3to1('GLY'), 'G')
        self.assertEqual(cls.seq_1to3('G'), 'GLY')
        self.assertEqual(cls.seq_3to1('HIS'), 'H')
        self.assertEqual(cls.seq_1to3('H'), 'HIS')
        self.assertEqual(cls.seq_3to1('ILE'), 'I')
        self.assertEqual(cls.seq_1to3('I'), 'ILE')
        self.assertEqual(cls.seq_3to1('LYS'), 'K')
        self.assertEqual(cls.seq_1to3('K'), 'LYS')
        self.assertEqual(cls.seq_3to1('LEU'), 'L')
        self.assertEqual(cls.seq_1to3('L'), 'LEU')
        self.assertEqual(cls.seq_3to1('MET'), 'M')
        self.assertEqual(cls.seq_1to3('M'), 'MET')
        self.assertEqual(cls.seq_3to1('ASN'), 'N')
        self.assertEqual(cls.seq_1to3('N'), 'ASN')
        self.assertEqual(cls.seq_3to1('PRO'), 'P')
        self.assertEqual(cls.seq_1to3('P'), 'PRO')
        self.assertEqual(cls.seq_3to1('GLN'), 'Q')
        self.assertEqual(cls.seq_1to3('Q'), 'GLN')
        self.assertEqual(cls.seq_3to1('ARG'), 'R')
        self.assertEqual(cls.seq_1to3('R'), 'ARG')
        self.assertEqual(cls.seq_3to1('SER'), 'S')
        self.assertEqual(cls.seq_1to3('S'), 'SER')
        self.assertEqual(cls.seq_3to1('THR'), 'T')
        self.assertEqual(cls.seq_1to3('T'), 'THR')
        self.assertEqual(cls.seq_3to1('TRP'), 'W')
        self.assertEqual(cls.seq_1to3('W'), 'TRP')
        self.assertEqual(cls.seq_3to1('TYR'), 'Y')
        self.assertEqual(cls.seq_1to3('Y'), 'TYR')
        self.assertEqual(cls.seq_3to1('VAL'), 'V')
        self.assertEqual(cls.seq_1to3('V'), 'VAL')
        cls = monomer.Nucleotide
        self.assertEqual(cls.seq_3to1('  A'), 'A')
        self.assertEqual(cls.seq_1to3('A'), '  A')
        self.assertEqual(cls.seq_3to1('  U'), 'U')
        self.assertEqual(cls.seq_1to3('U'), '  U')
        self.assertEqual(cls.seq_3to1('  G'), 'G')
        self.assertEqual(cls.seq_1to3('G'), '  G')
        self.assertEqual(cls.seq_3to1('  C'), 'C')
        self.assertEqual(cls.seq_1to3('C'), '  C')
        self.assertEqual(cls.seq_3to1(' DA'), 'A')
        self.assertEqual(cls.seq_3to1(' DT'), 'T')
        self.assertEqual(cls.seq_1to3('T'), ' DT')
        self.assertEqual(cls.seq_3to1(' DG'), 'G')
        self.assertEqual(cls.seq_3to1(' DC'), 'C')


class PatchStc(object):

    def __init__(self):
        self._monomers = []


if __name__ == '__main__':
    cla, sys.argv[1:] = parse_args()

    unittest.main()
