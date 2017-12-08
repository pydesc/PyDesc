# Copyright 2017 Pawel Daniluk
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
Unit tests for compat/dsc.py.

Usage:
    python compat_dsc_test.py [-v] [--fast]

    or

    python -m unittest [-v] compat_dsc_test

created: 05.05.2014, Pawel Daniluk
"""

import unittest
import random

import itertools

import os.path

from collections import defaultdict

import tests.syntax_check as syntax_check
from tests.syntax_check import testing, test_name_append, notest

import Bio.PDB

import warnings

import pydesc.compat.dsc as dsc
import pydesc.structure as structure
import tests

syntax_check.module = dsc

notest(dsc.split_num)

data_dir = tests.__path__[0] + '/data/dsc/'

fast = False

TestSyntax = syntax_check.module_syntax()


# pylint: disable=C0111, R0201, R0912, R0914

def assert_points_equal(self, p1, p2):
    self.assertAlmostEqual(p1.x, p2.x, places=3)
    self.assertAlmostEqual(p1.y, p2.y, places=3)
    self.assertAlmostEqual(p1.z, p2.z, places=3)


def assert_monomers_equal(self, m, cm):
    self.assertEqual(cm.ind, m.ind)
    self.assertEqual(cm.type, cm.types_dict[type(m)])
    self.assertEqual(cm.type_name, type(m).__name__)

    if m.next_monomer:
        self.assertEqual(cm.next_ind, m.next_monomer.ind)
    else:
        self.assertEqual(cm.next_ind, 0)

    self.assertEqual(cm.n_points, len(m.indicators))

    for cp, p in zip(list(cm.points[0:cm.n_points]), m.representation):
        assert_points_equal(self, cp, p)

    for cpn, pn in zip(list(cm.point_names[0:cm.n_points]), m.indicators):
        self.assertEqual(cpn, pn)  # broken in trunk see #43


def load_structure_to_class(cls, strname):
    if strname in struct_dict:
        cls.pdb_structure, cls.converter, cls.model, cls.struct = struct_dict[strname]
    else:
        cls.pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(strname, data_dir + '%s.pdb' % strname)
        cls.converter = numberconverter.NumberConverter(cls.pdb_structure)
        cls.model = random.choice(cls.pdb_structure)
        with warnings.catch_warnings(record=True):
            cls.struct = structure.Structure(cls.model, cls.converter)

        struct_dict[strname] = (cls.pdb_structure, cls.converter, cls.model, cls.struct)


def make_dscfiletest(name):
    """ Create and return a DSCFileTest test case for a given DSC file. """

    @testing(dsc.DSCFile)
    @test_name_append(name)
    class DSCFileTest(unittest.TestCase):

        @testing(dsc.DSCFile.__init__)
        def test_init(self):
            with open(os.path.join(data_dir, name), 'r') as file_obj:
                dscfile = dsc.DSCFile(file_obj)

                contact_tuples = list(itertools.chain(*dscfile.contacts.values()))

                self.assertEqual(len(dscfile.tuples), len(dscfile.main) + len(contact_tuples))

                for tuple_ in dscfile.tuples:
                    if len(tuple_) in [3, 4]:
                        self.assertEqual(len(tuple_[2]), 5)
                        self.assertFalse(' ' in tuple_[2])
                    elif len(tuple_) == 14:
                        self.assertEqual(len(tuple_[3]), 5)
                        self.assertFalse(' ' in tuple_[3])
                    else:
                        self.fail("Wrong tuple length: %s" % str(tuple_))

                for num, list_ in dscfile.contacts.items():
                    for tuple_ in list_:
                        self.assertEqual(num, int(tuple_[0]))

                for tuple_ in dscfile.main:
                    try:
                        dummy = int(tuple_[0])
                        dummy = map(float, tuple_[5:])
                    except:
                        self.fail("Conversions fail: %s" % str(tuple_))

    return DSCFileTest


def make_dscnumberconvertertest(name):
    """ Create and return a DSCNumberConverterTest test case for a given DSC file. """

    @testing(dsc.DSCNumberConverter)
    @test_name_append(name)
    class DSCNumberConverterTest(unittest.TestCase):

        @classmethod
        def setUpClass(cls):
            with open(os.path.join(data_dir, name), 'r') as file_obj:
                cls.dscfile = dsc.DSCFile(file_obj)

        @testing(dsc.DSCNumberConverter.__init__)
        def test_init(self):
            dscnumberconverter = dsc.DSCNumberConverter(self.dscfile)

            self.assertEqual(len(dscnumberconverter.dict_ind_to_pdb), len(self.dscfile.main))
            self.assertEqual(len(dscnumberconverter.dict_pdb_to_ind), len(self.dscfile.main))

            for tuple_ in self.dscfile.main:
                num = int(tuple_[0])

                pdb_id = dscnumberconverter.get_pdb_id(num)
                num1 = dscnumberconverter.get_ind(pdb_id)

                self.assertEqual(num, num1)

        @classmethod
        def tearDownClass(cls):
            del cls.dscfile

    return DSCNumberConverterTest


def make_dscstructuretest(name):
    """ Create and return a DSCStructureTest test case for a given DSC file. """

    @testing(dsc.DSCStructure)
    @testing(dsc.DSCResidue)
    @testing(dsc.DSCNucleotide)
    @testing(dsc.DSCIon)
    @testing(dsc.DSCMonomer)
    @test_name_append(name)
    class DSCStructureTest(unittest.TestCase):

        @classmethod
        def setUpClass(cls):
            with open(os.path.join(data_dir, name), 'r') as file_obj:
                cls.dscfile = dsc.DSCFile(file_obj)

        @testing(dsc.DSCResidue.__init__)
        def test_residue_init(self):
            for tuple_ in self.dscfile.main:
                dscresidue = dsc.DSCResidue(tuple_)

                indicators = dscresidue.indicators
                representation = dscresidue.representation

                self.assertEqual(len(indicators), len(representation))

                self.assertEqual(dscresidue.ind, int(tuple_[0]))

        @testing(dsc.DSCNucleotide.__init__)
        def test_nucleotide_init(self):
            for tuple_ in self.dscfile.main:
                dscnucleotide = dsc.DSCNucleotide(tuple_)

                indicators = dscnucleotide.indicators
                representation = dscnucleotide.representation

                self.assertEqual(len(indicators), len(representation))

                self.assertEqual(dscnucleotide.ind, int(tuple_[0]))

        @testing(dsc.DSCIon.__init__)
        def test_nucleotide_init(self):
            for tuple_ in self.dscfile.main:
                dscion = dsc.DSCIon(tuple_)

                indicators = dscion.indicators
                representation = dscion.representation

                self.assertEqual(len(indicators), len(representation))

                self.assertEqual(dscion.ind, int(tuple_[0]))

        @testing(dsc.DSCMonomer.create_dscmonomer)
        def test_monomer_create_dscmonomer(self):
            for tuple_ in self.dscfile.main:
                dscmonomer = dsc.DSCMonomer.create_dscmonomer(tuple_)
                types = {14: dsc.DSCResidue, 17: dsc.DSCNucleotide, 8: dsc.DSCIon}

                self.assertTrue(isinstance(dscmonomer, types[len(tuple_)]))

        @testing(dsc.DSCStructure.__init__)
        @testing(dsc.DSCContactCriterion.__init__)
        def test_structure_init(self):
            dscstructure = dsc.DSCStructure(self.dscfile)

            repr(dscstructure)

            self.assertEqual(len(dscstructure), len(self.dscfile.main))

            for mer in dscstructure:
                if mer.next_monomer:
                    self.assertEqual(mer, mer.next_monomer.previous_monomer)
                if mer.previous_monomer:
                    self.assertEqual(mer, mer.previous_monomer.next_monomer)

            for desc_ind in self.dscfile.contacts.keys():
                try:
                    dummy = structure.Element.build(dscstructure[desc_ind])
                except:
                    self.fail("Couldn't build an element around residue with contacts: %s" % str(self.dscfile.main[desc_ind]))

        @testing(dsc.DSCContactCriterion)
        @testing(dsc.DSCContactCriterion.is_in_contact)
        def test_structure_descriptors(self):
            dscstructure = dsc.DSCStructure(self.dscfile)

            for desc_ind in self.dscfile.contacts.keys():
                if len(self.dscfile.contacts[desc_ind]) > 0:
                    try:
                        desc = structure.AbstractDescriptor.build(structure.ElementChainable(dscstructure[desc_ind]))
                    except:
                        self.fail("Couldn't build a descriptor around residue with contacts: %s" % str(self.dscfile.main[desc_ind - 1]))

                    self.assertIsInstance(desc, structure.ProteinDescriptor)
                    self.assertEqual(len(desc.contacts), len(self.dscfile.contacts[desc_ind]))

                    for contact_obj in desc.contacts:
                        other_ind = contact_obj.get_other_element(desc.central_element).central_monomer.ind
                        tuple_ = [t for t in self.dscfile.contacts[desc_ind] if dsc.split_num(t[1])[0] == other_ind][0]

                        if tuple_[-1] == 'O':
                            self.assertEqual(contact_obj.value, 1)
                        else:
                            self.assertEqual(contact_obj.value, 2)

                    self.assertEqual(len(list(desc)), len(list(desc)))

        @testing(dsc.DSCStructure.select)
        def test_structure_select(self):
            dscstructure = dsc.DSCStructure(self.dscfile)
            selection_obj = dscstructure.select()

            spec_obj = selection_obj.specify(dscstructure)

            self.assertEqual(len(list(spec_obj)), len(dscstructure))

        @classmethod
        def tearDownClass(cls):
            del cls.dscfile

    return DSCStructureTest


def load_tests(loader, standard_tests, pattern):
    """ Add tests created by make_* functions for all structures. Return a complete TestSuite. """
    files = ['d1nkla_.dsc', 'd1qdma1.dsc']

    basic = unittest.TestSuite()

    for name in files:
        basic.addTests(loader.loadTestsFromTestCase(make_dscfiletest(name)))
        basic.addTests(loader.loadTestsFromTestCase(make_dscnumberconvertertest(name)))
        basic.addTests(loader.loadTestsFromTestCase(make_dscstructuretest(name)))

    standard_tests.addTests(basic)

    return standard_tests

if __name__ == '__main__':
    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
