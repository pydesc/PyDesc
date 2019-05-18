# Copyright 2017 Tymoteusz Oleniecki
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
Unit tests_legacy for contacts.py.

Usage:
    python cmap_test.py [-v] [--fast]

    or

    python -m unittest [-v] contacts_test

Tymoteusz
"""

import unittest
import random
import inspect

import syntax_check
from syntax_check import notest, testing, test, test_name_append

import Bio.PDB

import warnings

import pydesc.structure as structure
import pydesc.contacts as contacts
import pydesc.config as config
import pydesc.numberconverter as numberconverter
import tests_legacy
from pydesc.warnexcept import WrongMonomerType

config.ConfigManager.warnings_and_exceptions.class_filters.set("LocalCopyAccess", "ignore")

syntax_check.module = contacts

fast = False

TestSyntax = syntax_check.module_syntax()

test(contacts.ContactCriterion)
test(contacts.DihedralAngleCriterion)
test(contacts.DihedralAngleCriterion.__init__)
test(contacts.DihedralAngleCriterion.is_in_contact)
test(contacts.HorizontalBisectorDistanceCriterion)
test(contacts.HorizontalBisectorDistanceCriterion.__init__)
test(contacts.HorizontalBisectorDistanceCriterion.is_in_contact)
test(contacts.PointsDistanceCriterion)
test(contacts.PointsDistanceCriterion.__init__)
test(contacts.SetDistanceCriterion)
test(contacts.SetDistanceCriterion.__init__)
test(contacts.SetDistanceCriterion.is_in_contact)
test(contacts.VerticalBisectorDistanceCriterion)
test(contacts.VerticalBisectorDistanceCriterion.__init__)
test(contacts.VerticalBisectorDistanceCriterion.is_in_contact)
test(contacts.CombinedContact)
test(contacts.CombinedContact.__init__)

notest(contacts.ContactsExclusiveDisjunction)

#~ notest(contacts.ContactMap)

data_dir = tests_legacy.__path__[0] + '/data/test_structures/'

# pylint: disable=C0111,R0912


def make_contactbasictest(strname, crit_cls):
    """Create and return a ContactBasicTest testcase for a given structure and contact criterion. """

    @testing(crit_cls)
    @test_name_append(crit_cls, strname)
    class ContactBasicTest(unittest.TestCase):
        name = strname

        @classmethod
        def setUpClass(cls):
            cls.pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(cls.name, data_dir + '%s.pdb' % cls.name)
            cls.converter = numberconverter.NumberConverter(cls.pdb_structure)
            cls.model = random.choice(cls.pdb_structure)
            with warnings.catch_warnings(record=True):
                cls.struct = structure.Structure(cls.model, cls.converter)

        @testing(crit_cls.__init__)
        @testing(contacts.ContactCriterion.criteria)
        def test_init(self):
            crit_cls()

        @testing(getattr(crit_cls, 'is_in_contact', None))
        @testing(contacts.for_monomer_type_only)
        def test_is_in_contact(self):
            if fast:
                return
            with warnings.catch_warnings(record=True):
                for m1 in self.struct:
                    for m2 in self.struct:
                        try:
                            self.assertEqual(crit_cls().is_in_contact(m1, m2), crit_cls().is_in_contact(m2, m1))
                        except:
                            self.assertRaises(WrongMonomerType)

        @classmethod
        def tearDownClass(cls):
            del cls.struct
            del cls.model
            del cls.converter
            del cls.pdb_structure


    return ContactBasicTest


def make_combcontactbasictest(crit_cls):
    """Create and return a CombContactBasicTest testcase for a given structure and contact criterion. """

    class TestCrit(contacts.ContactCriterion):
        call_count = 0

        def __init__(self, val):
            self.bound = val

        def is_in_contact(self, m1, m2, lazy = 'dummy'):
            TestCrit.call_count += 1
            if m1 == m2:
                if m1 == self.bound:
                    return 1
                elif m1 > self.bound:
                    return 2
            return 0

        def _is_in_contact(self, m1, m2, lazy = 'dummy'):
            return self.is_in_contact(m1, m2, lazy)


    class TestCrit1(contacts.ContactCriterion):
        call_count = 0

        def __init__(self, v, w):
            self.v = v
            self.w = w

        def is_in_contact(self, m1, m2, lazy = 'dummy'):
            if (m1 == self.v and m2 == self.w) or (m1 == self.w and m2 == self.v):
                return 2
            return 0

        def _is_in_contact(self, m1, m2, lazy = 'dummy'):
            return self.is_in_contact(m1, m2, lazy)

    @testing(crit_cls)
    @test_name_append(crit_cls)
    class CombContactBasicTest(unittest.TestCase):

        @testing(contacts.CombinedContact.criteria)
        @testing(crit_cls.__init__)
        def test_init(self):
            crit_cls(*map(TestCrit, range(9)))

        @testing(contacts.ContactCriterion.__eq__)
        def test_CC_eq(self):
            a = TestCrit1(1,2)
            b = TestCrit1(1,2)
            self.assertEqual(a, b)
            c = TestCrit1(1,3)
            self.assertNotEqual(a, c)
            self.assertNotEqual(b, c)
            d = contacts.ContactsAlternative(a, c)
            e = contacts.ContactsAlternative(a, c)
            self.assertEqual(d, e)
            f = contacts.ContactsAlternative(a, b)
            self.assertNotEqual(e, f)
            g = contacts.ContactsConjunction(a, b)
            self.assertNotEqual(f, g)

        @testing(getattr(crit_cls, 'is_in_contact', None))
        def test_is_in_contact(self):
            crit = crit_cls(*map(TestCrit, range(1, 10)))

            for i in range(10):
                TestCrit.call_count = 0
                res = crit.is_in_contact(i, i, lazy=False)

                if crit_cls == contacts.ContactsAlternative:
                    if i == 0:
                        exp_res = 0
                    elif i == 1:
                        exp_res = 1
                    else:
                        exp_res = 2
                elif crit_cls == contacts.ContactsConjunction:
                    if i == 10:
                        exp_res = 2
                    elif i == 9:
                        exp_res = 1
                    else:
                        exp_res = 0

                exp_count = 9

                self.assertEqual(res, exp_res, "res: %d exp_res: %d i: %d" % (res, exp_res, i))
                self.assertEqual(TestCrit.call_count, exp_count, "call_count: %d exp_count: %d i: %d" % (TestCrit.call_count, exp_count, i))

        @testing(getattr(crit_cls, 'is_in_contact', None))
        def test_is_in_contact_lazy(self):
            crit = crit_cls(*map(TestCrit, range(1, 10)))

            for i in range(10):
                TestCrit.call_count = 0
                res = crit.is_in_contact(i, i, lazy=True)

                if crit_cls == contacts.ContactsAlternative:
                    if i == 0:
                        exp_res = 0
                        exp_count = 9
                    elif i == 1:
                        exp_res = 1
                        exp_count = 9
                    else:
                        exp_res = 2
                        exp_count = 1
                elif crit_cls == contacts.ContactsConjunction:
                    if i == 10:
                        exp_res = 2
                        exp_count = 9
                    elif i == 9:
                        exp_res = 1
                        exp_count = 9
                    else:
                        exp_res = 0
                        exp_count = i + 1

                self.assertEqual(res, exp_res, "res: %d exp_res: %d i: %d" % (res, exp_res, i))
                self.assertEqual(TestCrit.call_count, exp_count, "call_count: %d exp_count: %d i: %d" % (TestCrit.call_count, exp_count, i))

        @testing(getattr(crit_cls, 'is_in_contact', None))
        def test_is_in_contact_lazy1(self):
            crit = crit_cls(*map(TestCrit, reversed(range(1, 10))))

            for i in range(10):
                TestCrit.call_count = 0
                res = crit.is_in_contact(i, i, lazy=True)

                if crit_cls == contacts.ContactsAlternative:
                    if i == 0:
                        exp_res = 0
                        exp_count = 9
                    elif i == 1:
                        exp_res = 1
                        exp_count = 9
                    else:
                        exp_res = 2
                        exp_count = 11 - i
                elif crit_cls == contacts.ContactsConjunction:
                    if i == 10:
                        exp_res = 2
                        exp_count = 9
                    elif i == 9:
                        exp_res = 1
                        exp_count = 9
                    else:
                        exp_res = 0
                        exp_count = 1

                self.assertEqual(res, exp_res, "res: %d exp_res: %d i: %d" % (res, exp_res, i))
                self.assertEqual(TestCrit.call_count, exp_count, "call_count: %d exp_count: %d i: %d" % (TestCrit.call_count, exp_count, i))

        if crit_cls == contacts.ContactsAlternative:
            @testing(contacts.ContactsAlternative.get_validating_sub_criterion)
            def test_get_validating_subcriterion(self):
                mers = [1,2,3,4,5,6]
                tc = {(1,2): TestCrit1(1,2), (3,4): TestCrit1(3,4), (5,6): TestCrit1(5,6)}
                crit = crit_cls(crit_cls(tc[(1,2)], tc[(3,4)]), tc[(5,6)])
                for i in mers:
                    for j in mers:
                        try:
                            sbc = crit.get_validating_sub_criterion(i, j)
                        except ValueError:
                            self.assertTrue((i, j) not in tc and (j, i) not in tc, "mers %i and %i has no contact while they should")
                        else:
                            self.assertTrue((i == sbc.w and j == sbc.v) or (j == sbc.w and i == sbc.v))

    return CombContactBasicTest


def make_desctest(str_name):
    @testing(contacts.DescriptorCriterion)
    class DescCrit(unittest.TestCase):
        @classmethod
        def setUpClass(cls):
            cls.s = structure.StructureLoader().load_structure(str_name)[0]
            cls.s.set_contact_map(contacts.RcContact())
            cls.ds = filter(bool, structure.AbstractDescriptor.create_descriptors(cls.s))

        @testing(contacts.DescriptorCriterion.__init__)
        @testing(contacts.DescriptorCriterion.is_in_contact)
        def test_DescriptorCriterion_is_in_contact(self):
            ms = list(self.s)
            for desc in self.ds[random.choice([1,2,3,4,5])::5]:
                crit = contacts.DescriptorCriterion(desc)
                cm = []
                for c in desc.contacts:
                    cm.append(c.elements[0].central_monomer)
                    cm.append(c.elements[1].central_monomer)
                for m1 in ms:
                    for m2 in ms:
                        if m1 == m2:
                            continue
                        val = crit.is_in_contact(m1, m2)
                        if any(i == desc.central_element.central_monomer for i in (m1, m2)) and all(i in cm for i in (m1, m2)):
                            cmv = self.s.contact_map.get_contact_value(m1, m2)
                            self.assertEqual(cmv, val, "%s and %s Contact value is %i, while value in contact map is %i" % (str(m1), str(m2), val, cmv))
                        else:
                            self.assertEqual(val, 0)

        @classmethod
        def tearDownClass(cls):
            del cls.ds
            del cls.s
    return DescCrit


def load_tests(loader, standard_tests, pattern):
    """ Add tests_legacy created by make_* functions for all structures. Returns a complete TestSuite. """

    structures = ['1asz', '1gax', '1no5', '1pxq', '2dlc', '2lp2',
                  '3ftk', '3g88', '3lgb', '3m6x', '3npn', '3tk0', '3umy']

    if fast:
        structures = ['3npn', '1no5']

    def is_non_abstract_crit(c):
        if isinstance(c, type) and c.__module__ == 'pydesc.contacts' and issubclass(c, contacts.ContactCriterion) and len(c.__abstractmethods__) == 0:
            aspec = inspect.getargspec(c.__init__)
            if len(aspec.args) == 1 and aspec.varargs is None and aspec.keywords is None:
                return True

        return False

    crit_classes = filter(is_non_abstract_crit, contacts.__dict__.values())

    basic = unittest.TestSuite()

    for name in structures:
        for crit in crit_classes:
            basic.addTests(loader.loadTestsFromTestCase(make_contactbasictest(name, crit)))
        basic.addTests(loader.loadTestsFromTestCase(make_contactbasictest(name, contacts.CaCbxContact)))
        basic.addTests(loader.loadTestsFromTestCase(make_desctest(name)))
    test(contacts.CaCbxSubtractionCriterion)
    test(contacts.CaCbxSubtractionCriterion.__init__)
    test(contacts.CaCbxSubtractionCriterion.is_in_contact)
    # TO these are tested during tests_legacy added above while CaCbxContact finction calls CaCbxSubtractionCriterion methods
    test(contacts.VectorDistanceCriterion)
    test(contacts.VectorDistanceCriterion.__init__)
    # TO tested for subclass


    basic.addTests(loader.loadTestsFromTestCase(make_combcontactbasictest(contacts.ContactsAlternative)))
    basic.addTests(loader.loadTestsFromTestCase(make_combcontactbasictest(contacts.ContactsConjunction)))

    standard_tests.addTests(basic)

    return standard_tests

if __name__ == '__main__':
    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
