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
Unit tests_legacy for cydesc/maps.py.

Usage:
    python cydesc_test.py [-v] [--fast]

    or

    python -m unittest [-v] cydesc_test

created: 16.04.2014, Pawel Daniluk
"""

import unittest
import random
import ctypes
import operator

from collections import defaultdict

import tests_legacy.syntax_check as syntax_check
from tests_legacy.syntax_check import testing, test_name_append, notest

import Bio.PDB

import warnings

import pydesc.structure as structure
import pydesc.numberconverter as numberconverter
import pydesc.cydesc as cydesc
import pydesc.contacts as contacts
import tests_legacy

# This is not a constant.
libcydesc_test = cydesc.load_library('cydesc_test')  # pylint: disable=C0103

libcydesc_test.CPointTest_array_cpy.restype = ctypes.POINTER(cydesc.CPoint)
libcydesc_test.CMerTest_array_cpy.restype = ctypes.POINTER(cydesc.CMer)
libcydesc_test.CSegTest_array_cpy.restype = ctypes.POINTER(cydesc.CSeg)

syntax_check.module = cydesc

notest(cydesc.load_library)

data_dir = tests_legacy.__path__[0] + '/data/test_structures/'

fast = False

TestSyntax = syntax_check.module_syntax()

# pylint: disable=C0111, R0201, R0912, R0914

struct_dict = {}


def assert_points_equal(self, p1, p2):
    self.assertAlmostEqual(p1.x, p2.x, places=3)
    self.assertAlmostEqual(p1.y, p2.y, places=3)
    self.assertAlmostEqual(p1.z, p2.z, places=3)


def assert_monomers_equal(self, m, cm):
    self.assertEqual(cm.ind, m.ind)
    self.assertEqual(cm.type, cm.types_dict[type(m)])
    self.assertEqual(cm.type_name, type(m).__name__)

    if m.next_mer:
        self.assertEqual(cm.next_ind, m.next_mer.ind)
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


def del_structure_from_class(cls):
    del cls.pdb_structure
    del cls.converter
    del cls.model
    del cls.struct


@testing(cydesc.CPoint)
class CPointTest(unittest.TestCase):

    @testing(cydesc.CPoint.__init__)
    def test_init(self):
        c = syntax_check.randomcoord()

        point = cydesc.CPoint(c)

        self.assertAlmostEqual(c[0], point.x, places=3)
        self.assertAlmostEqual(c[1], point.y, places=3)
        self.assertAlmostEqual(c[2], point.z, places=3)

        repr(point)

    @testing(cydesc.CPoint.__iter__)
    def test_iter(self):
        c = syntax_check.randomcoord()

        point = cydesc.CPoint(c)

        self.assertEqual((point.x, point.y, point.z), tuple(point))

    def test_del(self):
        point = cydesc.CPoint(syntax_check.randomcoord())
        del point

    def test_array(self):
        n = 10
        xl, yl, zl = [syntax_check.randomcoord(n) for dummy in range(3)]

        points = map(cydesc.CPoint, zip(xl, yl, zl))
        point_arr = (cydesc.CPoint * n)(*points)

        float_arr_br = lambda x: ctypes.byref((ctypes.c_float * n)(*x))

        res = libcydesc_test.CPointTest_array(n, float_arr_br(xl), float_arr_br(yl), float_arr_br(zl), point_arr)

        self.assertEquals(res, 0)

    def test_array_cpy(self):
        n = 10
        xl, yl, zl = [syntax_check.randomcoord(n) for dummy in range(3)]

        points = map(cydesc.CPoint, zip(xl, yl, zl))
        point_arr = (cydesc.CPoint * n)(*points)

        res_ptr = libcydesc_test.CPointTest_array_cpy(n, point_arr)

        self.assertEqual(map(tuple, points), map(tuple, res_ptr[0:n]))

        libcydesc_test.free_(res_ptr)


@testing(cydesc.CSeg)
class CSegTest(unittest.TestCase):

    @testing(cydesc.CSeg.__iter__)
    def test_array_cpy(self):
        n = 10
        sl, el = [random.sample(range(100), n) for dummy in range(2)]

        segs = map(cydesc.CSeg, sl, el)
        seg_arr = (cydesc.CSeg * n)(*segs)

        res_ptr = libcydesc_test.CSegTest_array_cpy(n, seg_arr)

        self.assertEqual(map(tuple, segs), map(tuple, res_ptr[0:n]))

        libcydesc_test.free_(res_ptr)


def make_cmertest(strname):
    """ Create and return a CMerTest test case for a given structure. """

    @testing(cydesc.CMer)
    @test_name_append(strname)
    class CMerTest(unittest.TestCase):

        @classmethod
        def setUpClass(cls):
            load_structure_to_class(cls, strname)

        @testing(cydesc.CMer.__init__)
        def test_init(self):
            for m in self.struct:
                cm = cydesc.CMer(m)

                assert_monomers_equal(self, m, cm)
                repr(cm)

        def test_array_cpy(self):
            n = len(self.struct)
            mer_arr = (cydesc.CMer * n)(*map(cydesc.CMer, self.struct))

            res_ptr = libcydesc_test.CMerTest_array_cpy(n, mer_arr)

            for m1, m2 in zip(mer_arr, res_ptr[0:n]):
                def attr_eq(name):
                    self.assertEqual(getattr(m1, name), getattr(m2, name))

                map(attr_eq, ['ind', 'type', 'type_name', 'next_ind', 'n_points'])

                for p1, p2 in zip(list(m1.points[0:m1.n_points]), list(m2.points[0:m2.n_points])):
                    assert_points_equal(self, p1, p2)

                for n1, n2 in zip(list(m1.point_names[0:m1.n_points]), list(m2.point_names[0:m2.n_points])):
                    self.assertEqual(n1, n2)

            libcydesc_test.free_(res_ptr)

        @classmethod
        def tearDownClass(cls):
            del_structure_from_class(cls)

    return CMerTest


def make_cstructuretest(strname):
    """ Create and return a CStructureTest test case for a given structure. """

    @testing(cydesc.CStructure)
    @test_name_append(strname)
    class CStructureTest(unittest.TestCase):

        @classmethod
        def setUpClass(cls):
            load_structure_to_class(cls, strname)

        @testing(cydesc.CStructure.__init__)
        def test_init(self):
            cs = cydesc.CStructure(self.struct)

            self.assertEqual(self.struct.name, cs.name)
            self.assertEqual(cs.n_monomers, len(self.struct))

            for m, cm in zip(self.struct, cs.monomers):
                assert_monomers_equal(self, m, cm)

            def renumber(seg):
                return [list(self.struct)[i].ind for i in list(seg)]

            sum_len = 0
            for seg in cs.segs[0:cs.n_segs]:
                start, end = seg
                frag = self.struct[start:end]
                sum_len += len(frag)

                if seg.start != seg.end:
                    self.assertIsInstance(frag, structure.Segment)

            self.assertEqual(len(self.struct), sum_len)

            for s1, s2 in zip(cs.segs[0:cs.n_segs], cs.segs[1:cs.n_segs]):
                end1 = s1.end
                start2 = s2.start
                frag = self.struct[end1:start2]
                self.assertIsInstance(frag, structure.PartialStructure)

            repr(cs)

        def test_indices(self):
            cs = cydesc.CStructure(self.struct)
            res = libcydesc_test.CStructureTest_indices(ctypes.byref(cs))
            self.assertEqual(res, 0)

        def test_adjusted_number(self):
            cs = cydesc.CStructure(self.struct)

            for dummy in range(100):
                start, end = sorted(map(operator.attrgetter('ind'), random.sample(list(self.struct), 2)))
                res = libcydesc_test.CStructureTest_adjusted_number(ctypes.byref(cs), start, end)

                self.assertEqual(res, self.struct[start:end].adjusted_number())

        @classmethod
        def tearDownClass(cls):
            del_structure_from_class(cls)

    return CStructureTest


def make_ccontactmaptest(strname):
    """ Create and return a CContactTest test case for a given structure. """

    @testing(cydesc.CContactMap)
    @testing(cydesc.CContact)
    @test_name_append(strname)
    class CContactMapTest(unittest.TestCase):

        @classmethod
        def setUpClass(cls):
            load_structure_to_class(cls, strname)
            crit = contacts.ContactsAlternative(contacts.CaCbxContact(), contacts.RingCenterContact())
            cls.struct.set_contact_map(crit)

        @testing(cydesc.CContactMap.__init__)
        @testing(cydesc.CContact.__init__)
        def test_init(self):
            cs = cydesc.CStructure(self.struct)
            ccm = cydesc.CContactMap(cs, self.struct.contact_map)

            n_contacts = sum(map(len, self.struct.contact_map.contacts.values()))

            self.assertEqual(ctypes.addressof(ccm.structure.contents), ctypes.addressof(cs))
            self.assertEqual(ccm.n_contacts, n_contacts)

            visited = {}

            for c in ccm.contacts[0:n_contacts]:
                self.assertEqual(c.val, self.struct.contact_map.get_contact_value(c.mer1, c.mer2))
                if (c.mer1, c.mer2) in visited:
                    self.fail("Contact duplicated")
                visited[(c.mer1, c.mer2)] = 1
                repr(c)

            repr(ccm)

        def test_del(self):
            cs = cydesc.CStructure(self.struct)
            ccm = cydesc.CContactMap(cs, self.struct.contact_map)

            del ccm
            # del cs

        def test_indices(self):
            def callback(mer1, mer2, val):
                # This is just for prettiness. It won't fail the test.
                self.assertEqual(self.struct.contact_map.get_contact_value(mer1, mer2), val)

                # This is ugly. Failing assertion inside a callback doesn't
                # prevent it from returning to C. Unittest in such a case is
                # unable to catch an exception and fail the test. We have to
                # pass the failure to C, check for it, and pass it back to
                # Python.
                if self.struct.contact_map.get_contact_value(mer1, mer2) != val:
                    return 1
                else:
                    return 0

            CBFUNC = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int)

            cs = cydesc.CStructure(self.struct)
            ccm = cydesc.CContactMap(cs, self.struct.contact_map)
            res = libcydesc_test.CContactMapTest_indices(ctypes.byref(ccm), CBFUNC(callback))
            self.assertEqual(res, 0)

        @classmethod
        def tearDownClass(cls):
            del_structure_from_class(cls)

    return CContactMapTest


def make_cdescriptortest(strname):
    """ Create and return a CDescriptorTest test case for a given structure. """

    @testing(cydesc.CDescriptor)
    @testing(cydesc.CElement)
    @test_name_append(strname)
    class CDescriptorTest(unittest.TestCase):

        @classmethod
        def setUpClass(cls):
            load_structure_to_class(cls, strname)
            crit = contacts.ContactsAlternative(contacts.CaCbxContact(), contacts.RingCenterContact())
            cls.struct.set_contact_map(crit)

            cls.descriptors = []

            for m in cls.struct:
                try:
                    cls.descriptors.append(structure.AbstractDescriptor.build(structure.Element.build(m)))
                except (TypeError, ValueError):
                    pass

        @testing(cydesc.CDescriptor.__init__)
        @testing(cydesc.CElement.__init__)
        def test_init(self):
            for desc in self.descriptors:
                cd = cydesc.CDescriptor(desc)
                ccm = cd.contact_map.contents

                self.assertEqual(ccm.n_contacts / 2, len(desc.contacts))

                contact_dict = defaultdict(list)

                for contact in desc.contacts:
                    contact_dict[contact.elements[0].central_monomer.ind].append(contact.elements[1].central_monomer.ind)
                    contact_dict[contact.elements[1].central_monomer.ind].append(contact.elements[0].central_monomer.ind)

                for el in cd.elements[0:cd.n_elements]:
                    self.assertIn(el.center, contact_dict.keys())

                self.assertEqual(cd.n_elements, len(contact_dict))

                cov = set([m.ind for el in cd.elements[0:cd.n_elements] for m in self.struct[el.start:el.end]])

                self.assertSetEqual(cov, set([m.ind for m in desc]))

                res = libcydesc_test.CDescriptorTest_element_map(ctypes.byref(cd))
                self.assertEqual(res, 0)

                repr(cd)

        @classmethod
        def tearDownClass(cls):
            del cls.descriptors
            del_structure_from_class(cls)

    return CDescriptorTest


@testing(cydesc.CInDelMeta)
@testing(cydesc.CInDelMeta.__call__)
@testing(cydesc.CInDelMeta.__new__)
@testing(cydesc.use_library)
class CInDelMetaTest(unittest.TestCase):
    class CTest(ctypes.Structure):
        __metaclass__ = cydesc.use_library(libcydesc_test)(cydesc.CInDelMeta)

        _fields_ = [('f1', ctypes.c_int),
                    ('f2', ctypes.c_int)]

    def run_test_lifecycle(self, cls):
        getcount = lambda: (libcydesc_test.get_init_CTest_calls(), libcydesc_test.get_del_CTest_calls())

        c1 = getcount()

        test = cls(1, 2)

        c2 = getcount()

        self.assertEqual(test.f1, 1)
        self.assertEqual(test.f2, 2)

        del test

        c3 = getcount()

        self.assertEqual((c2[0] - c1[0], c2[1] - c1[1]), (1, 0))
        self.assertEqual((c3[0] - c2[0], c3[1] - c2[1]), (0, 1))

    def test_lifecycle(self):
        self.run_test_lifecycle(self.CTest)

    def test_lifecycle_subclass(self):
        class CTest1(self.CTest):
            pass
        self.run_test_lifecycle(CTest1)

    def test_lifecycle_del(self):
        del_called = []

        class CTest1(self.CTest):
            def __del__(self):
                del_called.append(True)

        self.run_test_lifecycle(CTest1)
        self.assertEqual(del_called, [True])

    def test_lifecycle_del1(self):
        del_called = []

        class CTestA(self.CTest):
            def __del__(self):
                del_called.append(True)

        class CTest1(CTestA):
            pass

        self.run_test_lifecycle(CTest1)
        self.assertEqual(del_called, [True])


def tearDownModule():
    global struct_dict
    struct_dict = {}


def load_tests(loader, standard_tests, pattern):
    """ Add tests_legacy created by make_* functions for all structures. Return a complete TestSuite. """
    structures = ['1asz', '1gax', '1no5', '1pxq', '2dlc', '2lp2',
                  '3ftk', '3g88', '3lgb', '3m6x', '3npn', '3tk0', '3umy']

    if fast:
        structures = ['3umy']

    basic = unittest.TestSuite()

    for name in structures:
        basic.addTests(loader.loadTestsFromTestCase(make_cmertest(name)))
        basic.addTests(loader.loadTestsFromTestCase(make_cstructuretest(name)))
        basic.addTests(loader.loadTestsFromTestCase(make_ccontactmaptest(name)))
        basic.addTests(loader.loadTestsFromTestCase(make_cdescriptortest(name)))

    standard_tests.addTests(basic)

    return standard_tests

if __name__ == '__main__':
    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
