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
Unit tests_legacy for cydesc/overfit.py.

Usage:
    python overfit_test.py [-v] [--skip-slow] [--fast]

    or

    python -m unittest [-v] overfit_test

created: 19.04.2014, Pawel Daniluk
"""

import unittest
import random
import itertools
import functools
import operator
import math
import ctypes

import syntax_check
from syntax_check import notest, testing

import Bio.PDB

import warnings

import pydesc.structure as structure
import pydesc.numberconverter as numberconverter
import pydesc.geometry as geometry
import pydesc.cydesc as cydesc
import pydesc.cydesc.overfit as overfit
import tests_legacy

syntax_check.module = overfit

notest(overfit.c_float_4)
notest(overfit.Overfit.token_context)

data_dir = tests_legacy.__path__[0] + '/data/test_structures/'

skip_slow = False
fast = False

TestSyntax = syntax_check.module_syntax()

# pylint: disable=C0111, R0201, R0912


def assert_sums_equal(self, sums1, sums2):
    for field_name, field_type in overfit.t_overfit_sums._fields_:
        if field_type is overfit.c_float_4:
            for i in range(4):
                v1, v2 = [getattr(s, field_name)[i] for s in [sums1, sums2]]
                self.assertAlmostEqual(v1, v2, places=3, msg="%s[%d]: %s != %s" % (field_name, i, str(v1), str(v2)))

        else:
            v1, v2 = [getattr(s, field_name) for s in [sums1, sums2]]
            self.assertAlmostEqual(v1, v2, places=3, msg="%s: %s != %s" % (field_name, str(v1), str(v2)))


@testing(overfit.Overfit)
class OverfitTest(unittest.TestCase):
    struct_names = ['1asz', '1gax', '1no5', '1pxq', '2dlc', '2lp2',
                    '3ftk', '3g88', '3lgb', '3m6x', '3npn', '3tk0', '3umy']

    @classmethod
    def setUpClass(cls):
        cls.structs = []

        if not skip_slow:
            for strname in cls.struct_names:
                pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(strname, data_dir + '%s.pdb' % strname)
                converter = numberconverter.NumberConverter(pdb_structure)
                model = random.choice(pdb_structure)
                with warnings.catch_warnings(record=True):
                    try:
                        cls.structs.append(structure.Structure(model, converter))
                    except:
                        print "Failed to load %s" % strname

            if len(cls.structs) == 0:
                print "All structures failed to load!!!"

        cls.points1 = [syntax_check.randomcoord() for dummy in range(100)]
        cls.points2 = [syntax_check.randomcoord() for dummy in range(100)]

    @classmethod
    def tearDownClass(cls):
        del cls.structs
        del cls.points1
        del cls.points2

    @testing(overfit.Overfit.__init__)
    def test_init(self):
        overfit.Overfit()

    @testing(overfit.Overfit.__init__)
    def test_zero(self):
        overfit_obj = overfit.Overfit()

        sums = overfit_obj.get_sums()

        assert_sums_equal(self, sums, overfit.t_overfit_sums())

    @testing(overfit.Overfit.add_sums)
    @testing(overfit.Overfit.add_point)
    @testing(overfit.Overfit.get_sums)
    @testing(overfit.t_overfit_sums)
    @testing(overfit.t_overfit_sums.__add__)
    def test_get_add_sums(self):

        sum_list = []


        for dummy in range(100):
            overfit_obj = overfit.Overfit()
            overfit_obj1 = overfit.Overfit()
            part_list = []
            for point1, point2 in zip(self.points1, self.points2):
                if random.choice([True, False]):
                    sums = overfit_obj.get_sums()
                    overfit_obj = overfit.Overfit()
                    overfit_obj.add_sums(sums)

                    part_list.append(overfit_obj1.get_sums())
                    overfit_obj1 = overfit.Overfit()

                overfit_obj.add_point(point1, point2)
                overfit_obj1.add_point(point1, point2)

            part_list.append(overfit_obj1.get_sums())
            sum_list.append(overfit_obj.get_sums())
            assert_sums_equal(self, sum_list[-1], functools.reduce(operator.add, part_list))

        for sums1, sums2 in zip(sum_list, sum_list[1:]):
            assert_sums_equal(self, sums1, sums2)


    @testing(overfit.Overfit.overfit)
    def test_identity(self):
        for dummy in range(10):
            point_list = [self.points1[n] for n in random.sample(range(100), random.randint(3, 50))]

            overfit_obj = overfit.Overfit()
            overfit_obj.add(point_list, point_list)

            rmsd, mat = overfit_obj.overfit()

            self.assertAlmostEqual(rmsd, 0, places=3)

            one = geometry.TRTMatrix()
            map(self.assertAlmostEqual, itertools.chain(*mat.rotation_matrix.tolist()), itertools.chain(*one.rotation_matrix.tolist()))
            map(self.assertAlmostEqual, mat.prerotational_translation_vector.tolist(), one.prerotational_translation_vector.tolist())
            map(self.assertAlmostEqual, mat.translation_vector.tolist(), one.translation_vector.tolist())

    @testing(overfit.Overfit.overfit)
    @testing(overfit.t_transrot)
    @testing(overfit.t_transrot.to_trtmatrix)
    def test_geom_center(self):
        for dummy in range(10):
            n_points = random.randint(3, 50)
            point_list1 = [self.points1[n] for n in random.sample(range(100), n_points)]
            point_list2 = [self.points2[n] for n in random.sample(range(100), n_points)]

            overfit_obj = overfit.Overfit()
            overfit_obj.add(point_list1, point_list2)

            dummy_rmsd, mat = overfit_obj.overfit()

            center1, center2 = [[sum(map(operator.itemgetter(pos), plist)) / n_points for pos in range(3)] for plist in (point_list1, point_list2)]

            map(functools.partial(self.assertAlmostEqual, places=4), center1, mat.transform(*center2))

    @testing(overfit.Overfit.overfit)
    @testing(overfit.Overfit.add)
    @testing(overfit.t_transrot.to_trtmatrix)
    @testing(overfit.overfit)
    def test_rmsd(self):
        for dummy in range(10):
            n_points = random.randint(3, 50)
            point_list1 = [self.points1[n] for n in random.sample(range(100), n_points)]
            point_list2 = [self.points2[n] for n in random.sample(range(100), n_points)]

            overfit_obj = overfit.Overfit()
            overfit_obj.add(point_list1, point_list2)

            rmsd, mat = overfit_obj.overfit()

            rmsd1, mat1 = overfit.overfit(point_list1, point_list2)

            self.assertAlmostEqual(rmsd, rmsd1)
            map(self.assertAlmostEqual, itertools.chain(*mat.rotation_matrix.tolist()), itertools.chain(*mat1.rotation_matrix.tolist()))
            map(self.assertAlmostEqual, mat.prerotational_translation_vector, mat1.prerotational_translation_vector)
            map(self.assertAlmostEqual, mat.translation_vector, mat1.translation_vector)

            sqr = lambda x: x ** 2
            dist2 = lambda point1, point2: sum([sqr(c1 - c2) for c1, c2 in zip(point1, point2)])

            myrmsd = math.sqrt(sum(map(dist2, point_list1, [mat.transform(*p) for p in point_list2])) / n_points)

            self.assertAlmostEqual(rmsd, myrmsd, places=3)

    @unittest.skipIf(skip_slow, "Too slow.")
    @testing(overfit.Overfit.add_mer)
    def test_add_mer(self):
        for dummy in range(10):
            overfit_obj1 = overfit.Overfit()
            overfit_obj2 = overfit.Overfit()

            for dummy in range(100):
                struct = random.choice(self.structs)
                mer = random.choice(list(struct))

                overfit_obj1.add_mer(mer, mer)

                for point in mer.representation:
                    overfit_obj2.add_point(point, point)

                sums1 = overfit_obj1.get_sums()
                sums2 = overfit_obj2.get_sums()

                assert_sums_equal(self, sums1, sums2)

    @unittest.skipIf(skip_slow, "Too slow.")
    @testing(overfit.Overfit.add_structure)
    def test_add_struct(self):
        for dummy in range(10):
            overfit_obj1 = overfit.Overfit()
            overfit_obj2 = overfit.Overfit()

            for dummy in range(100):
                struct1, struct2 = random.sample(self.structs, 2)

                seg_len = random.randint(5, 20)

                try:
                    start1 = random.choice(list(struct1[0:-seg_len]))
                    start2 = random.choice(list(struct2[0:-seg_len]))

                    sub1 = struct1[start1.ind:start1.ind + seg_len]
                    sub2 = struct2[start2.ind:start2.ind + seg_len]
                except:
                    continue

                will_fail = False

                for mer1, mer2 in zip(*map(list, [sub1, sub2])):
                    if len(mer1.indicators) == len(mer2.indicators):
                        overfit_obj2.add_mer(mer1, mer2)
                    else:
                        self.assertRaises(TypeError, overfit_obj2.add_mer, mer1, mer2)
                        will_fail = True

                if will_fail:
                    self.assertRaises(TypeError, overfit_obj1.add_structure, sub1, sub2)
                    break

                overfit_obj1.add_structure(sub1, sub2)

                sums1 = overfit_obj1.get_sums()
                sums2 = overfit_obj2.get_sums()

                assert_sums_equal(self, sums1, sums2)


if __name__ == '__main__':
    if syntax_check.rip_argv('--skip-slow'):
        skip_slow = True

    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
