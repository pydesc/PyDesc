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
Unit tests_legacy for cydesc/fitdesc.py.

Usage:
    python fitdesc_test.py [-v] [--skip-slow] [--fast]

    or

    python -m unittest [-v] fitdesc_test

created: 22.04.2014, Pawel Daniluk
"""

import unittest
import random
import itertools
import functools
import operator
import multiprocessing


import tests_legacy.syntax_check as syntax_check
from tests_legacy.syntax_check import notest, testing

import Bio.PDB

import warnings

import pydesc.structure as structure
import pydesc.numberconverter as numberconverter
import pydesc.geometry as geometry
#import pydesc.cydesc as cydesc
import pydesc.cydesc.fitdesc as fitdesc
import tests_legacy


syntax_check.module = fitdesc

notest(fitdesc.FitDesc.motif)
notest(fitdesc.FitDesc.structure)

data_dir = tests_legacy.__path__[0] + '/data/test_structures/'

skip_slow = False
fast = False

TestSyntax = syntax_check.module_syntax()

# pylint: disable=C0111, R0201, R0912


class TimeoutError(Exception):
    pass


def timeout(seconds=5, error_message="Timeout"):
    def decorator(func):
        def wrapper(*args, **kwargs):
            process = multiprocessing.Process(None, func, None, args, kwargs)
            process.start()
            process.join(seconds)
            if process.is_alive():
                process.terminate()
                raise TimeoutError(error_message)

        return functools.wraps(func)(wrapper)
    return decorator


@testing(fitdesc.FitDesc)
class FitDescTest(unittest.TestCase):
    struct_names = ['1asz', '1gax', '1no5', '1pxq', '2dlc', '2lp2',
                    '3ftk', '3g88', '3lgb', '3m6x', '3npn', '3tk0', '3umy']

    struct_names = ['1asz', '1no5', '1pxq', '3ftk', '3m6x', '3tk0', '3umy']

    longMessage = True

    RMSD_threshold = 5

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
        del cls.points1
        del cls.points2
        del cls.structs

    def assertIdent(self, motif, res, msg):
        almost_equal = functools.partial(self.assertAlmostEqual, msg=msg)

        rmsd, pair_al, mat = res
        almost_equal(rmsd, 0, places=3)

        one = geometry.TRTMatrix()
        map(almost_equal, itertools.chain(*mat.rotation_matrix.tolist()), itertools.chain(*one.rotation_matrix.tolist()))
        map(almost_equal, mat.prerotational_translation_vector.tolist(), one.prerotational_translation_vector.tolist())
        map(almost_equal, mat.translation_vector.tolist(), one.translation_vector.tolist())

        self.assertListEqual(*map(list, zip(*list(pair_al))), msg=msg)
        self.assertListEqual(list(motif), list(zip(*list(pair_al))[0]), msg=msg)

    def random_segment(self, par_struct=None):
        while True:
            if par_struct is None:
                struct = random.choice(self.structs)
            else:
                struct = par_struct

            seg_len = random.randint(5, 20)

            try:
                start = random.choice(list(struct[0:-seg_len]))
                sub = struct[start.ind:start.ind + seg_len]

                return (sub, struct)
            except IndexError:
                pass

    def random_multiseg(self, par_struct=None):
        while True:
            if par_struct is None:
                struct = random.choice(self.structs)
            else:
                struct = par_struct

            n_seg = random.randint(2, 4)

            segs = []

            for dummy in range(n_seg):
                try:
                    segs.append(self.random_segment(struct)[0])
                except:
                    raise

            if len(segs) == n_seg:
                try:
                    sub = reduce(operator.add, segs)
                except:
                    raise

                return (sub, struct)

    def identity_template(self, random_struct):
        sub, struct = random_struct()

        fitdesc_obj = fitdesc.FitDesc(sub, sub)

        msg = "%s %s" % (str(struct), str(sub))

        # A super ugly hack to fail whenever fitdesc takes too long.
        # If fitdesc run in a separate process (which can be terminated) ends
        # in time, it is run second time to collect the result.
        try:
            timeout(30)(fitdesc_obj.fitdesc)(self.RMSD_threshold)
        except TimeoutError:
            self.fail("Timeout: %s" % msg)

        res = fitdesc_obj.fitdesc(self.RMSD_threshold)

        self.assertGreaterEqual(len(res), 1, msg=msg)
        self.assertIdent(sub, res[0], msg)

    def find_in_struct_template(self, random_struct):
        sub, struct = random_struct()

        msg = "%s %s" % (str(struct), str(sub))

        fitdesc_obj = fitdesc.FitDesc(sub, struct)

        # A super ugly hack to fail whenever fitdesc takes too long.
        # If fitdesc run in a separate process (which can be terminated) ends
        # in time, it is run second time to collect the result.
        try:
            res = timeout(30)(fitdesc_obj.fitdesc)(self.RMSD_threshold)
        except TimeoutError:
            self.fail("Timeout: %s" % msg)

        res = fitdesc_obj.fitdesc(self.RMSD_threshold)

        ok = False

        for r in res:
            try:
                self.assertIdent(sub, r, msg)

                ok = True
            except AssertionError:
                pass

        self.assertTrue(ok, msg=msg + "Motif not found")

    @testing(fitdesc.FitDesc.__init__)
    def test_init(self):
        fitdesc.FitDesc()

    @testing(fitdesc.FitDesc.fitdesc)
    @testing(fitdesc.t_fitdesc_result)
    @testing(fitdesc.t_fitdesc_result.unpack)
    def test_identity(self):
        for dummy in range(100):
            self.identity_template(self.random_segment)

    def test_find_in_struct(self):
        for dummy in range(100):
            self.find_in_struct_template(self.random_segment)

    def test_identity_multiseg(self):
        for dummy in range(50):
            self.identity_template(self.random_multiseg)

    def test_find_in_struct_multiseg(self):
        for dummy in range(50):
            self.find_in_struct_template(self.random_multiseg)

    @testing(fitdesc.fitdesc)
    def test_fitdesc_function(self):
        sub, struct = self.random_multiseg()

        msg = "%s %s" % (str(struct), str(sub))

        fitdesc_obj = fitdesc.FitDesc(sub, struct)

        # A super ugly hack to fail whenever fitdesc takes too long.
        # If fitdesc run in a separate process (which can be terminated) ends
        # in time, it is run second time to collect the result.
        try:
            res = timeout(30)(fitdesc_obj.fitdesc)(self.RMSD_threshold)
        except TimeoutError:
            self.fail("Timeout: %s" % msg)

        res = fitdesc_obj.fitdesc(self.RMSD_threshold)

        res1 = fitdesc.fitdesc(sub, struct, self.RMSD_threshold)

        self.assertEqual(len(res), len(res1))

if __name__ == '__main__':
    if syntax_check.rip_argv('--skip-slow'):
        skip_slow = True

    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
