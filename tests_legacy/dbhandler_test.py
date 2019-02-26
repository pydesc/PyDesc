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
Unit tests_legacy for dbhandler.py.

Usage:
    python dbhandler_test.py [-v] [--fast]

    or

    python -m unittest [-v] dbhandler_test

Tymoteusz
"""
import syntax_check
from syntax_check import notest, test
import unittest

import os
import itertools
from collections import defaultdict

import pydesc.dbhandler as dbhandler
from pydesc.config import ConfigManager
import tests_legacy

syntax_check.module = dbhandler

TestSyntax = syntax_check.module_syntax()

notest(dbhandler.add_db_dir)
notest(dbhandler.validate_id)
notest(dbhandler.InvalidID)
notest(dbhandler.DBHandler)

data_dir = tests_legacy.__path__[0] + '/data/dbtest/'

fast = False

fast = True

# pylint: disable=C0111


class HandlerTestMixin(object):

    def setUp(self):
        self.handler = self.cls()

        test(self.cls)
        test(self.cls.__init__)
        try:
            notest(self.cls.is_file_valid)
        except AttributeError:
            pass

    def test_is_id_valid(self):
        test(self.cls.is_id_valid)

        for i in self.valid_ids:
            self.assertTrue(self.handler.is_id_valid(i), "Id %s should be valid." % i)

        for i in self.invalid_ids:
            self.assertFalse(self.handler.is_id_valid(i), "Id %s should be invalid." % i)

    def test_get_file_location(self):
        test(self.cls.get_file_location)

        for i in self.valid_ids:
            loc = self.handler.get_file_location(i)
            self.assertTrue(loc.startswith(ConfigManager.dbhandler.cachedir))

        for i in self.invalid_ids + self.nonexistent_ids:
            self.assertRaises(dbhandler.InvalidID, self.handler.get_file_location, i)

    def test_get_file_url(self):
        test(self.cls.get_file_url)

        for i in self.valid_ids:
            loc = self.handler.get_file_url(i)
            self.assertTrue(loc.startswith(ConfigManager.dbhandler.cachedir))

        for i in self.invalid_ids + self.nonexistent_ids:
            self.assertRaises(dbhandler.InvalidID, self.handler.get_file_url, i)

    def test_download_file(self):
        test(self.cls.download_file)

        for (i, ref) in zip(self.valid_ids, self.reference_files):
            self.handler.download_file(i)
            dwnld = self.handler.get_file_location(i)
            with open(dwnld, 'r') as f:
                buf1 = f.read()

            with open(ref, 'r') as f:
                buf2 = f.read()

            self.assertEqual(buf1, buf2, "Files %s and %s differ." % (dwnld, ref))

            if fast:
                break

        for i in self.invalid_ids:
            self.assertRaises(dbhandler.InvalidID, self.handler.download_file, i)

    def test_get_file(self):
        test(self.cls.get_file)

        for (i, ref) in zip(self.valid_ids, self.reference_files):

            dwnld = self.handler.get_file_location(i)
            with open(ref, 'r') as f:
                buf2 = f.read()

            buf1 = self.handler.get_file(i, 1).read()
            self.assertEqual(buf1, buf2, "Files %s and %s differ." % (dwnld, ref))

            buf1 = self.handler.get_file(i, 2).read()
            self.assertEqual(buf1, buf2, "Files %s and %s differ." % (dwnld, ref))

            os.remove(dwnld)

            self.assertRaises(Exception, self.handler.get_file, i, 2)

            buf1 = self.handler.get_file(i, 0).read()
            self.assertEqual(buf1, buf2, "Files %s and %s differ." % (dwnld, ref))

            os.remove(dwnld)
            with open(dwnld, 'w') as f:
                f.write('crap')

            buf1 = self.handler.get_file(i, 0).read()
            self.assertEqual(buf1, 'crap', "Files %s and %s differ." % (dwnld, ref))

            buf1 = self.handler.get_file(i, 2).read()
            self.assertEqual(buf1, 'crap', "Files %s and %s differ." % (dwnld, ref))

            buf1 = self.handler.get_file(i, 1).read()
            self.assertEqual(buf1, buf2, "Files %s and %s differ." % (dwnld, ref))

            if fast:
                break


class ScopTest(HandlerTestMixin, unittest.TestCase):
    cls = dbhandler.SCOPHandler

    valid_ids = ['d1nkla_', 'd3pi2a_', 'd1abwa2']
    nonexistent_ids = ['d114za_']
    invalid_ids = ['x1nkla2', 'dzxxx', '1asz']

    reference_files = [data_dir + 'SCOP/' + n for n in valid_ids]


class PDBTest(HandlerTestMixin, unittest.TestCase):
    cls = dbhandler.PDBHandler

    valid_ids = ['1nkl', '3pi2', '1abw']
    nonexistent_ids = ['114z']
    invalid_ids = ['zxxx']

    reference_files = [data_dir + 'PDB/' + n for n in valid_ids]


class CATHTest(HandlerTestMixin, unittest.TestCase):
    cls = dbhandler.CATHHandler

    valid_ids = ['1nklA00', '1abwA', '1nkl']
    nonexistent_ids = ['114zA00']
    invalid_ids = ['zxxx', '1abwA1']

    reference_files = [data_dir + 'CATH/' + n for n in valid_ids]


class BioUnitTest(HandlerTestMixin, unittest.TestCase):
    cls = dbhandler.BioUnitHandler

    valid_ids = ['3pi2/1', '3pi2/2', '3pi2/3']  # Each structure has to be given with all units in ascending order.
    nonexistent_ids = ['3pi2/4']
    invalid_ids = ['zxxx', '1abw/', '3pi2/x']

    reference_files = [data_dir + 'PDBUnit/' + n for n in valid_ids]

    def test_get_file_multi(self):
        valid_dict = defaultdict(list)

        for (i, f) in zip(self.valid_ids, self.reference_files):
            (name, num) = i.split('/')
            valid_dict[name].append((num, f))

        for n in valid_dict:
            dwnld_list = []
            buf2_list = []

            for (i, ref) in valid_dict[n]:
                dwnld_list.append(self.handler.get_file_location('%s/%s' % (n, i)))
                with open(ref, 'r') as f:
                    buf2_list.append(f.read())

            buf1_list = self.handler.get_file(i, 1).read()
            self.assertEqual(buf1_list, buf2_list, "Files differ (%s)." % (n))

            buf1_list = self.handler.get_file(i, 2).read()
            self.assertEqual(buf1_list, buf2_list, "Files differ (%s)." % (n))

            map(os.remove, dwnld_list)

            self.assertRaises(Exception, self.handler.get_file, i, 2)

            buf1_list = self.handler.get_file(i, 0).read()
            self.assertEqual(buf1_list, buf2_list, "Files differ (%s)." % (n))

            map(os.remove, dwnld_list)
            for dwnld in dwnld_list:
                with open(dwnld, 'w') as f:
                    f.write('crap')

            buf1_list = self.handler.get_file(i, 0).read()
            self.assertEqual(buf1_list, ['crap'] * len(valid_dict[n]), "Files differ (%s)." % (n))

            buf1_list = self.handler.get_file(i, 2).read()
            self.assertEqual(buf1_list, ['crap'] * len(valid_dict[n]), "Files differ (%s)." % (n))

            buf1_list = self.handler.get_file(i, 1).read()
            self.assertEqual(buf1_list, buf2_list, "Files differ (%s)." % (n))

            if fast:
                break


class MetaTest(HandlerTestMixin, unittest.TestCase):
    cls = dbhandler.MetaHandler

    prefdict = {ScopTest: 'scop://', PDBTest: 'pdb://', CATHTest: 'cath://', BioUnitTest: 'unit://'}

    def __init__(self, *args, **kwds):
        addpref = lambda attr: list(itertools.chain(*[[self.prefdict[c] + x for x in getattr(c, attr)] for c in self.prefdict]))

        self.valid_ids = addpref('valid_ids')
        self.nonexistent_ids = addpref('nonexistent_ids')
        self.invalid_ids = addpref('invalid_ids')

        self.reference_files = list(itertools.chain(*[c.reference_files for c in self.prefdict]))

        super(MetaTest, self).__init__(*args, **kwds)

if __name__ == '__main__':
    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
