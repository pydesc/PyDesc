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
Unit tests_legacy for config.py.

Usage:
    python config_test.py [-v]

    or

    python -m unittest [-v] config_test
"""
import syntax_check
from syntax_check import notest, testing, randomword, expando
import unittest

import random

import pydesc.config as config

syntax_check.module = config

TestSyntax = syntax_check.module_syntax()

notest(config.ConfigManager.show_config)


# pylint: disable=C0111

class BasicTests(unittest.TestCase):

    @testing(config.Param)
    @testing(config.Param.__init__)
    def test_param_init(self):
        name = randomword(8)
        val = expando()

        par = config.Param(name, val)

        self.assertEqual((name, val, None), (par.name, par.value, par.default_value))

        name = randomword(8)
        val = expando()
        def_val = expando()

        par = config.Param(name, val, def_val)

        self.assertEqual((name, val, def_val), (par.name, par.value, par.default_value))

    @testing(config.Branch)
    @testing(config.Branch.__init__)
    def test_branch_init(self):
        name = randomword(8)

        br = config.Branch(name)

        self.assertEqual(name, br.name)


@testing(config.Branch)
class BranchTests(unittest.TestCase):

    def setUp(self):
        self.branch = config.Branch('branch')

    @testing(config.Branch.new_branch)
    @testing(config.Branch.del_branch)
    def test_new_branch(self):
        self.branch.new_branch('branch1')
        self.assertTrue(hasattr(self.branch, 'branch1'))

        self.assertEqual(self.branch.branch1.name, 'branch.branch1')

        self.branch.del_branch('branch1')
        self.assertFalse(hasattr(self.branch, 'branch1'))

    def test_del_branch(self):
        self.branch.new_branch('branch1')
        self.assertTrue(hasattr(self.branch, 'branch1'))

        del self.branch.branch1
        self.assertFalse(hasattr(self.branch, 'branch1'))

    def test_namechange(self):
        self.branch.new_branch('branch1')

        self.assertEqual(self.branch.branch1.name, 'branch.branch1')

        self.branch.name = 'new_branch'
        self.assertEqual(self.branch.branch1.name, 'new_branch.branch1')

    @testing(config.Branch.set)
    def test_setget(self):
        val = expando()

        self.branch.set('par1', val)
        self.assertTrue(hasattr(self.branch, 'par1'))
        self.assertIs(getattr(self.branch, 'par1'), val)
        self.assertIs(self.branch.par1, val)

        name = randomword(8)
        val = expando()

        self.branch.set(name, val)

        self.assertTrue(hasattr(self.branch, name))
        self.assertIs(getattr(self.branch, name), val)

    def test_classprop(self):
        val = expando()

        self.branch.set('par1', val)

        branch = config.Branch('br1')
        self.assertFalse(hasattr(branch, 'par1'))
        self.assertNotIn('par1', dir(branch))

    @testing(config.Branch.delete)
    def test_delete(self):
        val = expando()

        self.branch.set('par1', val)
        self.assertIs(self.branch.par1, val)

        self.branch.delete('par1')
        self.assertRaisesRegexp(AttributeError, "\'par1\'", lambda: self.branch.par1)

    @testing(config.Branch.set_default)
    def test_set_default(self):
        val = expando()

        self.branch.set_default('par1', val)
        self.assertIs(self.branch.par1, val)

        self.branch.delete('par1')
        self.assertTrue(hasattr(self.branch, 'par1'))
        self.assertIs(self.branch.par1, val)

    def test_delete_default(self):
        val = expando()
        val1 = expando()

        self.branch.set('par1', val, val1)
        self.assertIs(self.branch.par1, val)

        self.branch.delete('par1')
        self.assertTrue(hasattr(self.branch, 'par1'))
        self.assertIs(self.branch.par1, val1)

    def test_delete_del(self):
        val = expando()
        val1 = expando()

        self.branch.set('par1', val, val1)
        self.assertIs(self.branch.par1, val)

        del self.branch.par1
        self.assertTrue(hasattr(self.branch, 'par1'))
        self.assertIs(self.branch.par1, val1)

    def test_parbranch_mix(self):
        self.branch.set('par1', 'val', 'default')

        self.branch.new_branch('par1')
        self.branch.delete('par1')

        self.assertEqual(self.branch.par1, 'val')

    @testing(config.Branch.get)
    def test_parbranch_mix1(self):
        self.branch.new_branch('branch1')
        self.assertEqual(self.branch.get('branch1'), self.branch.branch1)

    def test_getnested(self):
        self.branch.new_branch('branch1')

        self.branch.set('branch1.par', 'val')
        self.assertTrue(hasattr(self.branch.branch1, 'par'))
        self.assertEqual(self.branch.branch1.par, 'val')

        self.branch.branch1.set('par', 'val1')
        self.assertEqual(self.branch.get('branch1.par'), 'val1')

    @unittest.skipUnless(hasattr(config.Branch, '_list_dict'), "config.Branch doesn't have _list_dict.")
    def test_list_dict(self):
        try:
            self.branch._list_dict()  # pylint: disable=E1120, W0212
        except TypeError:
            self.fail('_list_dict requires an unnecessary argument.')

#@unittest.SkipTest


@testing(config.ConfigManager)
class ConfigManagerTests(unittest.TestCase):

    def setUp(self):
        def populate(br, depth=0):
            while random.randint(0, 9) < 8:
                if not depth or random.randint(0, 5) < 3:
                    br.set(randomword(8), randomword(8), randomword(8))
                else:
                    name = randomword(8)
                    br.new_branch(name)
                    nbr = getattr(br, name)
                    populate(nbr, depth - 1)

        config.ConfigManager.new_branch('root')
        populate(config.ConfigManager.root, 5)

    @testing(config.ConfigManager.new_branch)
    @testing(config.ConfigManager.del_branch)
    def test_new_branch(self):
        config.ConfigManager.new_branch('branch1')
        self.assertTrue(hasattr(config.ConfigManager, 'branch1'))

        self.assertEqual(config.ConfigManager.branch1.name, 'ConfigManager.branch1')

        config.ConfigManager.del_branch('branch1')
        self.assertFalse(hasattr(config.ConfigManager, 'branch1'))

    @testing(config.ConfigManager.save_config)
    def test_save(self):
        fname = '/tmp/config_test_save'

        config.ConfigManager.save_config(fname)

    @testing(config.ConfigManager.save_config)
    def test_branch_sig(self):
        fname = '/tmp/config_test_branch_sig'

        config.ConfigManager.root.set('par', 'branch_sig')
        config.ConfigManager.save_config(fname)

    @testing(config.ConfigManager.load_config)
    def test_load(self):
        fname = '/tmp/config_test_load'

        config.ConfigManager.save_config(fname)

        olddict = config.ConfigManager._get_dict()
        del config.ConfigManager.root
        config.ConfigManager.load_config(fname)
        newdict = config.ConfigManager._get_dict()

        self.assertEqual(olddict, newdict)

    def tearDown(self):
        try:
            del config.ConfigManager.root
        except:
            pass


if __name__ == '__main__':
    unittest.main()
