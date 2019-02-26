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
Unit tests_legacy for warnexcept.py.

Usage:
    python warnexcept_test.py [-v]

    or

    python -m unittest [-v] warnexcept_test
"""
import syntax_check
from syntax_check import notest, testing
import warnings
import unittest
import random

import pydesc.warnexcept as warnexcept
from pydesc.config import ConfigManager

ConfigManager.warnings_and_exceptions.class_filters.set("CopyDownload", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("IncompleteChainableParticle", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("Info", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("DeprecationWarning", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("LocalCopyAccess", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("MonomerCreationFailed", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("NoConfiguration", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("UnknownParticleName", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("UserWarning", "always")
ConfigManager.warnings_and_exceptions.class_filters.set("WrongMonomerType", "always")

syntax_check.module = warnexcept

TestSyntax = syntax_check.module_syntax()

notest(warnexcept.DiscontinuityError)
notest(warnexcept.IncompleteParticle)
notest(warnexcept.WrongAtomDistances)
notest(warnexcept.WrongMonomerType)
notest(warnexcept.CannotCalculateContact)


@testing(warnexcept.WarnManager)
class WarnManagerTest(unittest.TestCase):


    @testing(warnexcept.WarnManager.__init__)
    @classmethod
    def setUpClass(cls):
        cls.wm = warnexcept.WarnManager('obj_1')

    @testing(warnexcept.WarnManager.__enter__)
    @testing(warnexcept.WarnManager.__exit__)
    def test_init(self):
        self.assertEquals(self.wm.obj, 'obj_1')
        with self.wm as w:
            2 + 2

    @testing(warnexcept.WarnManager.__call__)
    @testing(warnexcept.WarnManager.__enter__)
    @testing(warnexcept.WarnManager.__exit__)
    def test_call(self):
        wm2 = warnexcept.WarnManager('obj_2')
        for i in ['ctx%i' % j for j in range(10)]:
            with wm2(i) as man:
                self.assertEqual(wm2.last_context, i)
                self.assertEqual(man.last_context, i)
                self.assertEqual(man, wm2)

    @testing(warnexcept.WarnManager.__exit__)
    def test_raising_errors(self):
        def helping_func():
            with self.wm as man:
                print archimidiwirilopotoczereppapinczakowska
        self.assertRaises(NameError, lambda: helping_func())

    @testing(warnexcept.WarnManager.__exit__)
    @testing(warnexcept.WarnManager.__call__)
    def test_not_raising_warnings(self):
        with self.wm('no_warns'):
            raise Warning('lala')
        res = self.wm.exceptions['no_warns'][-1]
        self.assertTrue(issubclass(res[1], Warning))
        self.assertEqual(str(res[0]), 'lala')

    @testing(warnexcept.WarnManager.raise_all)
    def test_delayed_raising_warns(self):
        #~ hlp = {}
        with self.wm('warns'):
            for i, wr in enumerate([i for i in Warning.__subclasses__()]):
                warnings.simplefilter('always')
                warnings.warn('wr%i' % i, wr)
                #~ hlp[wr] = i
        with warnings.catch_warnings(record=True) as ctxwar:
            self.wm.raise_all('warns')
        self.assertEqual(len(ctxwar), len([i for i in Warning.__subclasses__()]))

    @testing(warnexcept.WarnManager.raise_all)
    @testing(warnexcept.WarnManager.__exit__)
    @testing(warnexcept.WarnManager.__init__)
    def test_multicontext_warns_raising(self):          # TODO: sprawdzanie typu warna
        lns = {}
        for i in range(10):
            with self.wm('ctx%i' % i):
                no = random.choice([1,2,3])
                lns[i] = no
                for j in range(no):
                    wr = random.choice([i for i in Warning.__subclasses__()])
                    warnings.simplefilter('always')
                    warnings.warn('bla', wr)
        for i in range(10):
            with warnings.catch_warnings(record=True) as ctxwar:
                self.wm.raise_all('ctx%i' % i, 'always')
            self.assertEqual(len(ctxwar), lns[i])



@testing(warnexcept.warn)
class WarnFxTest(unittest.TestCase):

    @testing(warnexcept.CopyDownload)
    @testing(warnexcept.Info)
    @testing(warnexcept.LocalCopyAccess)
#    @testing(warnexcept.MonomerCreationFailed)
    @testing(warnexcept.NoConfiguration)
    @testing(warnexcept.UnknownParticleName)
    @testing(warnexcept.IncompleteChainableParticle)
    @testing(warnexcept.PyDescWarn)
    @testing(warnexcept.PyDescWarn.__init__)
    def test_warn(self):
        for wrn in (warnexcept.CopyDownload, warnexcept.CopyDownload, warnexcept.Info, warnexcept.LocalCopyAccess, warnexcept.NoConfiguration, warnexcept.UnknownParticleName, warnexcept.IncompleteChainableParticle):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                warnexcept.warn('text', wrn)
        self.assertEqual(len(w), 1)
        self.assertTrue(issubclass(w[-1].category, warnexcept.IncompleteChainableParticle))
        self.assertEqual(str(w[-1].message), wrn.__name__ + ": text")


# pylint: disable=C0111

if __name__ == '__main__':
    unittest.main()
