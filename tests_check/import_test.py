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
Unit tests for importing python modules.

Should be updated each time a new module is added.

Usage:
    python import_test.py [-v] [--skip-slow] [--fast]

    or

    python -m unittest [-v] import_test

created: 5.05.2014, Pawel Daniluk
"""

import unittest
import sys
import imp


# pylint: disable=C0111,R0201

class ImportTest(unittest.TestCase):
    modules = ['pydesc', 'pydesc.config', 'pydesc.contacts',
               'pydesc.cydesc', 'pydesc.cydesc.fitdesc', 'pydesc.cydesc.overfit',
               'pydesc.compat', 'pydesc.compat.dsc',
               'pydesc.dbhandler', 'pydesc.geometry', 'pydesc.monomer', 'pydesc.numberconverter', 'pydesc.selection',
               'pydesc.structure', 'pydesc.warnexcept',
               'tests', 'tests.compat_dsc_test', 'tests.config_test', 'tests.contacts_test', 'tests.cydesc_test',
               'tests.dbhandler_test', 'tests.fitdesc_test', 'tests.geometry_test',  # 'tests.import_test',
               'tests.monomer_test', 'tests.numberconverter_test', 'tests.overfit_test', 'tests.selection_test',
               'tests.structure_test', 'tests.syntax_check', 'tests.warnexcept_test']

    def try_import(self, name):
        if '.' in name:
            parent_mod, mod_name = name.rsplit('.', 1)
            self.try_import(parent_mod)
            mod_path = sys.modules[parent_mod].__path__
        else:
            mod_name = name
            mod_path = None
            parent_mod = None

        # code taken from http://docs.python.org/2/library/imp.html
        try:
            return sys.modules[name]
        except KeyError:
            pass

        # If any of the following calls raises an exception,
        # there's a problem we can't handle -- let the caller handle it.

        file_obj, pathname, description = imp.find_module(mod_name, mod_path)

        try:
            res = imp.load_module(name, file_obj, pathname, description)
            if parent_mod is not None:
                setattr(sys.modules[parent_mod], mod_name, res)

            if hasattr(res, 'load_tests'):
                res.load_tests(unittest.defaultTestLoader, unittest.TestSuite(), None)

        except:
            self.failures.append(name)

        finally:
            # Since we may exit via an exception, close fp explicitly.
            if file_obj:
                file_obj.close()

    def del_all(self):
        for n in sys.modules.keys():
            if n.startswith('pydesc.'):
                del sys.modules[n]

    def test_imports(self):
        self.failures = []
        for name in self.modules:
            self.try_import(name)
            self.del_all()

        if len(self.failures) > 0:
            self.fail("The following modules failed to load: %s" % str(self.failures))
