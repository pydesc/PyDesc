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
Unit tests_legacy for selection.py.

Usage:
    python selection_test.py [-v] [--fast]

    or

    python -m unittest [-v] selection_test
"""

import unittest
import random

import syntax_check
from syntax_check import notest, testing, test_name_append

import Bio.PDB

import warnings

import pydesc.selection as selection
import pydesc.structure as structure
import pydesc.numberconverter as numberconverter
import tests_legacy

import operator
#TO

syntax_check.module = selection

notest(selection.Selection)

fast = False

#~ fast = True

TestSyntax = syntax_check.module_syntax()

notest(selection.Selection.distinguish_chains)

# pylint: disable=C0111,R0912

random_generators = {}

data_dir = tests_legacy.__path__[0] + '/data/test_structures/'


def reg_ran(sel_cls):
    def wrapper(f):
        random_generators[sel_cls] = f
        return f

    return wrapper


@reg_ran(selection.Set)
def random_selection_set(struct):
    n = len(struct)
    inds = sorted(random.sample(range(n), random.randint(1, n)))

    mers = [struct[i] for i in inds]
    pdb_ids = [m.get_pdb_id() for m in mers]

    return (selection.Set(pdb_ids), mers)


@reg_ran(selection.Range)
def random_selection_range(struct):
    n = len(struct)
    (i_s, i_e) = sorted(random.sample(range(n), 2))

    mers = list(struct[i_s:i_e])
    (pi_s, pi_e) = [struct[i].get_pdb_id() for i in (i_s, i_e)]

    return (selection.Range(pi_s, pi_e), mers)


@reg_ran(selection.ChainSelection)
def random_selection_chain(struct):
    chain = random.choice(struct.chains)

    mers = list(chain)

    return (selection.ChainSelection(chain.chain_char), mers)


@reg_ran(selection.MonomerName)
def random_selection_monomername(struct):
    names = set(m.name for m in struct)
    name = random.choice(list(names))

    mers = [m for m in struct if m.name == name]

    return (selection.MonomerName(name), mers)


@reg_ran(selection.MonomerType)
def random_selection_monomertype(struct):
    types = set(m.__class__ for m in struct)
    t = random.choice(list(types))

    mers = [m for m in struct if isinstance(m, t)]

    return (selection.MonomerType(t), mers)


@reg_ran(selection.Everything)
def random_selection_everything(struct):
    mers = list(struct)

    return (selection.Everything(), mers)


@reg_ran(selection.Nothing)
def random_selection_nothing(struct):
    return (selection.Nothing(), [])


def random_selection(struct, sel_cls):
    return random_generators[sel_cls](struct)

def random_set(str):
    i1 = random.randint(0, len(str) - 1)
    return str[i1: random.randint(i1, len(str))]

def make_selectionbasictest(strname, sel_cls):
    """Create and return a SelectionTest testcase for a given structure."""
    @testing(selection.Selection)
    @testing(sel_cls)
    @test_name_append(sel_cls, strname)
    class SelectioniBasicTest(unittest.TestCase):
        name = strname

        @classmethod
        def setUpClass(cls):
            cls.pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(cls.name, data_dir + '%s.pdb' % cls.name)
            cls.converter = numberconverter.NumberConverter(cls.pdb_structure)
            cls.model = random.choice(cls.pdb_structure)
            with warnings.catch_warnings(record=True):
                cls.struct = structure.Structure(cls.model, cls.converter)

        @classmethod
        def tearDownClass(cls):
            del cls.struct
            del cls.model
            del cls.converter
            del cls.pdb_structure

        @testing(sel_cls.__init__)
        def test_init(self):
            for dummy in range(100):
                random_selection(self.struct, sel_cls)

        @testing(sel_cls.specify)
        def test_specify(self):
            for dummy in range(100):
                (sel, mers) = random_selection(self.struct, sel_cls)
                sel_set = set(sel.specify(self.struct))
                mer_set = set([m.get_pdb_id() for m in mers])
                self.assertSetEqual(sel_set, mer_set)

        @testing(sel_cls.create_structure)
        def test_create_structure(self):
            for dummy in range(100):
                (sel, mers) = random_selection(self.struct, sel_cls)
                #~ new_struct = set(sel.create_structure(self.struct))
                #TO
                new_struct = sel.create_structure(self.struct)
                #~ self.assertEqual(list(new_struct), mers)
                #TO created structure contains newly created mers, thus listed structures will never be equal
                #~ self.assertEqual(map(operator.attrgetter('pdb_residue'), new_struct), map(operator.attrgetter('pdb_residue'), mers))
                #TO
                self.assertEqual(map(operator.attrgetter('pdb_residue'), new_struct._monomers), map(operator.attrgetter('pdb_residue'), mers))

        if '__iter__' in sel_cls.__dict__:
            @testing(sel_cls.__iter__)
            def test_iter(self):
                (sel, mers) = random_selection(self.struct, sel_cls)
                mers = map(operator.attrgetter('ind'), mers)
                for i in sel:
                    self.assertTrue(self.struct.converter.get_ind(i) in mers)

        if 'create_new_structure' in sel_cls.__dict__:
            @testing(sel_cls.create_new_structure)
            def test_create_new_structure(self):
                (sel, mers) = random_selection(self.struct, sel_cls)
                nstr = sel.create_new_structure(self.struct)
                ids = list(sel)
                for i in nstr:
                    self.assertTrue(i.get_pdb_id() in ids)
                    self.assertNotEqual(i, self.struct[i.ind])

    return SelectioniBasicTest

def make_combinedselectiontest(strname, sel_cls):
    """Create and return a SelectionTest testcase for a given structure."""
    @testing(selection.CombinedSelection)
    @testing(sel_cls)
    @test_name_append(sel_cls, strname)
    class CombinedSelectionsTest(unittest.TestCase):
        name = strname

        @classmethod
        def setUpClass(cls):
            cls.basic_sel_cls = [random_selection_set, random_selection_range, random_selection_chain, random_selection_monomername, random_selection_monomertype, random_selection_everything, random_selection_nothing]
            cls.pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(cls.name, data_dir + '%s.pdb' % cls.name)
            cls.converter = numberconverter.NumberConverter(cls.pdb_structure)
            cls.model = random.choice(cls.pdb_structure)
            with warnings.catch_warnings(record=True):
                cls.struct = structure.Structure(cls.model, cls.converter)

        @classmethod
        def tearDownClass(cls):
            del cls.struct
            del cls.model
            del cls.converter
            del cls.pdb_structure
            del cls.basic_sel_cls

        @testing(sel_cls.__init__)
        @testing(sel_cls.__iter__)
        def test_init(self):
            sel = sel_cls([random.choice(self.basic_sel_cls)(self.struct)[0] for i in range(5)])
            self.assertEqual(5, len(list(iter(sel))))
            self.assertTrue(all(isinstance(i, selection.Selection) for i in sel))

        @testing(sel_cls.specify)
        def test_specify(self):
            mers1 = []
            sels = []
            for i in range(3):
                sel = random.choice(self.basic_sel_cls)(self.struct)
                if sel_cls == selection.SelectionsUnion:
                    mers1.extend(sel[1])
                elif sel_cls == selection.SelectionsIntersection:
                    mers1.append(set(sel[1]))
                else:       # relativecomplement
                    mers1.append(sel[1])
                sels.append(sel[0])
            if sel_cls == selection.SelectionsRelativeComplement:
                mers = [i for i in mers1[0] if i not in mers1[1]]
            elif sel_cls == selection.SelectionsIntersection:
                mers = set.intersection(*mers1)
            else:   # union
                mers = list(set(mers1))
            sel = sel_cls(sels)
            s_set = sel.specify(self.struct)
            spec_mers = sorted([self.struct.converter.get_ind(i) for i in s_set.ids])
            mers = sorted(map(operator.attrgetter('ind'), mers))
            if mers != spec_mers:
                import pdb; pdb.set_trace()
            self.assertEqual(spec_mers, mers, "%s: %s contains %s, while it should be %s" % (str(sel), ";".join(map(str, sel.selections)), str(spec_mers), str()))


    return CombinedSelectionsTest

def load_tests(loader, standard_tests, pattern):
    """ Add tests_legacy created by make_* functions for all structures. Returns a complete TestSuite. """
    structures = ['1asz', '1gax', '1no5', '1pxq', '2dlc', '2lp2',
                  '3ftk', '3g88', '3lgb', '3m6x', '3npn', '3tk0', '3umy']

    if fast:
        structures = ['3umy']

    basic = unittest.TestSuite()

    for name in structures:
        basic.addTests(loader.loadTestsFromTestCase(make_selectionbasictest(name, selection.Set)))
        basic.addTests(loader.loadTestsFromTestCase(make_selectionbasictest(name, selection.Range)))
        basic.addTests(loader.loadTestsFromTestCase(make_selectionbasictest(name, selection.ChainSelection)))
        basic.addTests(loader.loadTestsFromTestCase(make_selectionbasictest(name, selection.MonomerName)))
        basic.addTests(loader.loadTestsFromTestCase(make_selectionbasictest(name, selection.MonomerType)))
        basic.addTests(loader.loadTestsFromTestCase(make_selectionbasictest(name, selection.Everything)))
        basic.addTests(loader.loadTestsFromTestCase(make_selectionbasictest(name, selection.Nothing)))
        basic.addTests(loader.loadTestsFromTestCase(make_combinedselectiontest(name, selection.SelectionsUnion)))
        basic.addTests(loader.loadTestsFromTestCase(make_combinedselectiontest(name, selection.SelectionsIntersection)))
        basic.addTests(loader.loadTestsFromTestCase(make_combinedselectiontest(name, selection.SelectionsRelativeComplement)))

    standard_tests.addTests(basic)

    return standard_tests

if __name__ == '__main__':
    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
