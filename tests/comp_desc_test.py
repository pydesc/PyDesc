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

import unittest
import operator

import syntax_check
from syntax_check import notest, testing

import Bio.PDB

import pydesc.cydesc.compdesc as cpd
import pydesc.structure as stc

from pydesc.config import ConfigManager


ConfigManager.warnings_and_exceptions.set("quiet", True)

syntax_check.module = cpd

TestSyntax = syntax_check.module_syntax()

def create_basic_tests(pair, fast=False):

    @testing(cpd.CompDesc)
    @testing(cpd.compdesc)
    class BasicCompdescTests(unittest.TestCase):

        s1n = pair[0]
        s2n = pair[1]

        @classmethod
        def setUpClass(cls):
            cls.ldr = stc.StructureLoader()

            cls.s1 = cls.ldr.load_structure(cls.s1n)[0]
            cls.s2 = cls.ldr.load_structure(cls.s2n)[0]

            cls.s1.set_contact_map()
            cls.s2.set_contact_map()

            end = 10 if fast else None

            s1ds = stc.AbstractDescriptor.create_descriptors(cls.s1)
            s2ds = stc.AbstractDescriptor.create_descriptors(cls.s2)

            cls.rdic1 = {}
            cls.rdic2 = {}
            for d1, d2 in zip(s1ds, s2ds):
                try: cls.rdic1[d1.central_element.central_monomer.ind] = d1
                except AttributeError: pass
                try: cls.rdic2[d2.central_element.central_monomer.ind] = d2
                except AttributeError: pass

            cls.res = [(cpd.compdesc(d1, d2), d1, d2) for d1 in s1ds[:end] for d2 in s2ds[:end] if None not in (d1, d2)]

            if len(cls.res) == 0:
                raise Warning("No results from compdesc.")

        @classmethod
        def tearDownClass(cls):
            del cls.s1
            del cls.s2
            del cls.rdic1
            del cls.rdic1
            del cls.res

        @testing(cpd.compdesc)
        @testing(cpd.CompDesc)
        def test_symmetry(self):
            for r, d1, d2 in self.res:
                rs = cpd.compdesc(d2, d1)
                self.assertEqual(len(r), len(rs))
                for r1, r2 in zip(r, rs):
                    self.assertAlmostEqual(r1[0], r2[0])
                    self.assertEqual(r1[1], r2[1])

        def test_overall_rmsd(self):

            th = ConfigManager.compdesc.th_overall_rmsd
            rmsds = [sr[0] for r in self.res for sr in r[0] if len(sr) != 0]

            for rmsd in rmsds:
                self.assertLessEqual(rmsd, th)

        def test_aa_cover(self):
            for r, d1, d2 in self.res:
                for tr in r:
                    alg = tr[1]
                    for m1, m2 in alg:
                        els1 = [el for el in d1.elements if m1 in el]   #??? zamiast d1.element powinna byc lista elementow z debugera? to samo z contaktami nizej?
                        els2 = [el for el in d2.elements if m2 in el]
                        tup = (str(d1.central_element.central_monomer.get_pdb_id()), str(d2.central_element.central_monomer.get_pdb_id()))
                        msg = "%s returned as aligned in %s-%s desc alignment, but not covered by any element"
                        self.assertGreater(len(els1), 0, msg % ((str(m1.get_pdb_id()),) + tup))
                        self.assertGreater(len(els2), 0, msg % ((str(m2.get_pdb_id()),) + tup))
                        conels1 = reduce(operator.add, [c.elements for c in d1.contacts])
                        conels2 = reduce(operator.add, [c.elements for c in d2.contacts])
                        self.assertTrue(any(e in conels1 for e in els1))
                        self.assertTrue(any(e in conels2 for e in els2))

        def test_tresholds(self):
            for r, d1, d2 in self.res:
                th = ConfigManager.compdesc.th_overall_rmsd
                try:
                    ls = [len(r[0][1])]
                except IndexError:
                    ls = [0]
                for diff, ind in [(-1., 0), (1., 3), (2., 4)]:
                    ConfigManager.compdesc.th_overall_rmsd = th + diff
                    nr = cpd.compdesc(d1, d2)
                    try:
                        ls.insert(ind, len(nr[0][1]))
                    except IndexError:
                        # raised when there is no compdes results
                        ls.insert(ind, 0)

                for tl1, tl2 in zip(ls[:-1], ls[1:]):
                    self.assertGreaterEqual(tl2, tl1)

                ConfigManager.compdesc.th_overall_rmsd = th

        """
        def test_continuity(self):

            def mk_star_desc(desc, resd):
                ""Arguments:
                desc -- basic desc.
                resd -- dict of desc with mers inds as values.
                ""
                cons = dict((tuple(sorted([e.central_monomer.ind for e in c.elements])), c) for c in desc.contacts)
                for e in desc.elements:
                    if e == desc.central_element: continue
                    try:
                        for c in resd[e.central_monomer.ind].contacts:
                            inds = tuple(sorted([e.central_monomer.ind for e in c.elements]))
                            cons[inds] = c
                    except KeyError:
                        continue
                return stc.ProteinDescriptor(desc.central_element, cons.values())

            nd1 = []
            nd2 = []
            for dummy, d1, d2 in self.res:
                nd1.append(mk_star_desc(d1, self.rdic1))
                nd2.append(mk_star_desc(d2, self.rdic2))

            for d1 in nd1[10:]:
                for d2 in nd2[10:]:
                    print d1, d2
#                    import pdb; pdb.set_trace()
                    print (cpd.compdesc(d1, d2), d1, d2)
#            nres = [(cpd.compdesc(d1, d2), d1, d2) for d1 in nd1 for d2 in nd2]
        """

        @classmethod
        def tearDownClass(cls):
            del cls.res
            del cls.s1
            del cls.s2
            del cls.ldr

    return BasicCompdescTests

def load_tests(loader, standard_tests, pattern):

    create = unittest.TestSuite()

    pairs = [('scop://d1m55a_', 'scop://d1p4da_')]

    for pair in pairs:
        create.addTests(loader.loadTestsFromTestCase(create_basic_tests(pair, fast)))

    standard_tests.addTests(create)

    return standard_tests

if __name__ == '__main__':
    if syntax_check.rip_argv('--skip-slow'):
        skip_slow = True

    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
