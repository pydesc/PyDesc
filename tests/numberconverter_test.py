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
Unit tests for numberconverter.py.

Usage:
    python numberconverter_test.py [-v]

    or

    python -m unittest [-v] numberconverter_test
"""

import syntax_check
from syntax_check import testing, notest
import unittest

import Bio.PDB

import tests
import pydesc.numberconverter as numberconverter

# pylint: disable=C0111

syntax_check.module = numberconverter

TestSyntax = syntax_check.module_syntax()

notest(numberconverter.convert_to_id)
notest(numberconverter.NumberConverter.get_list_of_inds)

data_dir = tests.__path__[0] + '/data/test_structures/'


@testing(numberconverter.NumberConverter)
class NumberConverterTest(unittest.TestCase):
    structures = ['1asz', '1gax', '1no5', '1pxq', '2dlc', '2lp2', '3ftk', '3g88', '3lgb', '3m6x', '3npn', '3tk0', '3umy']

    @classmethod
    def setUpClass(cls):
        cls.biopdb = {}
        cls.res_count = {}
        for s in cls.structures:
            cls.biopdb[s] = Bio.PDB.PDBParser(QUIET=True).get_structure(s, data_dir + '%s.pdb' % s)
            #~ cls.res_count[s] = len([res for res in cls.biopdb[s].get_residues() if not res.get_full_id()[3][0].startswith('W')])
            # TO test fails in case of NMR structure because all NMR models contributes to res_count
            cls.res_count[s] = len([res for res in cls.biopdb[s][0].get_residues() if not res.get_full_id()[3][0].startswith('W')])

    @classmethod
    def tearDownClass(cls):
        del cls.res_count
        del cls.biopdb

    @testing(numberconverter.NumberConverter.__init__)
    @testing(numberconverter.perform_Smith_Waterman)
    @testing(numberconverter.build_SW_matrix)
    @testing(numberconverter.go_backwards)
    @testing(numberconverter.PDB_id)
    @testing(numberconverter.PDB_id.__init__)
    @testing(numberconverter.PDB_id.from_pdb_residue)
    def test_init_file(self):
        for s in self.structures:
            #~ numberconverter.NumberConverter(s)
            # TO passing string to nc is not allowed
            numberconverter.NumberConverter(Bio.PDB.PDBParser(QUIET=True).get_structure(s, data_dir + '%s.pdb' % s))

    @testing(numberconverter.NumberConverter.__init__)
    def test_init_biopdb(self):
        for s in self.biopdb.values():
            numberconverter.NumberConverter(s)

    @testing(numberconverter.NumberConverter.get_ind)
    def test_get_ind(self):
        for (sn, s) in self.biopdb.items():
            nc = numberconverter.NumberConverter(s)

            for pdb_m in s:
                # TO returning twice the same ind is not an error in case of NMR structures
                # so I think test should check if same ind is not returned twice for one model
                # not for structure.

                l = []

                for res in pdb_m.get_residues():
                    res_id = res.get_full_id()
                    if not res_id[3][0].startswith('W'):
                        try:
                            pdb_id = list(res.get_id())
                            pdb_id[0] = res.parent.id
                            pdb_id[2] = pdb_id[2] if pdb_id[2] != ' ' else None
                            i = nc.get_ind(tuple(pdb_id))
                        except ValueError:
                            self.fail("Structure %s residue %s not found." % (sn, str(res_id)))
    #                    self.assertEqual(res_id, nc.get_pdb_id(i))
                        if i in l:
                            self.fail("Structure %s PDB index %d returned twice." % (sn, i))

                        l.append(i)

    @testing(numberconverter.NumberConverter.get_pdb_id)
    @testing(numberconverter.PDB_id.__str__)
    def test_get_pdb_id(self):
        for (sn, s) in self.biopdb.items():
            nc = numberconverter.NumberConverter(s)

            l = []
            sl = []

            for i in range(self.res_count[sn]):
                try:
                    pdb_id = nc.get_pdb_id(i + 1)
                except KeyError:
                    self.fail("Structure %s residue %d not found." % (sn, i + 1))

                if pdb_id in l:
                    self.fail("Structure %s PDB id %s returned twice." % (sn, str(pdb_id)))

                l.append(pdb_id)
                sl.append(str(pdb_id))

    @testing(numberconverter.NumberConverter.get_ind)
    @testing(numberconverter.NumberConverter.get_pdb_id)
    def test_get_reversibility(self):
        for (sn, s) in self.biopdb.items():
            nc = numberconverter.NumberConverter(s)

            for i in range(self.res_count[sn]):
                try:
                    pdb_id = nc.get_pdb_id(i + 1)
                except KeyError:
                    self.fail("Structure %s residue %d not found." % (sn, i + 1))

                try:
                    nc.get_ind(pdb_id)
                except ValueError:
                    self.fail("Structure %s residue %s not found." % (sn, str(pdb_id)))

    def try_pdb_id(self, pid):
        self.assertIsInstance(pid.chain, str)
        self.assertTrue(pid.icode is None or (isinstance(pid.icode, str) and len(pid.icode) == 1))
        self.assertIsInstance(pid.ind, int)

    @testing(numberconverter.PDB_id.chain)
    @testing(numberconverter.PDB_id.icode)
    @testing(numberconverter.PDB_id.ind)
    def test_PDB_id_properties(self):
        for s in self.biopdb.values():
            nc = numberconverter.NumberConverter(s)
            for i in nc.dict_ind_to_pdb.values():
                self.try_pdb_id(i)

    @testing(numberconverter.PDB_id.from_string)
    def test_PDB_id_from_string(self):
        for ch in ['a', 'A', ';', '%']:
            for i in range(201):
                for ic in [' ', 'A', 'B']:
                    self.try_pdb_id(numberconverter.PDB_id.from_string(ch+str(i)+ic))


if __name__ == '__main__':
    unittest.main()
