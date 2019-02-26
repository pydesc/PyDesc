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
Unit tests_legacy for alignment.py.

Usage:
    python alignment_test.py [-v] [--fast]

    or

    python -m unittest [-v] alignments_test
"""

import unittest

import syntax_check
from syntax_check import notest, testing, test_name_append

import pydesc.alignment as alignment
import pydesc.monomer as monomer
import pydesc.structure as structure
import pydesc.config as config
import pydesc.warnexcept as warnexcept
import tests_legacy

config.ConfigManager.warnings_and_exceptions.class_filters.set("LocalCopyAccess", "ignore")
config.ConfigManager.warnings_and_exceptions.class_filters.set("UnknownParticleName", "ignore")
config.ConfigManager.warnings_and_exceptions.class_filters.set("IncompleteChainableParticle", "ignore")
config.ConfigManager.warnings_and_exceptions.class_filters.set("Info", "ignore")


data_dir = tests_legacy.__path__[0] + '/data/test_alignments/'
file_extensions = ['csv', 'fasta', 'pal', 'xml']

syntax_check.module = alignment

TestSyntax = syntax_check.module_syntax()


@testing(alignment.Alignment)
class AlignmentTest(unittest.TestCase):

    @testing(alignment.AlignmentLoader)
    @testing(alignment.AlignmentLoader.load_alignment)
    @testing(alignment.AlignmentLoader.load_structure)
    @testing(alignment.AlignmentLoader.load_xml)
    @testing(alignment.AlignmentLoader.load_csv)
    @testing(alignment.AlignmentLoader.load_pal)
    @testing(alignment.AlignmentLoader.load_fasta)
    @testing(alignment.AlignmentLoader.__init__)
    @classmethod
    def setUpClass(cls):
        cls.loader = alignment.AlignmentLoader()
        cls.als = []
        for ex in file_extensions:
            cls.als.append(cls.loader.load_alignment(data_dir + 'ory.' + ex))

    @testing(alignment.Alignment.build_from_pair_alignments)
    @testing(alignment.Alignment.build_from_list_of_mers)
    @testing(alignment.PairAlignment.__init__)
    @testing(alignment.AlignedTuples.__init__)
    @testing(alignment.MergedPairAlignments.__init__)
    def test_build(self):

        class TestStructure(structure.AbstractStructure):
            def __init__(self):
                pass

        class TestMer(monomer.Monomer):
            def __init__(self, ind, struc):
                self.ind = ind
                self.structure = struc
            def _valid_indicators(instance):
                pass

        ts1 = TestStructure()
        ts2 = TestStructure()
        ts1._monomers = [TestMer(i, ts1) for i in range(100)]
        ts2._monomers = [TestMer(i, ts2) for i in range(100)]
        pa1 = alignment.PairAlignment([ts1, ts2], zip(ts1._monomers[:50], ts2._monomers[50:]))

    @testing(alignment.Alignment.get_selections)
    def test_get_selections(self):
        for al in self.als:
            mers = zip(*al.aligned_mers)
            sels = al.get_selections()
            for i, struc in enumerate(al.structures):
                us = sorted(list(sels[i].create_structure(struc)))
                mrs = sorted(list(set(filter(lambda x: x != '-', mers[i]))))
                self.assertEqual(us, mrs)

    @classmethod
    def tearDownClass(cls):
        del cls.als
        del cls.loader

@testing(alignment.PairAlignment)
class PairAlignmentTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.loader = alignment.AlignmentLoader()
        cls.als = []
        for ex in file_extensions:
            temp = cls.loader.load_alignment(data_dir + 'ory.' + ex).alignments
            temp2 = [alignment.PairAlignment(i.structures, i.aligned_mers) for i in temp]
            cls.als.extend(zip(temp2, temp))

    @testing(alignment.PairAlignment.__eq__)
    @testing(alignment.PairAlignment.__gt__)
    @testing(alignment.PairAlignment.__ge__)
    @testing(alignment.PairAlignment._prepare_comparison)
    def test_comparison(self):
        for al, alcp in self.als:
            self.assertTrue(al is not alcp)
            self.assertTrue(al == alcp)
            self.assertFalse(al > alcp)
            self.assertTrue(al >= alcp)
            self.assertFalse(al < alcp)
            self.assertTrue(al <= alcp)

            alcp.aligned_mers = alcp.aligned_mers[:-1]

            self.assertTrue(al > alcp)
            self.assertTrue(al >= alcp)
            self.assertFalse(al < alcp)
            self.assertFalse(al <= alcp)

            self.assertTrue(alcp <= al)
            self.assertFalse(alcp >= al)
            self.assertFalse(al < alcp)
            self.assertTrue(al > alcp)

    @classmethod
    def tearDownClass(cls):
        del cls.als
        del cls.loader

if __name__ == '__main__':
    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
