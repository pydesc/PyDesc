# Copyright 2017 Pawel Daniluk, Tymoteusz Oleniecki
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
Unit tests_legacy for structure.py.

Usage:
    python structure_test.py [-v] [--skip-slow] [--fast]

    or

    python -m unittest [-v] structure_test
"""

import unittest
import tests_legacy
from tests_legacy.syntax_check import module_syntax
from tests_legacy.syntax_check import parse_args

import random
import operator
import Bio.PDB
import warnings
import sys
import os

import pydesc.structure as structure
import pydesc.selection as selection
import pydesc.monomer as monomer
import pydesc.numberconverter as numberconverter
import pydesc.warnexcept as warnexcept
import pydesc.config as config

config.ConfigManager.warnings_and_exceptions.class_filters.set("NoConfiguration", "ignore")
config.ConfigManager.warnings_and_exceptions.class_filters.set("LocalCopyAccess", "ignore")
config.ConfigManager.warnings_and_exceptions.class_filters.set("Info", "ignore")
config.ConfigManager.warnings_and_exceptions.class_filters.set("UnknownParticleName", "ignore")
config.ConfigManager.warnings_and_exceptions.class_filters.set("IncompleteChainableParticle", "ignore")

warnexcept.set_filters()

data_dir = os.path.join(os.path.abspath(os.path.dirname(tests_legacy.__file__)), 'data/test_structures/')
dcd_data_dir = os.path.join(os.path.abspath(os.path.dirname(tests_legacy.__file__)), 'data/dcd/')


class StructureLoaderTest(unittest.TestCase):

    """ Testcase for StructureLoader class. """

    structures = ['1asz', '1gax', '1no5', '1pxq', '2dlc', '2lp2',
                  '3ftk', '3g88', '3lgb', '3m6x', '3npn', '3tk0', '3umy']

    def test_loader(self):
        """ Test StructureLoader.load_structure. """
        loader = structure.StructureLoader()
        for name in self.structures:
            with warnings.catch_warnings(record=True):
                struct = loader.load_structure(name)
            f = loader.handler.get_file(name)[0]
            pdb_models = Bio.PDB.PDBParser(QUIET=True).get_structure(name, f)

            for s in struct:
                self.assertIsInstance(s, structure.Structure)

            self.assertEqual(len(struct), len(pdb_models))



def make_trajectorytest(strname):
    """Create and return a StructureTest testcase for a structures linked with dcd trajectories."""

    sl = structure.StructureLoader()

    class TrajectoryTest(unittest.TestCase):
        name = strname
        @classmethod
        def setUpClass(cls):
            cls.pdb_structure = sl.load_structure("mdm2", path=dcd_data_dir + cls.name + ".pdb")[0]

        def test_linking_trajectory(self):
            self.assertIsNone(self.pdb_structure.frame)
            self.pdb_structure.link_dcd_file(dcd_data_dir + self.name + ".dcd")
            fr = self.pdb_structure.frame
            self.pdb_structure.next_frame()
            self.assertEqual(fr + 1, self.pdb_structure.frame)
            self.assertIsNotNone(self.pdb_structure.frame)
            self.pdb_structure.disconnect_trajectory()
            self.assertIsNone(self.pdb_structure.frame)

    return TrajectoryTest


def make_structure_basic_test(str_name, mer_type):
    """Create and return a StructureBasicTest test case for a given structure."""

    m1 = mer_type

    class StructureBasicTest(unittest.TestCase):
        name = str_name
        mer_type = m1

        def test_init(self):
            pdbV_structure = Bio.PDB.PDBParser(
                QUIET=True).get_structure(self.name, os.path.join(data_dir,
                                                                  "%s_only" % self.mer_type,
                                                                  '%s.pdb' % self.name)
                                          )
            converter = numberconverter.NumberConverter(pdb_structure)
            for model in pdb_structure:
                with warnings.catch_warnings(record=True) as wlist:
                    warnings.simplefilter("always")
                    structure.Structure(model, converter, )
                    syntax_check.warning_message(self, wlist)

    return StructureBasicTest


def make_structuretest(strname):
    """Create and return a StructureTest testcase for a given structure."""
    class StructureTest(unittest.TestCase):
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

        def test_len(self):
            #~ self.assertEqual(len(self.struct), sum(map(len, self.model)))
            #TO water particles were taken under consideration in last version
            self.assertEqual(len(self.struct), sum(map(len, [filter(lambda x: True if x.get_id()[0] != 'W' else False, chain) for chain in self.model])))

        def test_ss(self):
            self.struct.set_secondary_structure()
            ss = self.struct.get_secondary_structure()
            sss = self.struct.get_simple_secondary_structure()
            self.assertEqual(len(self.struct), len(ss))
            self.assertEqual(len(self.struct), len(sss))
            self.assertLessEqual(len(ss), 6)
            self.assertLessEqual(len(sss), 3)
            self.assertTrue(all(i in ('E', 'G', 'H', '-', 'S', 'T') for i in set(ss)))
            self.assertTrue(all(i in ('E', 'C', 'H') for i in set(sss)))

        def test_iter(self):
            l = list(iter(self.struct))

            self.assertEqual(len(l), len(self.struct))
            self.assertEqual(len(l), len(set(l)))

        def test_monomer_ind(self):
            for m in self.struct:
                pdb_id = numberconverter.PDB_id.create_from_pdb_residue(m.pdb_residue)
                self.assertEqual(m.ind, self.converter.get_ind(tuple(pdb_id)))

        def test_getitem(self):
            for m in self.struct:
                self.assertEqual(self.struct[m.ind], m)

        def test_continuity(self):
            for mer in self.struct:
                if isinstance(mer, monomer.MonomerChainable):
                    if mer.next_mer:
                        self.assertEqual(mer, mer.next_mer.previous_mer)
                    if mer.previous_mer:
                        self.assertEqual(mer, mer.previous_mer.next_mer)

        def test_getsliceunbound(self):
            whole = list(self.struct)
            pref = []
            for m in self.struct:
                slc = self.struct[m.ind:]
                slc1 = self.struct[:m.ind]
                self.assertItemsEqual(whole, pref + list(slc))
                self.assertEqual(len(whole), len(pref) + len(slc))
                self.assertItemsEqual(pref + [m], list(slc1))
                # TO structure slices take slice.stop mer!
                self.assertEqual(len(pref) + 1, len(slc1))
                pref.append(m)

        @unittest.skipIf(cla.fast, "Takes too long")
        def test_getslicebound(self):
            step = 20
            step1 = 10
            whole = list(self.struct)
            for sg in range(0, len(whole), step):
                s = sg + random.randint(0, min(step - 1, len(whole) - sg - 1))
                for eg in range(s, len(whole), step1):
                    e = eg + random.randint(0, min(step1 - 1, len(whole) - eg - 1))
                    slc = self.struct[whole[s].ind:whole[e].ind]
                    self.assertItemsEqual(whole[s:e + 1], list(slc))

                    seg = True
                    for m1, m2 in zip(list(slc)[:-1], list(slc)[1:]):
                        with warnings.catch_warnings(record=True):
                            #~ if not m1.is_next(m2):
                                #~ seg = False
                                #~ self.assertIsNone(slc.next_mer(m1))
                            #~ else:
                                #~ self.assertEqual(slc.next_mer(m1), m2)
                            # TO is_next is deprecated, so:
                            if len(slc) == 1:
                                seg = True if isinstance(slc[0], monomer.MonomerChainable) else False
                                break
                            if m1.next_mer == m2:
                                self.assertEqual(slc.next_mer(m1), m2)
                            else:
                                seg = False
                                self.assertIsNone(slc.next_mer(m1))
                    if seg:
                        self.assertIsInstance(slc, structure.Segment, "Slice is Userstructure, should be Segment: structure %s, slice: from %s to %s" % (str(self.struct), str(slc[0]), str(slc[-1])))

        def test_add_in(self):
            whole = list(self.struct)

            slices = []
            for dummy_i in range(10):
                s = random.randint(0, len(whole) - 1)
                e = random.randint(s, len(whole) - 1)
                slices.append(self.struct[whole[s].ind:whole[e].ind])

            for dummy_i in range(20):
                k = random.randint(2, 10)
                sample = random.sample(slices, k)
                res = reduce(operator.add, sample)
                resset = reduce(operator.or_, map(set, sample))

                self.assertEqual(len(set(res)), len(list(res)), "Duplicate mers after add.")
                self.assertSetEqual(set(res), resset)

                seg = True
                for m1, m2 in zip(list(res)[:-1], list(res)[1:]):
                    with warnings.catch_warnings(record=True):
                        #~ if not m1.is_next(m2):
                            #~ seg = False
                            #~ break
                        # TO is_next is deprecated
                        if not m1.next_mer == m2:
                            seg = False
                            break
                if seg:
                    self.assertIsInstance(res, structure.Segment)

        @unittest.skipIf(cla.fast, "Takes too long")
        def test_contactmap_and_contact_values(self):
            self.assertRaises(AttributeError, self.struct.get_contact_map)
            with warnings.catch_warnings(record=True) as wlist:
                warnings.simplefilter("always")
                self.struct.set_contact_map()

                syntax_check.warning_message(self, wlist)

            self.struct.get_contact_map()

            inds = map(operator.attrgetter('ind'), self.struct._monomers)

            config.ConfigManager.element.set("element_chainable_length", 3)

            for ind1 in inds:
                for ind2 in inds:
                    if ind1 == ind2 or self.struct[ind1].next_mer == None or self.struct[ind2].next_mer == None or self.struct[ind1].previous_mer == None or self.struct[ind2].previous_mer == None:
                        continue
                    c = structure.Contact(*[structure.ElementChainable(self.struct[i]) for i in (ind1, ind2)])
                    if ind2 in self.struct.contact_map.contacts[ind1]:
                        self.assertEqual(c.value, self.struct.contact_map.contacts[ind1][ind2])
                    else:
                        self.assertEqual(c.value, 0)

        def test_getsequence(self):
            for ch in self.struct.chains:
                seq = ch.get_sequence()
                chainable = filter(lambda x: isinstance(x, monomer.MonomerChainable), ch)
                self.assertEqual(len(seq), len(chainable), "Wrong len of sequence of %s chain" % ch.chain_char)

        def test_select(self):
            self.struct.select()

        def test_segment(self):
            whole = list(self.struct)

            nseg = 0
            cnt = 0

            while cnt < 1000 and nseg < 100:
                s = random.randint(0, len(whole) - 2)
                e = random.randint(s + 1, len(whole) - 1)
                slc = self.struct[whole[s].ind:whole[e].ind]

                cnt += 1
                if isinstance(slc, structure.Segment):
                    nseg += 1

                    slc.select()
                    self.assertEqual(whole[s], slc.start)
                    self.assertEqual(whole[e], slc.end)

            self.assertTrue(nseg, "Random segment generation failed to generate any segments.")

        def test_element(self):
            chs = [list(ch) for ch in self.struct.chains]
            whole = reduce(operator.add, chs)
            valerr = False
            msg = []

            config.ConfigManager.element.set("element_chainable_length", 5)

            for (i, m) in enumerate(whole):
                if isinstance(m, monomer.MonomerChainable):
                    willfail = False
                    mer_ch = max([ch for ch in chs if m in ch])
                    index = mer_ch.index(m)
                    if index < 2 or index >= len(mer_ch) - 2:
                        willfail = True
                    else:
                        try:
                            s = mer_ch[index - 2]
                            e = mer_ch[index + 2]
                        except IndexError:
                            willfail = True
                        else:
                            for m1, m2 in zip(mer_ch[(index - 2):(index + 2)], mer_ch[(index - 1):(index + 3)]):
                                with warnings.catch_warnings(record=True):
                                    if m1.next_mer == m2:
                                        continue
                                    willfail = True
                                    break

                    tests_performed['ElementChainable'] = True
                    if willfail:
                        try:
                            self.assertRaises(ValueError, lambda: structure.ElementChainable(m))
                        except:
                            msg.append((sys.exc_info(), m))
                            valerr = True
                    else:
                        el = structure.ElementChainable(m)
                        elb = structure.Element.build(m)
                        self.assertEqual(type(el), type(elb))
                        self.assertEqual(len(el), config.ConfigManager.element.element_chainable_length)
                        self.assertEqual(el.start, s)
                        self.assertEqual(el.end, e)

                elif isinstance(m, monomer.MonomerOther):
                    tests_performed['ElementOther'] = True
                    el = structure.ElementOther(m)

            if valerr:
                self.fail("Unexpected exceptions caught for mer %i:\n %s" % (msg[0][1].ind, str(msg[0][0])))

        @unittest.expectedFailure
        def test_contact_future(self):
            # test fails because of lack of pydesc.contacts imported
            # their are not imported in order to avoid false-positive result
            # of set_contact_map test in presence of import errors
            whole = list(self.struct)

            for i in range(100):
                r1 = random.choice(whole)
                r2 = random.choice(whole)

                if r1 == r2:
                    continue

                try:
                    el1 = structure.ElementChainable(r1)
                    el2 = structure.ElementChainable(r2)
                except:
                    continue

                def check(c, v, crit):
                    self.assertEqual((tuple(map(lambda x: x.central_monomer, c.elements)), c.value, c.criterion), ((r1, r2), v, crit))

                c = structure.Contact(r1, r2)
                check(c, None, None)

                val = syntax_check.expando()
                c = structure.Contact(r1, r2, value=val)
                check(c, val, None)

                cc = contacts.ContactCriterion()
                c = structure.Contact(r1, r2, contact_criterion=cc)
                check(c, None, cc)

                self.assertEqual(c - el1, el2)
                self.assertEqual(c - el2, el1)

                self.assertTrue(isinstance(c.select(), selection.SelectionsUnion) or isinstance(c.select(), selection.Range))

        def test_contact(self):
            whole = list(self.struct)

            for i in range(100):
                r1, r2 = random.sample(whole, 2)

                try:
                    el1 = structure.ElementChainable(r1)
                    el2 = structure.ElementChainable(r2)
                except:  # pylint: disable=W0702
                    continue

                def check(c, v, crit):
                    #~ self.assertEqual((tuple(map(lambda x: x.central_monomer, c.elements)), c.value, c.criterion), ((r1, r2), v, crit))
                    # TO contact no longer have criterion attr, and value property doesn't work in this test due to lack of mers' derived_from attr named contact_map
                    # which is coused by test design
                    self.assertEqual(sorted(map(operator.attrgetter('central_monomer'), c.elements)), sorted([r1, r2]))

                c = structure.Contact(el1, el2)
                check(c, None, None)

                self.assertEqual(c.get_other_element(el1).central_monomer.ind, el2.central_monomer.ind)
                self.assertEqual(c.get_other_element(el2).central_monomer.ind, el1.central_monomer.ind)

                sc = c.select()
                self.assertEqual(len(sc.specify(c).ids), len(c))

        def test_adjusted_number(self):
            whole = list(self.struct)

            cnt = 0

            while cnt < 100:
                s = random.randint(0, len(whole) - 2)
                e = random.randint(s + 1, len(whole) - 1)
                slc = self.struct[whole[s].ind:whole[e].ind]

                cnt += 1
                if isinstance(slc, structure.Segment):
                    adj = slc.adjusted_number()
                    inds = map(operator.attrgetter('ind'), slc)
                    for i in inds[1:]:
                        self.assertTrue(slc[inds[0]:i].adjusted_number() <= adj)
                    for i in inds[:-1]:
                        self.assertTrue(slc[i:inds[-1]].adjusted_number() <= adj)
                else:
                    conts = 0
                    for i, j in zip(slc[:-2], slc[-len(slc) + 1:]):
                        if i.next_mer == j:
                            continue
                        if isinstance(i, monomer.MonomerChainable):
                            conts += 1
                    self.assertTrue(conts <= slc.adjusted_number())

        def test_pdbstring(self):
            stg = self.struct.create_pdb_string()
            self.struct.save_pdb('stc_test.pdb')
            nstc = structure.StructureLoader().load_structure('test', path='stc_test.pdb')[0]
            stg2 = nstc.create_pdb_string()
            self.assertEqual(stg.read(), stg2.read())

    return StructureTest


def make_descriptortest(strname):
    """Create and return a StructureTest testcase for a given structure."""

    @unittest.skipIf(cla.fast, "Takes too long")
    class DescriptorTest(unittest.TestCase):
        name = strname

        @classmethod
        def setUpClass(cls):
            cls.pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(cls.name, data_dir + '%s.pdb' % cls.name)
            cls.converter = numberconverter.NumberConverter(cls.pdb_structure)
            cls.model = random.choice(cls.pdb_structure)
            with warnings.catch_warnings(record=True):
                cls.struct = structure.Structure(cls.model, cls.converter)
                #~ crit = contacts.ContactsAlternative(contacts.CaContact(), contacts.CbxContact(), contacts.RaContact())
                #~ cls.struct.set_contact_map(crit)
                cls.struct.set_contact_map()

        @classmethod
        def tearDownClass(cls):
            del cls.struct
            del cls.model
            del cls.converter
            del cls.pdb_structure

        @unittest.expectedFailure
        def test_descriptor_future(self):
            whole = list(self.struct)
            valerr = False
            msg = []

            for (i, m) in enumerate(whole):
                if isinstance(m, monomer.Nucleotide):
                    desc_class = structure.NucleotideDescriptor
                    name = 'NucleotideDescriptor'
                elif isinstance(m, monomer.Residue):
                    desc_class = structure.ProteinDescriptor
                    name = 'ProteinDescriptor'
                else:
                    continue

                willfail = False

                try:
                    structure.ElementChainable(m)
                except:
                    willfail = True

                tests_performed[name] = True
                if willfail:
                    try:
                        self.assertRaises(warnexcept.DiscontinuityError, lambda: desc_class.build(m))
                    except:
                        msg.append(sys.exc_info())
                        valerr = True
                    continue

                desc_class.build(m)

            if valerr:
                self.fail("Unexpected exceptions caught:\n %s" % str(msg))

        def test_descriptor(self):
            whole = list(self.struct)

            for (i, m) in enumerate(whole):
                if isinstance(m, monomer.Nucleotide):
                    desc_class = structure.NucleotideDescriptor
                    name = 'NucleotideDescriptor'
                elif isinstance(m, monomer.Residue):
                    desc_class = structure.ProteinDescriptor
                    name = 'ProteinDescriptor'
                else:
                    continue

                try:
                    el = structure.ElementChainable(m)
                except:
                    continue

                tests_performed[name] = True

                if self.struct.contact_map.get_monomer_contacts(m) == []:
                    self.assertRaises(ValueError, lambda: desc_class.build(el))
                    continue
                staph = False
                for i in self.struct.contact_map.get_monomer_contacts(m):
                    try:
                        structure.Contact(i[0], i[1], self.struct)
                    except:
                        staph = True
                        break
                    # TO sometimes there are contcts in contact map, but mers contacted with central monomer are to close to edge to create elements
                if staph:
                    continue
                desc = desc_class.build(el, self.struct.contact_map)
                absdesc = structure.AbstractDescriptor.build(el)
                self.assertEqual(type(desc), type(absdesc))
                self.assertEqual(len(desc), len(absdesc))

                def el_in_desc(el, desc):
                    for m in el:
                        self.assertIn(m, desc)

                el_in_desc(el, desc)

                self.assertEquals(list(desc.central_element), list(el))

                secondary_contacts = filter(lambda x: False if el.central_monomer.ind in [i.central_monomer.ind for i in x.elements] else True, desc.contacts)
                primary_inds = [i.central_monomer.ind for i in reduce(operator.add, [i.elements for i in desc.contacts if i not in secondary_contacts])]
                for c in desc.contacts:
                    try:
                        el_in_desc(c.get_other_element(el), desc)
                    except ValueError:
                        self.assertTrue(c in secondary_contacts, '%s contains elements settled by mers that are not in contact with central monomer' % str(desc))
                        self.assertTrue(any(el_.central_monomer.ind in primary_inds for el_ in c.elements), '%s elements are not on a list of primary elements of %s' % (str(c), str(desc)))

                for s in desc.segments:
                    self.assertIsInstance(s, structure.Segment)
                    el_in_desc(s, desc)

                for c in desc.contacts:
                    for el in c.elements:
                        if isinstance(el, structure.ElementChainable):
                            ok = False
                            for s in desc.segments:
                                if set(el).issubset(s):
                                    ok = True
                                    break
                            self.assertTrue(ok, "Element %s is not contained in any segment" % str(el))

                if len(desc.contacts):
                    self.assertSetEqual(set(desc._monomers), set(reduce(operator.add, desc.contacts)))
                else:
                    self.assertSetEqual(set(desc._monomers), set(el))
                if isinstance(desc, structure.ProteinDescriptor):
                    self.assertSetEqual(set(desc._monomers), set(reduce(operator.add, desc.segments)))
                else:
                    self.assertSetEqual(set(filter(lambda i: isinstance(i, monomer.MonomerChainable), desc._monomers)), set(reduce(operator.add, desc.segments)))

                self.assertEqual(len(list(desc)), len(set(desc)))


                desc.select()

    return DescriptorTest


def load_tests(loader, standard_tests, pattern):
    """ Add tests_legacy created by make_* functions for all structures. Return a complete TestSuite. """
    structures = {
        'dna': ['411D', '1EL2', '279D', '1X2Y', '1AG5', '2QEG', '400D', '1II1', '1EXL', '179D', '1DUF', '1D89',
                  '1VT9', '2N0Q', '2K0V', '2LSZ', '1NAJ', '1ZYH', '1F6I'],
        'prots': ['5MPV', '4ZTD', '5ERB', '5X55', '4ONK', '3NPU', '5IFH', '1A24', '3J96', '2JRM', '3PSC', '3G67',
               '4NJ6', '2BLL', '3BIP', '4LTT', '4YYN', '2LJP', '5LF9', '3FV6', '1KAN'],
        'rna': ['3J2E', '2H49', '1OW9', '2JXQ', '4A4S', '2MFD', '1FL8', '2KBP', '1KP7', '2MEQ', '2EUY', '2A0P', '1WKS', '2LPS', '2JYJ', '2AGN', '1KIS', '1A3M', '1SDR']
    }
    dcd_structures = ['mdm2']

    if cla.fast:
        structures = {
            'dna': ['411D', '1EL2', '279D', ],
            'prots': ['5MPV', '4ZTD', '5ERB', ],
            'rna': ['3J2E', '2H49', '1OW9', ]
        }

    basic = unittest.TestSuite()

    for mer_type, files_lst in structures.items():
        for file in files_lst:
            basic.addTests(loader.loadTestsFromTestCase(make_structure_basic_test(file, mer_type)))
            basic.addTests(loader.loadTestsFromTestCase(make_structuretest(file)))
            basic.addTests(loader.loadTestsFromTestCase(make_descriptortest(file)))

    for name in dcd_structures:
        basic.addTests(loader.loadTestsFromTestCase(make_trajectorytest(name)))

    standard_tests.addTests(basic)

    standard_tests.addTests(loader.loadTestsFromTestCase(module_syntax(structure)))

    return standard_tests

if __name__ == '__main__':
    cla, sys.argv[1:] = parse_args()

    unittest.main()
