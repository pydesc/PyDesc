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
Unit tests_legacy for cydesc/compdesc.py.

Usage:
    python compdesc_test.py [-v] [--fast]

    or

    python -m unittest [-v] compdesc_test

created: 12.05.2014, Pawel Daniluk
"""

import unittest
import ctypes
import itertools
import math


import tests_legacy.syntax_check as syntax_check
from tests_legacy.syntax_check import testing, test_name_append, notest

import pydesc.structure as structure
import pydesc.numberconverter as numberconverter
import pydesc.cydesc as cydesc
import pydesc.cydesc.compdesc as compdesc
import pydesc.compat.dsc as dsc
import tests_legacy


class t_compdesc_debug_libdesc(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.

    """
    Class corresponding to t_compdesc_debug in desc.h (in libdesc).

    Contains intermediary results for debug purposes.
    """
    _fields_ = [('max_al_len', ctypes.c_int),
                ('n_el1', ctypes.c_int),
                ('n_el2', ctypes.c_int),
                ('elts1', ctypes.POINTER(ctypes.c_int)),
                ('elts2', ctypes.POINTER(ctypes.c_int)),
                ('th_n_el', ctypes.c_float),
                ('stage1', ctypes.c_int),
                ('n_seg_adj1', ctypes.c_int),
                ('n_seg_adj2', ctypes.c_int),
                ('th_n_seg_adj', ctypes.c_float),
                ('stage2', ctypes.c_int),
                ('zero_el_rmsd', ctypes.c_float),
                ('th_zero_el_rmsd', ctypes.c_float),
                ('stage3', ctypes.c_int),
                ('th_el_rmsd', ctypes.c_float),
                ('th_pair_rmsd', ctypes.c_float),
                ('n_comp', ctypes.c_int),
                ('comp1', ctypes.POINTER(ctypes.c_int)),
                ('comp2', ctypes.POINTER(ctypes.c_int)),
                ('el_rmsd', ctypes.POINTER(ctypes.c_float)),
                ('pair_rmsd', ctypes.POINTER(ctypes.c_float)),
                ('n_edges', ctypes.c_int),
                ('edge_list1', ctypes.POINTER(ctypes.c_int)),
                ('edge_list2', ctypes.POINTER(ctypes.c_int)),
                ('stage4', ctypes.c_int),
                ('edge_graph', ctypes.POINTER(ctypes.POINTER(ctypes.c_char))),
                ('n_sets', ctypes.c_int),
                ('indep_sets', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
                ('max_allowed', ctypes.c_int),
                ('n_allowed', ctypes.POINTER(ctypes.c_int)),
                ('allowed_comb', ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_int)))),
                ('n_align', ctypes.POINTER(ctypes.c_int)),
                ('allowed_align', ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_int)))),
                ('stage5', ctypes.c_int),
                ('n_full_align', ctypes.c_int),
                ('full_align', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
                ('full_align_cov', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
                ('full_align_rmsd', ctypes.POINTER(ctypes.c_float)),
                ('stage6', ctypes.c_int),
                ('full_align_crop_val_run1', ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
                ('full_align_crop_max_run1', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
                ('full_align_run1', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
                ('full_align_rmsd_run1', ctypes.POINTER(ctypes.c_float)),
                ('n_res_align', ctypes.c_int),
                ('res_align', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
                ('stage7', ctypes.c_int)]

    def print_full_align(self, n):
        self.print_alignment(self.full_align[n])

    def print_res_align(self, n):
        self.print_alignment(self.res_align[n])

    def print_full_align_run1(self, n):
        self.print_alignment(self.full_align_run1[n])

    def print_alignment(self, al):
        pairs = zip(al[0:self.max_al_len], al[self.max_al_len:2 * self.max_al_len])

        print " ".join(["%4d" % p1 for (p1, p2) in pairs if p1 != 0 or p2 != 0])
        print " ".join(["%4d" % p2 for (p1, p2) in pairs if p1 != 0 or p2 != 0])


# This is not a constant.
libcompdesc_test = cydesc.load_library('compdesc_test')  # pylint: disable=C0103

libcompdesc_test.set_desc1.argtypes = [ctypes.c_char_p]
libcompdesc_test.set_desc2.argtypes = [ctypes.c_char_p]
libcompdesc_test.my_load_domain.argtypes = [ctypes.c_char_p]
libcompdesc_test.unpack_pars.argtypes = [ctypes.POINTER(compdesc.t_compdesc_parameters)]
libcompdesc_test.compdesc_ref.argtypes = [ctypes.c_int, ctypes.POINTER(compdesc.t_compdesc_parameters), ctypes.POINTER(ctypes.POINTER(compdesc.t_compdesc_result)), ctypes.POINTER(ctypes.POINTER(t_compdesc_debug_libdesc))]
libcompdesc_test.average_coords.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_float))), ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_int)))]

libcompdesc_test.unpack_pars(compdesc.t_compdesc_parameters())

syntax_check.module = compdesc

data_dir = tests_legacy.__path__[0] + '/data/dsc/'

fast = False

TestSyntax = syntax_check.module_syntax()

# pylint: disable=C0111, R0201, R0912, R0914

struct_dict = {}


def average_coords(struct):
    coord_p = ctypes.POINTER(ctypes.POINTER(ctypes.c_float))()
    seg_adj_p = ctypes.POINTER(ctypes.POINTER(ctypes.c_int))()
    libcompdesc_test.average_coords(struct.dsc_name, coord_p, seg_adj_p)

    return extract_array(coord_p, (len(struct), 7)), extract_array(seg_adj_p, (len(struct), len(struct)))


def compdesc_ref(desc1, desc2, max_res=0, debug=False):
    """
    Computes alignments of the desc2 to the desc1 using libdesc.

    Parameters:
        * desc1, desc2 - Descriptor objects equivalent to descriptors loaded in libdesc.
        * max_res -- max number of results to be returned (0 means any)
        * debug -- if set an instance of t_compdesc_debug is returned along with the resuls.

    Results are sorted by size (decreasing) and RMSD (increasing).
    """

    pars = compdesc.t_compdesc_parameters()
    results_p = ctypes.POINTER(compdesc.t_compdesc_result)()
    debug_p = ctypes.POINTER(t_compdesc_debug_libdesc)()

    n_results = libcompdesc_test.compdesc_ref(max_res, pars, results_p, debug_p)

    results = [res.unpack(desc1, desc2) for res in results_p[0:n_results]]

    compdesc.libcompdesc.free_results(results_p, n_results)

    if debug:
        return results, debug_p.contents
    else:
        return results

loaded_structs = {}


def load_struct(str_name):
    dscfile = dsc.DSCFile(open('%s/%s.dsc' % (data_dir, str_name)))
    struct = dsc.DSCStructure(dscfile)
    struct.dsc_name = str_name
    loaded_structs[str_name] = struct
    libcompdesc_test.my_load_domain('file://%s/%s.dsc@%s' % (data_dir, str_name, str_name))
    return struct


def prep_desc(desc_name):
    str_name, pdb_num = desc_name.split('#')

    pdb_id = numberconverter.PDB_id.from_string(pdb_num)

    try:
        struct = loaded_structs[str_name]
    except KeyError:
        struct = load_struct(str_name)

    ind = struct.converter.get_ind(pdb_id)
    desc = structure.AbstractDescriptor.build(structure.Element.build(struct[ind]))

    return desc


def extract_array(arr, sizes):
    if sizes == ():
        return arr
    else:
        return [extract_array(row, sizes[1:]) for row in arr[0:sizes[0]]]


def build_descs(structure_obj):

    res = []

    for mer in structure_obj:
        try:
            desc = structure.AbstractDescriptor.build(structure.Element.build(mer))
        except:
            pass
        else:
            res.append(desc)

    return res


def make_compdesctest(struct1, struct2):
    """ Create and return a CompDescTest test case for given descriptors. """

    def fix_pdb_id(id_string):
        return id_string.replace('?', '')

    @test_name_append(struct1.dsc_name, struct2.dsc_name)
    class CompDescTest(unittest.TestCase):
        def assertSubalignments(self, res1, res2, msg=''):
            report = ""

            if len(res1) != len(res2):
                report += "Different numbers of results PyDesc:%d, libdesc:%d.\n" % (len(res1), len(res2))

            for r2 in res2:
                rmsd2, al2, tr2 = r2
                ok = 0
                for r1 in res1:
                    rmsd1, al1, tr1 = r1
                    if al2 <= al1:
                        ok = 1
                        if al2 < al1:
                            report += "Alignment %d contained in %d but not equal.\n" % (res2.index(r2), res1.index(r1))
                        break

                if not ok:
                    self.fail(msg + "\n" + report + "Alignment %d not found.\n" % res2.index(r2))

            if report != "":
                print report

        @unittest.skip("Not relevant any more.")
        def test_seg_adj(self):
            map(self.do_seg_adj, (struct1, struct2))

        def do_seg_adj(self, struct):
            coords = [tuple(m.ca) + tuple(m.backbone_average) + (m.adjusted_length(), ) for m in struct]
            seg_adj = [[struct[m1.ind:m2.ind].adjusted_number() if isinstance(struct[m1.ind:m2.ind], structure.Segment) else 0 for m2 in struct] for m1 in struct]

            coords1, seg_adj1 = average_coords(struct)

            rnd2 = lambda x: None if x is None or math.isnan(x) else round(x, 2)

            coords = [map(rnd2, r) for r in coords]
            coords1 = [map(rnd2, r) for r in coords1]

            for l1, l2 in zip(seg_adj, seg_adj1):
                for n1, n2 in zip(l1, l2):
                    if n1 != 0 and n1 != n2:
                        self.assertTupleEqual(*zip(*[t for t in zip(l1, l2) if t[0] != 0]))

        def test_comp(self):
            descs1, descs2 = map(build_descs, (struct1, struct2))

            for desc1 in descs1[:]:
                print "%s#%s / %s" % (desc1.derived_from.dsc_name, fix_pdb_id(str(desc1.central_element.central_monomer.get_pdb_id())), fix_pdb_id(str(descs1[-1].central_element.central_monomer.get_pdb_id())))
                for desc2 in descs2[:]:
                    # print descs1.index(desc1), descs2.index(desc2)
                    self.do_comp(desc1, desc2)

        def do_comp(self, desc1, desc2):
            self.maxDiff = None
            self.longMessage = True

            desc1_name = "%s#%s" % (desc1.derived_from.dsc_name, fix_pdb_id(str(desc1.central_element.central_monomer.get_pdb_id())))
            desc2_name = "%s#%s" % (desc2.derived_from.dsc_name, fix_pdb_id(str(desc2.central_element.central_monomer.get_pdb_id())))

            # print desc1_name, desc2_name

            libcompdesc_test.set_desc1(desc1_name)
            libcompdesc_test.set_desc2(desc2_name)

            res1, deb1, raw_res1, cdesc1, cdesc2 = compdesc.compdesc(desc1, desc2, debug=True)
            res2, deb2 = compdesc_ref(desc1, desc2, debug=True)
            # if deb1.stage2:
            #     import ipdb; ipdb.set_trace() # BREAKPOINT

            msg = '%s vs %s' % (desc1_name, desc2_name)

            try:
                self.assertSubalignments(res1, res2)

            except AssertionError:
                import ipdb; ipdb.set_trace() # BREAKPOINT
                raise


    return CompDescTest


def load_tests(loader, standard_tests, pattern):
    """ Add tests_legacy created by make_* functions for all structures. Return a complete TestSuite. """

    pairs = [('d1nkla_', 'd1qdma1'),
             ('d1an9a1', 'd1npxa1'),
             ('d1ay9b_', 'd1b12a_'),
             ('d1b5ta_', 'd1k87a2'),
             ('d1crla_', 'd1edea_'),
             ('d1d5fa_', 'd1nd7a_'),
             ('d1dlia1', 'd1mv8a1'),
             ('d1gbga_', 'd1ovwa_'),
             ('d1ggga_', 'd1wdna_'),
             ('d1gsaa1', 'd2hgsa1'),
             ('d1hava_', 'd1kxfa_'),
             ('d1hcyb2', 'd1lnlb1'),
             ('d1jj7a_', 'd1lvga_'),
             ('d1jwyb_', 'd1puja_'),
             ('d1jwyb_', 'd1u0la2'),
             ('d1kiaa_', 'd1nw5a_'),
             ('d1l5ba_', 'd1l5ea_'),
             ('d1nlsa_', 'd2bqpa_'),
             ('d1nw5a_', 'd2adma1'),
             ('d1qasa2', 'd1rsya_'),
             ('d1qq5a_', 'd3chya_'),
             ('d2adma1', 'd2hmyb_'),
             ('d2bbma_', 'd4clna_')]

    if fast:
        pairs = pairs[:1]
        # pairs = [#('d1dlia1', 'd1mv8a1'),
        #          # ('d1b5ta_', 'd1k87a2'),
        #          # ('d1d5fa_', 'd1nd7a_'),
        #          ('d1hcyb2', 'd1lnlb1')
        #          ]
    basic = unittest.TestSuite()

    for pair in pairs:
        struct1, struct2 = map(load_struct, pair)

        basic.addTests(loader.loadTestsFromTestCase(make_compdesctest(struct1, struct2)))

    standard_tests.addTests(basic)

    return standard_tests

if __name__ == '__main__':
    if syntax_check.rip_argv('--fast'):
        fast = True

    unittest.main()
