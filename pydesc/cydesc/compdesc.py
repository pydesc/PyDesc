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
Interface to compdesc.


created: 08.05.2014 - Pawel Daniluk
"""

import ctypes
from ctypes import POINTER

import pydesc.alignment as alignment

import pydesc.cydesc as cydesc
import pydesc.cydesc.overfit as overfit
from pydesc.config import ConfigManager

ConfigManager.new_branch("compdesc")
ConfigManager.compdesc.set_default("ca_contact_distance", 6.0)
ConfigManager.compdesc.set_default("th_zero_el_rmsd", 2.0)
ConfigManager.compdesc.set_default("th_el_rmsd", 2.5)
ConfigManager.compdesc.set_default("th_pair_rmsd", 2.5)
ConfigManager.compdesc.set_default("force_order", 0)
ConfigManager.compdesc.set_default("th_mer_align", 0.33)
ConfigManager.compdesc.set_default("th_el_align", 0.5)
ConfigManager.compdesc.set_default("th_con_align", 0)
ConfigManager.compdesc.set_default("th_align_size_equiv", 3)
ConfigManager.compdesc.set_default("th_overall_rmsd", 2.5)


# This is not a constant.
libcompdesc = cydesc.load_library("compdesc")  # pylint: disable=C0103


class t_compdesc_parameters(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.

    """
    Class corresponding to t_compdesc_paramerers in compdesc.h.

    Contains parameters for descriptor comparison.
    """
    _fields_ = [
        ("th_zero_el_rmsd", ctypes.c_float),
        ("th_el_rmsd", ctypes.c_float),
        ("th_pair_rmsd", ctypes.c_float),
        ("force_order", ctypes.c_int),
        ("th_mer_align", ctypes.c_float),
        ("th_el_align", ctypes.c_float),
        ("th_con_align", ctypes.c_float),
        ("th_align_size_equiv", ctypes.c_float),
        ("th_overall_rmsd", ctypes.c_float),
    ]

    def __init__(self):  # pylint: disable=W0231
        """
        t_compdesc_parameters is created with sensible default values.
        """
        # __init__ supplied by ctypes.Structure should not be called, if there
        # is an __init__ supplied in a subclass.

        self.th_zero_el_rmsd = ConfigManager.compdesc.th_zero_el_rmsd
        self.th_el_rmsd = ConfigManager.compdesc.th_el_rmsd
        self.th_pair_rmsd = ConfigManager.compdesc.th_pair_rmsd
        self.force_order = ConfigManager.compdesc.force_order
        self.th_mer_align = ConfigManager.compdesc.th_mer_align
        self.th_el_align = ConfigManager.compdesc.th_el_align
        self.th_con_align = ConfigManager.compdesc.th_con_align
        self.th_align_size_equiv = ConfigManager.compdesc.th_align_size_equiv
        self.th_overall_rmsd = ConfigManager.compdesc.th_overall_rmsd


class t_map(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [
        ("len1", ctypes.c_int),
        ("len2", ctypes.c_int),
        ("map12", POINTER(ctypes.c_int)),
        ("map21", POINTER(ctypes.c_int)),
    ]


class t_desc_coverage(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [
        ("n_mer", ctypes.c_int),
        ("n_mer_all", ctypes.c_int),
        ("n_el", ctypes.c_int),
        ("n_el_all", ctypes.c_int),
        ("n_con", ctypes.c_int),
        ("n_con_all", ctypes.c_int),
        ("n_AA", ctypes.c_int),
        ("n_contacts", ctypes.c_int),
        ("desc_vec", POINTER(ctypes.c_int)),
        ("con_vec", POINTER(ctypes.c_int)),
    ]


class t_al_coverage(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [("desc1", t_desc_coverage), ("desc2", t_desc_coverage)]


class t_alignment(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [
        ("mer_map", POINTER(t_map)),
        ("el_map", POINTER(t_map)),
        ("contact_map", POINTER(t_map)),
        ("mer_cov_count", POINTER(ctypes.c_int)),
        ("el_cov_count", POINTER(ctypes.c_int)),
        ("desc1", POINTER(cydesc.CDescriptor)),
        ("desc2", POINTER(cydesc.CDescriptor)),
        ("n_AA", ctypes.c_int),
        ("coverage", POINTER(t_al_coverage)),
        ("RMSD", ctypes.c_float),
        ("TR", cydesc.overfit.t_transrot),
    ]


class t_elsim(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [
        ("el1", POINTER(cydesc.CElement)),
        ("el2", POINTER(cydesc.CElement)),
        ("RMSD", ctypes.c_float),
        ("overfit_sums", cydesc.overfit.t_overfit_sums),
        ("consim_cnt", ctypes.c_int),
        ("core_cnt", ctypes.c_int),
        ("pos", ctypes.c_int),
    ]


class t_consim(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [
        ("con1", POINTER(cydesc.CContact)),
        ("con2", POINTER(cydesc.CContact)),
        ("elsim1", POINTER(t_elsim)),
        ("elsim2", POINTER(t_elsim)),
        ("RMSD", ctypes.c_float),
        ("pos", ctypes.c_int),
    ]


class t_sim_components(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [
        ("elsim_list", POINTER(t_elsim)),
        ("consim_list", POINTER(t_consim)),
        ("consim_map", POINTER(POINTER(POINTER(t_consim)))),
        ("n_elsim", ctypes.c_int),
        ("n_consim", ctypes.c_int),
    ]


class t_sim(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [("elsim", POINTER(t_elsim)), ("consim", POINTER(t_consim))]


class t_sim_graph(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    _fields_ = [
        ("sim_list", POINTER(t_sim)),
        ("sim_graph", POINTER(POINTER(ctypes.c_char))),
        ("n_sim", ctypes.c_int),
    ]


class t_compdesc_debug(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.

    """
    Class corresponding to t_compdesc_debug in compdesc.h.

    Contains intermediary results for debug purposes.
    """
    _fields_ = [
        ("n_el1", ctypes.c_int),
        ("n_el2", ctypes.c_int),
        ("elts1", POINTER(ctypes.c_int)),
        ("elts2", POINTER(ctypes.c_int)),
        ("stage0", ctypes.c_int),
        ("prelim_cov", POINTER(t_al_coverage)),
        ("stage1", ctypes.c_int),
        ("zero_el_rmsd", ctypes.c_float),
        ("n_all_element_comp", ctypes.c_int),
        ("all_element_comp", POINTER(t_elsim)),
        ("n_all_contact_comp", ctypes.c_int),
        ("all_contact_comp", POINTER(t_consim)),
        ("sim_components", POINTER(t_sim_components)),
        ("sim_graph", POINTER(t_sim_graph)),
        ("stage2", ctypes.c_int),
        ("n_sets", ctypes.c_int),
        ("indep_sets", POINTER(POINTER(ctypes.c_int))),
        ("n_allowed", POINTER(ctypes.c_int)),
        ("allowed_comb_per_set", POINTER(POINTER(POINTER(ctypes.c_int)))),
        ("n_full_align", ctypes.c_int),
        ("allowed_comb", POINTER(POINTER(ctypes.c_int))),
        ("full_align", POINTER(t_alignment)),
        ("stage3", ctypes.c_int),
        ("full_align_crop_val_run1", POINTER(POINTER(ctypes.c_float))),
        ("full_align_crop_max_run1", POINTER(POINTER(ctypes.c_int))),
        ("full_align_run1", POINTER(t_alignment)),
        ("n_res_align", ctypes.c_int),
        ("res_align", POINTER(t_alignment)),
        ("stage4", ctypes.c_int),
    ]


class t_compdesc_result(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.

    """
    Class corresponding to t_compdesc_result in compdesc.h.

    Contains a single descriptor assignment.

    Probably this class should be merged with pydesc.cydesc.fitdest.t_fitdesc_result.
    """
    _fields_ = [
        ("n_monomers", ctypes.c_int),
        ("desc1_monomers", POINTER(ctypes.c_int)),
        ("desc2_monomers", POINTER(ctypes.c_int)),
        ("n_elements", ctypes.c_int),
        ("desc1_elements", POINTER(ctypes.c_int)),
        ("desc2_elements", POINTER(ctypes.c_int)),
        ("n_contacts", ctypes.c_int),
        ("desc1_contacts", POINTER(ctypes.c_int)),
        ("desc2_contacts", POINTER(ctypes.c_int)),
        ("RMSD", ctypes.c_float),
        ("TR", overfit.t_transrot),
    ]

    def unpack(self, desc1, desc2):
        """
        Returns a triple containing RMSD value, PairAlignment object and a TRTMatrix.

        desc1 and desc2 are AbstractDescriptor instances which the assignment refers to.
        """
        n_mers = self.n_monomers
        get_monomers = lambda struct, array: list(
            map(struct.__getitem__, array[0:n_mers])
        )
        pair_al = alignment.PairAlignment(
            [desc1, desc2],
            list(
                zip(
                    get_monomers(desc1, self.desc1_monomers),
                    get_monomers(desc2, self.desc2_monomers),
                )
            ),
        )
        return self.RMSD, pair_al, self.TR.to_trtmatrix()


if libcompdesc.has_debug():
    libcompdesc.compdesc.argtypes = [
        ctypes.POINTER(cydesc.CDescriptor),
        ctypes.POINTER(cydesc.CDescriptor),
        ctypes.c_int,
        ctypes.POINTER(t_compdesc_parameters),
        ctypes.POINTER(ctypes.POINTER(t_compdesc_result)),
        ctypes.POINTER(ctypes.POINTER(t_compdesc_debug)),
    ]
else:
    libcompdesc.compdesc.argtypes = [
        ctypes.POINTER(cydesc.CDescriptor),
        ctypes.POINTER(cydesc.CDescriptor),
        ctypes.c_int,
        ctypes.POINTER(t_compdesc_parameters),
        ctypes.POINTER(ctypes.POINTER(t_compdesc_result)),
    ]

libcompdesc.free_results.argtypes = [ctypes.POINTER(t_compdesc_result), ctypes.c_int]


def _desc_getter(num):
    """
    Returns a property getter for nth descriptor (num in [1,2])
    """

    def get_desc(self):
        """
        Descriptor being compared.
        """
        return getattr(self, "_desc%d" % num)

    return get_desc


def _desc_setter(num):
    """
    Returns a property setter for nth descriptor (num in [1,2])
    """

    def set_desc(self, desc_obj):
        """
        Descriptor setter. Creates a corresponding CDescriptor object.
        """

        setattr(self, "_desc%d" % num, desc_obj)
        setattr(
            self,
            "_c_desc%d" % num,
            None if desc_obj is None else cydesc.CDescriptor(desc_obj),
        )

    return set_desc


class CompDesc(object):

    """
    Computation of alignment between descriptors.


    Allows storing descriptors in order to avoid overhead of creating
    corresponding structures in C many times.
    """

    def __init__(self, desc1=None, desc2=None):
        """
        Initializator. Can set desc1 and desc2 attributes.
        """

        # Nicely pacifies Pylint
        self._desc1 = None
        self._c_desc1 = None
        self._desc2 = None
        self._c_desc2 = None

        self.desc1 = desc1
        self.desc2 = desc2

    desc1 = property(_desc_getter(1), _desc_setter(1))
    desc2 = property(_desc_getter(2), _desc_setter(2))

    def compdesc(self, max_res=0, debug=False):
        """
        Computes alignments of the desc2 to the desc1.

        Parameters:
            * max_res -- max number of results to be returned (0 means any)
            * debug -- if set an instance of t_compdesc_debug is returned along with the resuls.

        Results are sorted by size (decreasing) and RMSD (increasing).
        """

        if self._desc1 is None:
            raise Exception("desc1 not set.")

        if self._desc2 is None:
            raise Exception("desc2 not set.")

        if debug and not libcompdesc.has_debug():
            raise Exception("Compdesc is compiled without debug option.")

        pars = t_compdesc_parameters()
        results_p = ctypes.POINTER(t_compdesc_result)()

        if debug:
            debug_p = ctypes.POINTER(t_compdesc_debug)()
        else:
            debug_p = None

        n_results = libcompdesc.compdesc(
            self._c_desc1, self._c_desc2, max_res, pars, results_p, debug_p
        )

        results = [res.unpack(self.desc1, self.desc2) for res in results_p[0:n_results]]

        if debug:
            return (
                results,
                debug_p.contents,
                results_p[0:n_results],
                self._c_desc1,
                self._c_desc2,
            )
        else:
            libcompdesc.free_results(results_p, n_results)
            return results


def compdesc(desc1, desc2, max_res=0, debug=False):
    """
        Computes alignments of the desc2 to the desc1.

        Parameters:
            * max_res -- max number of results to be returned (0 means any)
            * debug -- if set an instance of t_compdesc_debug is returned along with the resuls.

        Results are sorted by size (decreasing) and RMSD (increasing).
    """
    compdesc_obj = CompDesc(desc1, desc2)
    return compdesc_obj.compdesc(max_res, debug)
