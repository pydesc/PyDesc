# Copyright 2017 Pawel Daniluk, Maciek Dziubinski
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
Interface to FitDesc.


created: 03.12.2013 - Maciek Dziubinski
overhauled: 7.04.2014 - Pawel Daniluk
"""

import ctypes
from ctypes import POINTER

import pydesc.alignment as alignment
import pydesc.cydesc as cydesc
import pydesc.cydesc.overfit as overfit
from pydesc.alignment.factory import AlignmentFactory
from pydesc.cydesc.ctopology import CStructure

libfitdesc = cydesc.load_library("fitdesc")


class t_fitdesc_result(ctypes.Structure):
    """
    Class corresponging to t_fitdesc_result in fitdesc.h.

    Contains a single fitdesc assignment.

    Probably this class should be merged with pydesc.cydesc.desc_comp.t_desc_comp_result.
    """

    _fields_ = [
        ("n_monomers", ctypes.c_int),
        ("motif_monomers", POINTER(ctypes.c_int)),
        ("structure_monomers", POINTER(ctypes.c_int)),
        ("RMSD", ctypes.c_float),
        ("TR", overfit.t_transrot),
    ]

    def unpack(self, motif, structure):
        """
        Returns a triple containing RMSD value, PairAlignment object and a TRTMatrix.

        motif and structure are AbstractStructure instances which the assignment refers to.
        """
        n_mers = self.n_monomers
        get_mers = lambda struct, array: list(map(struct.__getitem__, array[0:n_mers]))
        alignment_factory = AlignmentFactory()
        structures = motif.derived_from, structure.derived_from
        inds_lists = [
            self.motif_monomers[0 : self.n_monomers],
            self.structure_monomers[0 : self.n_monomers],
        ]
        pair_al = alignment_factory.create_from_list_of_inds(structures, inds_lists)
        return self.RMSD, pair_al, self.TR.to_trtmatrix()


libfitdesc.fitdesc.argtypes = [
    ctypes.POINTER(CStructure),
    ctypes.POINTER(CStructure),
    ctypes.c_float,
    ctypes.c_int,
    ctypes.POINTER(ctypes.POINTER(t_fitdesc_result)),
]
libfitdesc.free_results.argtypes = [ctypes.POINTER(t_fitdesc_result), ctypes.c_int]


class FitDesc(object):

    """
    Computation of assignment using the FitDesc algorithm.


    Allows storing a structure and motif in order to avoid overhead of creating
    corresponding structures in C many times.
    """

    def __init__(self, motif=None, structure=None):
        """
        Initializator. Can set motif and structure attributes.
        """
        self.motif = motif
        self.structure = structure

    @property
    def structure(self):
        """
        A structure the motif is being assigned to.
        """
        return self._structure

    @structure.setter
    def structure(self, struct):
        """
        Structure setter. Creates a corresponding C structure.
        """

        # This method is called in __init__
        # pylint: disable=W0201
        self._structure = struct
        self._c_structure = None if struct is None else CStructure(struct)

    @property
    def motif(self):
        """
        A motif being assigned to the structure.
        """
        return self._motif

    @motif.setter
    def motif(self, struct):
        """
        Motif setter. Creates a corresponding C structure.
        """

        # This method is called in __init__
        # pylint: disable=W0201
        self._motif = struct
        self._c_motif = None if struct is None else CStructure(struct)

    def fitdesc(self, max_rmsd, max_res=0):
        """
        Computes FitDesc assignments of the motif to the structure.

        Parameters:
            * r -- maximum allowed RMSD
            * max_res -- max number of results to be returned (0 means any)

        Results are sorted by increasing RMSD.
        """

        if self._structure is None or self._motif is None:
            raise Exception("Structure not set.")

        if self._motif is None:
            raise Exception("Motif not set.")

        results_p = ctypes.POINTER(t_fitdesc_result)()

        n_results = libfitdesc.fitdesc(
            self._c_motif, self._c_structure, max_rmsd, max_res, results_p
        )

        results = [
            res.unpack(self.motif, self.structure) for res in results_p[0:n_results]
        ]

        libfitdesc.free_results(results_p, n_results)

        return results
