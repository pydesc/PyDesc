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

# This is not a constant.
libfitdesc = cydesc.load_library("fitdesc")  # pylint: disable=C0103


class t_fitdesc_result(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.

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
        get_monomers = lambda struct, array: list(
            map(struct.__getitem__, array[0:n_mers])
        )
        pair_al = alignment.PairAlignment(
            [motif, structure],
            list(
                zip(
                    get_monomers(motif, self.motif_monomers),
                    get_monomers(structure, self.structure_monomers),
                )
            ),
        )
        return self.RMSD, pair_al, self.TR.to_trtmatrix()


libfitdesc.fitdesc.argtypes = [
    ctypes.POINTER(cydesc.CStructure),
    ctypes.POINTER(cydesc.CStructure),
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
        self._c_structure = None if struct is None else cydesc.CStructure(struct)

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
        self._c_motif = None if struct is None else cydesc.CStructure(struct)

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


def fitdesc(motif, struct, max_rmsd, max_res=0):
    """
        Computes FitDesc assignments of the motif to the structure.

        Parameters:
            * motif -- structure being assigned
            * struct -- structrure being assigned to
            * r -- maximum allowed RMSD
            * max_res -- max number of results to be returned (0 means any)

        Results are sorted by increasing RMSD.
    """
    fitdesc_obj = FitDesc(motif, struct)
    return fitdesc_obj.fitdesc(max_rmsd, max_res)
