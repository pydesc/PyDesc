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
Interface to liboverfit.


created: 27.03.2014 - Pawel Daniluk
"""

import contextlib
import ctypes
import operator
from ctypes import byref
from functools import reduce
from functools import wraps

import numpy

import pydesc.geometry as geometry
from pydesc.cydesc import load_library

_liboverfit = load_library("overfit")

_c_float_4 = ctypes.c_float * 4
_c_float_4.__doc__ = "Array of 4 floats."


class t_transrot(ctypes.Structure):
    """Class corresponding to t_transrot in overfit.h.

    Transformation is to be applied as follows:
        x' U + Tr

    *' denotes transposition

    Fourth vector element is required for proper memory alignment in C
    and should be ignored.

    """

    _fields_ = [("U", _c_float_4 * 3), ("Tr", _c_float_4)]

    def to_trtmatrix(self):
        """Returns geometry.TRTMatrix containing the same transformation."""
        res = geometry.TRTMatrix()

        for i in range(3):
            for j in range(3):
                res.rotation_matrix[i, j] = self.U[j][i]
            res.post_vector[i] = self.Tr[i]

        return res


class t_overfit_sums(ctypes.Structure):
    """Class corresponding to t_overfit_sums in overfit.h containing
    data accumulated in overfit."""

    _fields_ = [
        ("N", ctypes.c_int),
        ("startA", _c_float_4),
        ("startB", _c_float_4),
        ("A", _c_float_4),
        ("AA", _c_float_4),
        ("B", _c_float_4),
        ("BB", _c_float_4),
        ("AB0", _c_float_4),
        ("AB1", _c_float_4),
        ("AB2", _c_float_4),
    ]

    def __add__(self, sums2):
        """Adds two sets of accumulated sums.

        This method guarantees a correct addition
        (e.g. accounting for different origins).

        """
        overfit_obj = Overfit()
        overfit_obj.add_sums(self)
        overfit_obj.add_sums(sums2)
        return overfit_obj.get_sums()


_liboverfit.overfit_reset.restype = None
_liboverfit.overfit_release_token.restype = None

_liboverfit.overfit_sumadd_str.restype = None
_liboverfit.overfit_sumadd_str.argtypes = [ctypes.c_int, ctypes.POINTER(t_overfit_sums)]

_liboverfit.overfit_sumsave_str.restype = None
_liboverfit.overfit_sumsave_str.argtypes = [
    ctypes.c_int,
    ctypes.POINTER(t_overfit_sums),
]

_liboverfit.overfit_add.restype = None
_liboverfit.overfit_add.argtypes = [ctypes.c_int, _c_float_4, _c_float_4]

_liboverfit.fast_overfit.restype = ctypes.c_float
_liboverfit.fast_overfit.argtypes = [ctypes.c_int, ctypes.POINTER(t_transrot)]


def _ensure_token(func):
    """Decorator for methods in Overfit class. Executes a decorated method
    in context manager ensuring that token in liboverfit is correctly obtained
    and released.
    """

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        with self.token_context():
            return func(self, *args, **kwargs)

    return wrapper


class Overfit:
    """Stateful computation of RMSD using Kabsch algorithm.

    Enables accumulation of pairs of points, mers, structures etc.

    _sums attribute is updated only when token is released. Use get_sums method
    to access sums regardless of state.

    """

    def __init__(self):
        """Initially Overfit instance accumulator is empty.

        No additional zeroing is needed.
        Sums (empty) can be extracted right after creation.

        """
        self._sums = t_overfit_sums()
        self._token = None

    def _get_token(self):
        """Grabs a liboverfit token. Adds data stored in _sums to
        accumulator in liboverfit."""
        self._token = _liboverfit.overfit_get_token()
        _liboverfit.overfit_reset(self._token)
        _liboverfit.overfit_sumadd_str(self._token, self._sums)

    def _has_token(self):
        return self._token is not None

    def _release_token(self):
        """Stores liboverfit accumulator in _sums and releases a token."""
        _liboverfit.overfit_sumsave_str(self._token, byref(self._sums))
        _liboverfit.overfit_release_token(self._token)
        self._token = None

    @contextlib.contextmanager
    def token_context(self):
        """Context manager responsible for acquiring and releasing overfit tokens."""
        if not self._has_token():
            self._get_token()
        try:
            yield
        finally:
            if self._has_token():
                self._release_token()

    @_ensure_token
    def add_sums(self, sums):
        """Adds previously accumulated sums to overfit accumulator."""
        _liboverfit.overfit_sumadd_str(self._token, sums)

    @_ensure_token
    def get_sums(self):
        """Returns accumulated sums."""
        _liboverfit.overfit_sumsave_str(self._token, byref(self._sums))
        return self._sums

    @_ensure_token
    def add_points(self, point1, point2):
        """Adds a pair of points to overfit accumulator."""
        lpoint1 = list(point1)
        lpoint2 = list(point2)

        if not 3 <= len(lpoint1) <= 4 or not 3 <= len(lpoint2) <= 4:
            raise TypeError("Point must have 3 (or 4) float coordinates.")

        _liboverfit.overfit_add(
            self._token,
            _c_float_4(*list(lpoint1[:3] + [0])),
            _c_float_4(*list(lpoint2[:3] + [0])),
        )

    @_ensure_token
    def add_mers(self, mer1, mer2):
        """Adds a pair of mers to overfit accumulator."""
        lmer1 = mer1.representation
        lmer2 = mer2.representation

        if len(lmer1) != len(lmer2):
            raise TypeError("Mers have to have representations of the same length.")

        for point1, point2 in zip(lmer1, lmer2):
            self.add_points(point1, point2)

    @_ensure_token
    def add_structures(self, struct1, struct2):
        """Adds a pair of structures to overfit accumulator."""
        lstruct1 = list(struct1)
        lstruct2 = list(struct2)

        if len(lstruct1) != len(lstruct2):
            raise TypeError("Structures have to have same lengths")

        for mer1, mer2 in zip(lstruct1, lstruct2):
            self.add_mers(mer1, mer2)

    @_ensure_token
    def add_alignment(self, alignment_obj):
        stc1, stc2 = alignment_obj.structures
        for ind1, ind2 in alignment_obj.iter_rows():
            self.add_mers(stc1[ind1], stc2[ind2])

    @_ensure_token
    def overfit(self):
        """Computes RMSD and optimal superposition of accumulated pairs of points."""
        trot = t_transrot()
        rmsd = _liboverfit.fast_overfit(self._token, byref(trot))
        return rmsd, trot.to_trtmatrix()


class Multifit:
    """Computation of RMSDs for multiple structure imposition."""

    def __init__(self, n_structures):
        self.points = []
        self.n = n_structures

    def add_points(self, *points):
        if len(points) != self.n:
            msg = f"This multifit expects {self.n} points, got {len(points)}."
            raise TypeError(msg)
        points = tuple(map(tuple, points))
        self.points.append(points)

    def add_mers(self, *mers):
        for points in zip(*[[a.vector for a in i.representation] for i in mers]):
            self.add_points(*points)

    def add_structures(self, *structures):
        for mers in zip(*map(list, structures)):
            self.add_mers(*mers)

    def add_alignment(self, alignment):
        structures = alignment.structures
        for row in alignment.iter_rows():
            mers = [stc[ind] for stc, ind in zip(structures, row)]
            self.add_mers(*mers)

    def iter_once(self, matrices=None):
        if matrices is None:
            matrices = [geometry.TRTMatrix() for _ in range(self.n)]
        average_structure = []
        for points_row in self.points:
            points_row = [
                matrix.transform(vec=point)
                for point, matrix in zip(points_row, matrices)
            ]
            average_point = reduce(operator.add, points_row) / self.n
            average_structure.append(average_point)
        overfits = [Overfit() for _ in range(self.n)]
        for points_row, average_point in zip(self.points, average_structure):
            for point, overfit, matrix in zip(points_row, overfits, matrices):
                overfit.add_points(average_point, matrix.transform(vec=point))
        return [ovf.overfit() for ovf in overfits]

    def multifit(self):
        dev_previous = numpy.inf
        dev = 999999
        matrices_previous = [geometry.TRTMatrix() for i in self.points[0]]
        while dev_previous - dev > 0.1:
            dev_previous = dev
            rmsds, matrices = list(zip(*self.iter_once(matrices=matrices_previous)))
            dev = max([rmsd ** 2 * len(self.points) for rmsd in rmsds])
            matrices_previous = [
                i.combine(j) for i, j in zip(matrices_previous, matrices)
            ]
        return list(zip(rmsds, matrices_previous))
