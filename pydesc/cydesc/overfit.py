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
from ctypes import byref

import operator

import pydesc.cydesc as cydesc
import pydesc.geometry as geometry

import numpy
from functools import reduce

# This is not a constant.
liboverfit = cydesc.load_library('overfit')  # pylint: disable=C0103

# This is not a constant.
c_float_4 = ctypes.c_float * 4  # pylint: disable=C0103
c_float_4.__doc__ = "Array of 4 floats."


class t_transrot(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.
    """
    Class corresponding to t_transrot in overfit.h.


    Transformation is to be applied as follows:
        x' U + Tr

    *' denotes transposition

    Fourth vector element is required for proper memory alignment in C and should be ignored.
    """
    _fields_ = [('U', c_float_4 * 3),
                ('Tr', c_float_4)]

    def to_trtmatrix(self):
        """ Returns geometry.TRTMatrix containing the same transformation. """
        res = geometry.TRTMatrix()

        for i in range(3):
            for j in range(3):
                res.rotation_matrix[i, j] = self.U[j][i]
            res.post_vector[i] = self.Tr[i]

        return res


class t_overfit_sums(ctypes.Structure):  # pylint: disable=C0103
    # This class has to be named as the corresponding C structure.

    """
    Class corresponging to t_overfit_sums in overfit.h.

    Contains data accumulated in overfit.
    """

    _fields_ = [('N', ctypes.c_int),
                ('startA', c_float_4),
                ('startB', c_float_4),
                ('A', c_float_4),
                ('AA', c_float_4),
                ('B', c_float_4),
                ('BB', c_float_4),
                ('AB0', c_float_4),
                ('AB1', c_float_4),
                ('AB2', c_float_4)]

    def __add__(self, sums2):
        """
        Adds two sets of accumulated sums.

        This method guarantees a correct addition (e.g. accounting for different origins).
        """
        overfit_obj = Overfit()
        overfit_obj.add_sums(self)
        overfit_obj.add_sums(sums2)
        return overfit_obj.get_sums()


liboverfit.overfit_reset.restype = None
liboverfit.overfit_release_token.restype = None

liboverfit.overfit_sumadd_str.restype = None
liboverfit.overfit_sumadd_str.argtypes = [ctypes.c_int, ctypes.POINTER(t_overfit_sums)]

liboverfit.overfit_sumsave_str.restype = None
liboverfit.overfit_sumsave_str.argtypes = [ctypes.c_int, ctypes.POINTER(t_overfit_sums)]

liboverfit.overfit_add.restype = None
liboverfit.overfit_add.argtypes = [ctypes.c_int, c_float_4, c_float_4]

liboverfit.fast_overfit.restype = ctypes.c_float
liboverfit.fast_overfit.argtypes = [ctypes.c_int, ctypes.POINTER(t_transrot)]


def _ensure_token(func):
    """
    Decorator for methods in Overfit class. Executes a decorated method
    in context manager ensuring that token in liboverfit is correctly obtained
    and released.
    """
    def wrapper(self, *args, **kwargs):
        """%s
        This method is execued in context guaranteeing that liboverfit tokens
        are properly managed.
        """
        with self.token_context():
            return func(self, *args, **kwargs)

    wrapper.__doc__ = wrapper.__doc__ % func.__doc__
    return wrapper


class Overfit(object):
    """
    Stateful computation of RMSD using Kabsch algorithm.

    Enables accumulation of pairs of points, mers, structures etc.

    _sums attribute is updated only when token is released. Use get_sums method
    to access sums regardless of state.
    """

    def __init__(self):
        """
        Initially Overfit instance accumulator is empty.

        No additional zeroing is needed. Sums (empty) can be extracted right after creation.
        """
        self._sums = t_overfit_sums()
        self._token = None

    def _get_token(self):
        """
        Grabs a liboverfit token. Adds data stored in _sums to accumulator in liboverfit.
        """
        if self._token is None:
            self._token = liboverfit.overfit_get_token()
            liboverfit.overfit_reset(self._token)
            liboverfit.overfit_sumadd_str(self._token, self._sums)
            return True
        return False

    def _release_token(self):
        """
        Stores liboverfit accumulator in _sums and releases a token.
        """
        if self._token is not None:
            liboverfit.overfit_sumsave_str(self._token, byref(self._sums))
            liboverfit.overfit_release_token(self._token)
            self._token = None

    @contextlib.contextmanager
    def token_context(self):
        """
        Context manager responsible for acquiring and releasing overfit tokens.
        """
        got_token = self._get_token()
        try:
            yield
        finally:
            if got_token:
                self._release_token()

    @_ensure_token
    def add_sums(self, sums):
        """
        Adds previously accumulated sums to overfit accumulator.
        """
        liboverfit.overfit_sumadd_str(self._token, sums)

    @_ensure_token
    def get_sums(self):
        """
        Returns accumulated sums.
        """
        liboverfit.overfit_sumsave_str(self._token, byref(self._sums))
        return self._sums

    @_ensure_token
    def add_point(self, point1, point2):
        """
        Adds a pair of points to overfit accumulator.
        """
        lpoint1, lpoint2 = list(map(list, (point1, point2)))

        if not 3 <= len(lpoint1) <= 4 or not 3 <= len(lpoint2) <= 4:
            raise TypeError('Point must have 3 (or 4) float coordinates.')

        liboverfit.overfit_add(self._token, c_float_4(*list(lpoint1[:3] + [0])), c_float_4(*list(lpoint2[:3] + [0])))

    @_ensure_token
    def add_mer(self, mer1, mer2):
        """
        Adds a pair of mers to overfit accumulator.
        """
        lmer1 = mer1.representation
        lmer2 = mer2.representation

        if len(lmer1) != len(lmer2):
            raise TypeError('Mers have to have representations of the same length.')

        for point1, point2 in zip(lmer1, lmer2):
            self.add_point(point1, point2)

    @_ensure_token
    def add_structure(self, struct1, struct2):
        """
        Adds a pair of structures to overfit accumulator.
        """
        lstruct1 = list(struct1)
        lstruct2 = list(struct2)

        if len(lstruct1) != len(lstruct2):
            raise TypeError('Structures have to have same lengths')

        for mer1, mer2 in zip(lstruct1, lstruct2):
            self.add_mer(mer1, mer2)

    @_ensure_token
    def add_alignment(self, alignment_obj):
        list1, list2 = list(zip(*alignment_obj.aligned_mers))
        self.add_structure(list1, list2)    #??? use self.add as soon as it works

    @_ensure_token
    def add(self, list1, list2):
        """
        Adds a pair of arbitrarily nested iterables to overfit accumulator.

        They should be 'homeomorphic' and have contain points at the lowest level.
        """
        if len(list1) != len(list2):
            raise TypeError('Lists have to have same lengths')

        for obj1, obj2 in zip(list1, list2):
            try:
                self.add_point(obj1, obj2)
            except TypeError:
                self.add(obj1, obj2)

    @_ensure_token
    def overfit(self):
        """
        Computes RMSD and optimal superposition of accumulated pairs of points.
        """
        trot = t_transrot()

        rmsd = liboverfit.fast_overfit(self._token, byref(trot))

        return (rmsd, trot.to_trtmatrix())


class Multifit(object):
    """Computation of RMSDs for multiple structure imposistion."""

    def __init__(self):
        self.tokens = []

    def add_point(self, *args):
        if len(set(map(len, args)) - set([3,4])) != 0:
            raise TypeError('Point must have 3 (or 4) float coordinates.')
        if self.tokens != [] and len(args) != len(self.tokens[-1]):
            raise TypeError('Wrong number of points.')
        self.tokens.append(args)

    def add_mer(self, *args):
        for points in zip(*[[a.vector for a in i.representation] for i in args]):
            self.add_point(*points)

    def add_structure(self, *args):
        for mers in map(list, args):
            self.add_mer(*mers)

    def add_alignment(self, alg):
        """Adds representation of aligned mers.

        Argument:
        alg -- pydesc.alignment.MultipleAlignment instance.

        Only mers aligned for all structures will be taken.
        """
        for i in alg.aligned_mers:
            if str in set(map(type, i)): continue
            self.add_mer(*i)

    def iter_once(self, mtxs=None):
        if not mtxs:
            mtxs = [geometry.TRTMatrix() for i in self.tokens[0]]
        n = len(self.tokens[0])
        av_stc = [reduce(operator.add, [mtx.transform(vec=coords) for coords, mtx in zip(vec, mtxs)]) / n for vec in self.tokens]
        ovfs = [Overfit() for i in self.tokens[0]]
        for row, av_atm in zip(self.tokens, av_stc):
            for atm, ovf, mtx in zip(row, ovfs, mtxs):
                ovf.add_point(av_atm, mtx.transform(vec=atm))
        return [ovf.overfit() for ovf in ovfs]

    def multifit(self):
        dev_prv = numpy.inf
        dev = 999999
        mtxs_prv = [geometry.TRTMatrix() for i in self.tokens[0]]
        while dev_prv - dev > .1:
            dev_prv = dev
            rmsds, mtxs = list(zip(*self.iter_once(mtxs=mtxs_prv)))
            dev = max([rmsd ** 2 * len(self.tokens) for rmsd in rmsds])
            mtxs_prv = [i.combine(j) for i, j in zip(mtxs_prv, mtxs)]
            #~ print "ITER %f" % (dev_prv - dev)
        return list(zip(rmsds, mtxs_prv))

def overfit(list1, list2):
    """
    Computes RMSD and optimal superposition of iterables containing points.
    See Overfit.add and Overfit.overfit for more information.
    """
    overfit_obj = Overfit()
    with overfit_obj.token_context():
        overfit_obj.add(list1, list2)
        res = overfit_obj.overfit()

    return res
