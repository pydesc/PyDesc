# Copyright 2017 Maciej Dziubinski
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

"""Functions and interfaces enabling running alignment methods written in C.

created: 25.12.2013 - Maciej Dziubinski

"""

import ctypes
import ctypes.util
import os

import pkg_resources

import pydesc


def load_library(name):
    """
    Finds and loads a dynamic library. Returns an instance of ctypes.CDLL.


    If pydesc has been imported from a package it uses pkg_resources module to unpack and retrieve
    libraries. Otherwise it searches lib directory in the directory pydesc package is located.

    To provide portability it searches for files with extensions .dll, .so and .dylib.
    """

    extensions = [".dll", ".so", ".dylib"]

    try:
        req = pkg_resources.get_distribution("pydesc").as_requirement()
        if pydesc.__file__ != pkg_resources.resource_filename(req, "pydesc/maps.py"):
            req = None

    except pkg_resources.DistributionNotFound:
        req = None

    for ext in extensions:
        try:
            if req is not None:
                file_name = pkg_resources.resource_filename(
                    req, "/lib/lib%s%s" % (name, ext)
                )
            else:
                file_name = os.path.join(
                    os.path.dirname(os.path.dirname(pydesc.__file__)),
                    "lib/lib%s%s" % (name, ext),
                )
        except KeyError:
            continue

        if os.path.isfile(file_name):
            return ctypes.CDLL(file_name)

    raise Exception("Could not load lib%s." % name)


CYDESC_LIB = load_library("cydesc")


def use_library(lib):
    """Decorator for CInDelMeta metaclass setting a C library to
    be searched for functions."""

    def class_wrapper(mcs):
        """Binds CInDelMeta to %s."""

        class WrappingClass(mcs):
            """This class wraps CInDelMeta to change library object."""

            _lib = lib

        return WrappingClass

    class_wrapper.__doc__ = class_wrapper.__doc__ % lib._name
    return class_wrapper
