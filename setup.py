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

# !/usr/bin/env python

"""
PyDesc setup script

created: 27.03.2014 - Pawel Daniluk

To install PyDesc type
    python ./setup.py install

"""

import os
import os.path
from distutils.core import setup

try:
    import bitten.util.testrunner

    test = True
except ImportError:
    test = False

from shlib.build_shlib import SharedLibrary
from shlib.build_shlib import build_shlib
from shlib.install_shlib import develop
from shlib.install_shlib import install_lib

cmdclass = {"install_lib": install_lib, "build_shlib": build_shlib, "develop": develop}

if test:
    cmdclass["unittest"] = bitten.util.testrunner.unittest


def discover_libdesc(path):
    """ Searches for libdesc includes in given path. Returns a boolean value
    and a path or an empty string.
    """

    inc_path = os.path.join(path, "include")

    inc_files = ["desc.h", "compdesc_debug.h"]

    if all([os.path.isfile(os.path.join(inc_path, f)) for f in inc_files]):
        return True, path
    else:
        return False, ""


has_libdesc = False

try:
    has_libdesc, libdesc_path = discover_libdesc(os.environ["LIBDESC_PATH"])
except KeyError:
    locations = [os.path.join(os.environ["HOME"], "local"), "/usr/local", "/usr"]
    for loc in locations:
        has_libdesc, libdesc_path = discover_libdesc(loc)
        if has_libdesc:
            break

compile_args = [
    "-Wall",
    "-Wno-unknown-pragmas",
    "-std=gnu99",
    "-mmmx",
    "-msse",
    "-msse3",
    "-D_GNU_SOURCE",
    "-m64",
    "-fopenmp",
]

compdesc_sources = [
    "compdesc.c",
    "rmsd.c",
    "alignment.c",
    "01_prelim.c",
    "02_components.c",
    "03_build.c",
    "04_prune.c",
]
prepend_with_dir = lambda x: "src/compdesc/" + x
compdesc_sources = list(map(prepend_with_dir, compdesc_sources))

shlibs = [
    SharedLibrary(
        name="cydesc",
        sources=[
            "src/cstructures.c",
            "src/util/arrays.c",
            "src/util/bitops.c",
            "src/util/maps.c",
            "src/util/random.c",
        ],
        include_dirs=["src/include"],
        libraries=["m"],
        extra_compile_args=list(compile_args),
    ),
    SharedLibrary(
        name="overfit",
        sources=["src/overfit.c"],
        include_dirs=["src/include"],
        libraries=["m"],
        extra_compile_args=list(compile_args),
    ),
    SharedLibrary(
        name="fitdesc",
        sources=["src/fitdesc.c", "src/overfit.c"],
        include_dirs=["src/include"],
        libraries=["m"],
        extra_compile_args=list(compile_args),
    ),
    SharedLibrary(
        name="compdesc",
        sources=compdesc_sources
        + ["src/overfit.c", "src/util/maps.c", "src/util/arrays.c"],
        include_dirs=["src/include"],
        libraries=["m"],
        extra_compile_args=list(compile_args),
    ),
    SharedLibrary(
        name="cydesc_test",
        sources=["src/cstructures_test.c", "src/util/maps.c"],
        include_dirs=["src/include"],
        libraries=["m"],
        extra_compile_args=list(compile_args),
    ),
]

if has_libdesc:
    libdesc_inc_path = os.path.join(libdesc_path, "include")
    libdesc_lib_path = os.path.join(libdesc_path, "lib")
    shlibs.append(
        SharedLibrary(
            name="compdesc_test",
            sources=["src/compdesc_test.c", "src/util/arrays.c"],
            include_dirs=["src/include", libdesc_inc_path],
            libraries=[
                "m",
                "desc",
                "gen_avl",
                "overfit",
                "desc_comp_domains",
                "arrays",
                "png",
                "graph",
                "motzkin",
            ],
            extra_compile_args=compile_args,
            library_dirs=[libdesc_lib_path],
        )
    )

    for lib in shlibs:
        if lib.name == "compdesc":
            # lib.extra_compile_args.append('-DCOMPDESC_DEBUG')
            lib.include_dirs.append(libdesc_inc_path)

setup(
    name="PyDesc",
    version="0.x",
    description="Toolkit for analysis of biopolymer structures using"
    " local descriptors",
    author="Pawel Daniluk",
    author_email="pawel@bioexploratorium.pl",
    url="http://trac.dw/trac",
    packages=["pydesc", "pydesc.cydesc"],
    cmdclass=cmdclass,
    test_suite="tests_legacy",
    shlibs=shlibs,
    install_requires=["biopython", "mdtraj",],
)
