# Copyright 2017 Tymoteusz Oleniecki
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
"""Package providing convenience functions for dealing with structure creation."""

from pydesc.structure import StructureLoader


def get_structures(name):
    """Load structures using default behaviour.

    Default behaviour means:
    - try to get structure from cache dir
    - if it is not stored yet, try to download it

    Args:
        name(str): proper structure identifier. See pydes.dbhandler for more
        information.

    Returns:
        : list of structures. There will be more than one structure if structure
        contains more than one model (trajectory, biounit or NMR).

    """
    sl = StructureLoader()
    res = sl.load_structures(name)
    return res


def get_structures_from_file(path):
    """Load file from local path.

    Args:
        path(str): relative or absolute path to file (pdb or cif).

    Returns:
        : list of structures. There will be more than one structure if structure
        contains more than one model (trajectory, biounit or NMR).

    """
    sl = StructureLoader()
    res = sl.load_structures(path=path)
    return res
