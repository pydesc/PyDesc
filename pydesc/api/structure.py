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

from pathlib import Path

from pydesc.dbhandler import MetaHandler
from pydesc.structure import StructureLoader


def get_structures(name, common_converter=False):
    """Load structures using default behaviour.

    Default behaviour means:
    - try to get structure from cache dir
    - if it is not stored yet, try to download it

    Args:
        name(str): proper structure identifier. See pydes.dbhandler for more
            information.
        common_converter(bool): determines if all structures should share converter.
            False by default.

    Returns:
        : list of structures. There will be more than one structure if structure
        contains more than one model (trajectory, bio-unit or NMR).

    """
    sl = StructureLoader()
    with MetaHandler().open(name) as files:
        structures = sl.load_structures(files, common_converter=common_converter)
    return structures


def get_structures_from_file(path, common_converter=False):
    """Load file from local path.

    Args:
        path(str): relative or absolute path to file (pdb or cif).
        common_converter(bool): determines if all structures should share converter.
            False by default.

    Returns:
        : list of structures. There will be more than one structure if structure
        contains more than one model (trajectory, bio-unit or NMR).

    """
    path = Path(path)
    sl = StructureLoader()
    if not path.is_file():
        raise ValueError("Given path does not lead to a file.")
    with open(path) as file_:
        structures = sl.load_structures([file_], common_converter=common_converter)
    return structures
