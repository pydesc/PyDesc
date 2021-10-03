# Copyright 2017 Pawel Daniluk, Grzegorz Firlik, Tymoteusz Oleniecki
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

"""PDB and mmCIF file parsers."""

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser


class MetaParser:
    """Wrapper and composition of PDB- and mmCIF-parsers."""

    def __init__(self, *args, **kwargs):
        self.parsers = [PDBParser(*args, **kwargs), MMCIFParser(*args, **kwargs)]

    def get_structure(self, stc, file, *args, **kwargs):
        """Wrapper of get_structure method of both parsers."""
        for parser in self.parsers:
            try:
                return parser.get_structure(stc, file, *args, **kwargs)
            except ValueError:  # the only exception known to be raised when
                # proper mmCIF is passed to PDBParser
                file.seek(0)
                continue
        raise ValueError(
            "None of parsers could get %s structure. Tried: %s."
            % (stc, ", ".join([type(i).__name__ for i in self.parsers]))
        )
