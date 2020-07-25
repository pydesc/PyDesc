# Copyright 2020 Tymoteusz Oleniecki
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
"""Convenience functions for full atom representation of mers."""

from pydesc.chemistry.full_atom import Residue
from pydesc.chemistry.full_atom import calculate_residues_angles_vectorized
from pydesc.selection import AtomSetSubclass


def calculate_residues_angles(structure):
    """Calculate torsion angles of residues in given structure.

       Only subclasses of Residue are taken into account.

       Sets "angles" property in residues.
       Purpose of this method is to speed up calculation of angles -- its faster
       than calling property in every residue.

       Args:
           structure: structure instance (preferably with some Residues in it).

   """
    residues = AtomSetSubclass(Residue).create_structure(structure)
    calculate_residues_angles_vectorized(residues)
