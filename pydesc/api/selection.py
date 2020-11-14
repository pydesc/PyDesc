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
"""Functions facilitating usage of selections."""

from pydesc.selection import Set


def get_selection_from_sub_structure(substructure):
    """Return selection of PDB ids of mers from given substructure."""
    converter = substructure.derived_from.converter
    pdb_ids = [converter.get_pdb_id(mer.ind) for mer in substructure]
    selection = Set(pdb_ids)
    return selection
