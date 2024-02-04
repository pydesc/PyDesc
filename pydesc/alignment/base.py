# Copyright 2020 Tymoteusz 'vdhert' Oleniecki
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
"""Basic classes storing alignments."""

from collections import defaultdict
from functools import wraps

import numpy


def none_returns_none(function):
    @wraps(function)
    def wrapper(self, arg):
        if arg is None:
            return None
        return function(self, arg)

    return wrapper


class _Dash:
    """Singleton representing missing value in alignment."""

    def __repr__(self):
        return "<->"

    def __gt__(self, other):
        return False

    def __lt__(self, other):
        return True


DASH = _Dash()


class IncorrectAlignmentError(Exception):
    pass


class Alignment:
    """Structure storing alignment between mers of some structures."""

    def __init__(self, structures, table):
        if len(structures) < 2:
            raise IncorrectAlignmentError(
                "Cannot create an alignment with only one structure."
            )
        if len(structures) != table.shape[1]:
            raise IncorrectAlignmentError(
                "Number of structures do not match table shape."
            )
        self._structures = tuple(structures)
        self._table = table
        self._mer_map = {structure: {} for structure in structures}
        self._fill_mer_map()
        self.pdb_ids = PDBidGetter(self)

    def _fill_mer_map(self):
        mer_map = {structure: defaultdict(list) for structure in self._structures}
        for no, row in enumerate(self._table):
            for structure, ind in zip(self._structures, row):
                if ind is DASH:
                    continue
                mer_map[structure][ind].append(no)
        for structure in self._structures:
            stc_mer_map = mer_map[structure]
            for ind, occurs in stc_mer_map.items():
                self._mer_map[structure][ind] = numpy.array(occurs, dtype=numpy.uint32)

    def __getitem__(self, item):
        if isinstance(item, slice):
            start_index = self._get_min_index(item.start)
            stop_index = 1 + self._get_max_index(item.stop)
            arr = self._table[start_index:stop_index, :]
        else:
            indices = self._get_indices(item)
            arr = self._table[indices, :]
        return Alignment(self._structures, arr)

    def __iter__(self):
        return iter(self._table)

    @none_returns_none
    def _get_min_index(self, structure_with_ind):
        return self._get_indices(structure_with_ind).min()

    @none_returns_none
    def _get_max_index(self, structure_with_ind):
        return self._get_indices(structure_with_ind).max()

    @none_returns_none
    def _get_indices(self, structure_with_inds):
        structure, inds = structure_with_inds
        if isinstance(inds, int):
            inds = (inds,)
        return numpy.concatenate([self._mer_map[structure][ind] for ind in inds])

    def __len__(self):
        length, _ = self._table.shape
        return length

    def get_structures(self):
        return self._structures

    def get_inds_table(self):
        return numpy.array(self._table)

    def get_inds_map(self):
        return dict(self._mer_map)


class PDBidGetter:
    def __init__(self, alignment):
        self.alignment = alignment
        self._names_map = {
            structure.name: structure for structure in alignment.get_structures()
        }

    def __getitem__(self, item):
        if isinstance(item, slice):
            start = self._str_to_structure_with_ind(item.start)
            stop = self._str_to_structure_with_ind(item.stop)
            return self.alignment[start:stop]
        structure_with_ind = self._str_to_structure_with_ind(item)
        return self.alignment[structure_with_ind]

    def _str_to_structure_with_ind(self, item):
        if item is None:
            return None
        structure_name, mer_pdb_id = item.split(":", 1)
        try:
            structure = self._names_map[structure_name]
        except KeyError:
            names = ", ".join(self._names_map.keys())
            msg = f"{structure_name} is not a valid name for this alignment. Choose from: {names}"
            raise KeyError(msg)
        return structure, structure.pdb_ids[mer_pdb_id].ind
