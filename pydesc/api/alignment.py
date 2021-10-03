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
"""Functions facilitating usage of alignments."""

from pathlib import Path

from pydesc.alignment.base import DASH
from pydesc.alignment.loaders import CSVLoader
from pydesc.alignment.loaders import FASTALoader
from pydesc.alignment.loaders import PALLoader
from pydesc.api.selection import get_selection_from_sub_structure


def get_loader(path):
    """Get loader for alignment in file with given path.

    Args:
        path(str, Path): path to alignment file (fasta, tsv/csv or pal).

    Returns:
        alignment loader bonded with given path.

    """
    path = Path(str(path))
    extension = path.suffix
    class_dct = {
        ".pal": PALLoader,
        ".tsv": CSVLoader,
        ".csv": CSVLoader,
        ".fasta": FASTALoader,
    }
    klass = class_dct.get(extension, FASTALoader)
    with open(path) as file_handler:
        loader = klass(file_handler)
    return loader


def load_alignment(path, structures):
    """Load alignment from file in given path with given structures.

    Args:
        path(str, Path): path to alignment file (fasta, tsv/csv or pal).
        structures: list of structures or dict mapping labels to structures.

    Returns:
        alignment (multiple or pair).

    """
    loader = get_loader(path)
    if isinstance(structures, dict):
        alignment = loader.load_alignment_mapping(structures)
    else:
        alignment = loader.load_alignment(structures)
    return alignment


def get_selections(alignment):
    """Get sets of PDB ids of aligned mers.

    Args:
        alignment: instance of pydesc.alignments.base.AbstractAlignment subclass
            representing aligned structures.

    Returns:
        dict: with structures being keys and selections (Sets) being values.

    """
    dct = get_partial_structures(alignment)
    for structure, sub_structure in dct.items():
        selection = get_selection_from_sub_structure(sub_structure)
        dct[structure] = selection
    return dct


def get_partial_structures(alignment):
    """Get partial structures containing only aligned mers.

    Args:
        alignment: instance of pydesc.alignments.base.AbstractAlignment subclass
            representing aligned structures.

    Returns:
        dict: with structures being keys and partial structures being values.

    """
    structures = alignment.structures
    result = {}
    for col, structure in zip(alignment.iter_columns(), structures):
        inds = col[col != DASH]
        partial_structure = structure[inds]
        result[structure] = partial_structure
    return result
