from pydesc.alignment.base import DASH
from pydesc.api.selection import get_selection_from_sub_structure


def get_selections(alignment):
    dct = get_partial_structures(alignment)
    for structure, sub_structure in dct.items():
        selection = get_selection_from_sub_structure(sub_structure)
        dct[structure] = selection
    return dct


def get_partial_structures(alignment):
    structures = alignment.structures
    result = {}
    for col, structure in zip(alignment.iter_columns(), structures):
        inds = col[col != DASH]
        partial_structure = structure[inds]
        result[structure] = partial_structure
    return result
