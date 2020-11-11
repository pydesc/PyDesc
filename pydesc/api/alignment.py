from pydesc.selection import Set
from pydesc.alignment.base import DASH


def get_selections(alignment):
    selections = {}
    for col, structure in enumerate(alignment.structures):
        converter = structure.converter
        inds = alignment.inds[:, col].flatten()
        pdb_ids = [converter.get_pdb_id(ind) for ind in inds if ind is not DASH]
        selection = Set(pdb_ids)
        selections[structure] = selection
    return selections
