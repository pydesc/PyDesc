from pydesc.selection import Set


def get_selection_from_sub_structure(sub_structure):
    converter = sub_structure.derived_from.converter
    pdb_ids = [converter.get_pdb_id(mer.ind) for mer in sub_structure]
    selection = Set(pdb_ids)
    return selection
