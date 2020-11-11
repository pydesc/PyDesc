import numpy
from pydesc.alignment.base import DASH
from pydesc.api.alignment import get_selections
from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file


def test_selections_sars(alignments_dir):
    stc_path = alignments_dir / "structures" / "sars_pair.pdb"
    structures = get_structures_from_file(stc_path)
    path = alignments_dir / "pal" / "sars_pair.pal"
    loader = PALLoader(path)

    alignment = loader.load_alignment(structures)

    selections_dct = get_selections(alignment)

    for structure, expected_inds in zip(alignment.structures, alignment.iter_columns()):
        selection = selections_dct[structure]
        inds = selection.get_list_of_inds(structure)
        numpy.testing.assert_array_equal(expected_inds, inds)


def test_selections_3_kinases(alignments_dir):
    stc_path = alignments_dir / "structures" / "kinases3_dama.pdb"
    structures = get_structures_from_file(stc_path)
    path = alignments_dir / "pal" / "kinases3_dama.pal"
    loader = PALLoader(path)

    alignment = loader.load_alignment(structures)

    selections_dct = get_selections(alignment)

    for structure, expected_inds in zip(alignment.structures, alignment.iter_columns()):
        selection = selections_dct[structure]
        inds = selection.get_list_of_inds(structure)
        expected_inds = expected_inds[numpy.where(expected_inds != DASH)]
        numpy.testing.assert_array_equal(expected_inds, inds)
