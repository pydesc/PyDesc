import numpy
import pytest

from pydesc.alignment.base import DASH
from pydesc.alignment.loaders import PALLoader
from pydesc.api.alignment import get_partial_structures
from pydesc.api.alignment import get_selections
from pydesc.api.structure import get_structures_from_file


@pytest.fixture(scope="session")
def kinases(alignments_dir):
    stc_path = alignments_dir / "structures" / "kinases3_dama.pdb"
    structures = get_structures_from_file(stc_path)
    return structures


@pytest.fixture(scope="session")
def kinases_alignment(alignments_dir, kinases):
    path = alignments_dir / "pal" / "kinases3_dama.pal"
    with open(path) as fh:
        loader = PALLoader(fh)
        alignment = loader.load_alignment(kinases)
    alignment = alignment.close()
    return alignment


def test_selections_sars(alignments_dir):
    stc_path = alignments_dir / "structures" / "sars_pair.pdb"
    structures = get_structures_from_file(stc_path)
    path = alignments_dir / "pal" / "sars_pair.pal"

    with open(path) as file_:
        loader = PALLoader(file_)

    alignment = loader.load_alignment(structures)
    selections_dct = get_selections(alignment)

    for structure, expected_inds in zip(alignment.structures, alignment.iter_columns()):
        selection = selections_dct[structure]
        inds = selection.get_list_of_inds(structure)
        numpy.testing.assert_array_equal(expected_inds, inds)


def test_selections_3_kinases(kinases_alignment):
    selections_dct = get_selections(kinases_alignment)
    for structure, expected_inds in zip(
        kinases_alignment.structures, kinases_alignment.iter_columns()
    ):
        selection = selections_dct[structure]
        inds = selection.get_list_of_inds(structure)
        expected_inds = expected_inds[numpy.where(expected_inds != DASH)]
        numpy.testing.assert_array_equal(expected_inds, inds)


def test_get_partial_structures(kinases_alignment):
    dct = get_partial_structures(kinases_alignment)
    col_inds_map = kinases_alignment.get_structure_indices()

    for structure, sub_structure in dct.items():
        assert structure is sub_structure.derived_from
        col_ind = col_inds_map[structure]
        col = kinases_alignment.inds[:, col_ind]
        col_len = numpy.count_nonzero(col[col != DASH])
        assert len(sub_structure) == col_len
        assert len(sub_structure) > 46

    # order of columns is different
    lens = {len(item) for item in dct.values()}
    assert lens == {238, 47}
