import numpy
import pytest

from pydesc.alignment.base import DASH
from pydesc.alignment.loaders import CSVLoader
from pydesc.alignment.loaders import FASTALoader
from pydesc.alignment.loaders import PALLoader
from pydesc.api.alignment import get_loader
from pydesc.api.alignment import get_partial_structures
from pydesc.api.alignment import get_selections
from pydesc.api.alignment import load_alignment
from pydesc.api.structure import get_structures_from_file


@pytest.fixture(scope="session")
def kinases_al_pth(alignments_dir):
    path = alignments_dir / "pal" / "kinases3_dama.pal"
    return path


@pytest.fixture(scope="session")
def kinases(alignments_dir):
    stc_path = alignments_dir / "structures" / "kinases3_dama.pdb"
    structures = get_structures_from_file(stc_path)
    return structures


@pytest.fixture(scope="session")
def kinases_alignment(alignments_dir, kinases, kinases_al_pth):
    with open(kinases_al_pth) as fh:
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


@pytest.mark.parametrize(
    "path,loader_type",
    [
        ("a.pal", PALLoader),
        ("a.fasta", FASTALoader),
        ("a.bla", FASTALoader),
        ("a.csv", CSVLoader),
        ("a.tsv", CSVLoader),
    ],
)
def test_get_loader(path, loader_type, tmp_path):
    pth = tmp_path / path
    pth.write_text("10\ntest")
    loader = get_loader(str(pth))
    assert type(loader) == loader_type


def test_load_alignment(kinases, kinases_al_pth, kinases_alignment):
    alignment = load_alignment(kinases_al_pth, kinases).close()
    numpy.testing.assert_array_equal(alignment.inds, kinases_alignment.inds)
    assert alignment.structures == kinases_alignment.structures
