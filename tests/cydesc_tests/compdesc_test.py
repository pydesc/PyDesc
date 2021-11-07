import pytest

from pydesc.api.cmaps import calculate_contact_map
from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.criteria import get_gc_distance_criterion
from pydesc.api.descriptor import create_descriptor
from pydesc.api.descriptor import create_descriptors
from pydesc.api.structure import get_structures
from pydesc.api.structure import get_structures_from_file
from pydesc.chemistry.bbtrace import CATrace
from pydesc.chemistry.factories import BioPythonAtomSetFactory
from pydesc.cydesc.compdesc import CompDesc
from pydesc.structure import StructureLoader


@pytest.fixture(scope="session")
def criterion():
    return get_default_protein_criterion()


def test_same_desc_gold_standard_comparison(structures_dir, criterion):
    (stc,) = get_structures_from_file(structures_dir / "prots_only" / "3NPU.pdb")
    contact_map = calculate_contact_map(stc, criterion)

    desc1 = create_descriptor(stc, stc[42], contact_map)
    desc2 = create_descriptor(stc, stc[42], contact_map)

    fitter = CompDesc(desc1, desc2)
    res = fitter.compdesc()
    assert len(res) == 3
    best_res = res[0]
    best_rmsd, best_al, best_trt = best_res
    assert len(best_al) == 30
    for ind_stc1, ind_stc2 in best_al.iter_rows():
        assert ind_stc1 == ind_stc2

    res_rmsds = [rmsd for rmsd, _, _ in res]
    # gold standard comes from last seen working version of PyDesc in Python 2.7
    gold_standard_rmsds = (8.70581544631932e-08, 1.9098607301712036, 1.9098607301712036)
    for rmsd, expected_rmsd in zip(res_rmsds, gold_standard_rmsds):
        assert rmsd == expected_rmsd


def test_very_similar_desc_comp(structures_dir, criterion):
    (stc,) = get_structures_from_file(structures_dir / "prots_only" / "3NPU.pdb")
    contact_map = calculate_contact_map(stc, criterion)

    desc1 = create_descriptor(stc, stc[42], contact_map)
    desc2 = create_descriptor(stc, stc[291], contact_map)

    fitter = CompDesc(desc1, desc2)
    res = fitter.compdesc()

    res_rmsds, res_als, res_trts = zip(*res)
    # gold standard comes from last seen working version of PyDesc in Python 2.7
    gold_standard_rmsds = (0.19330592453479767, 1.8893914222717285, 1.9629560708999634)
    for rmsd, expected_rmsd in zip(res_rmsds, gold_standard_rmsds):
        assert rmsd == expected_rmsd

    gold_standard_al_lens = (30, 27, 27)
    for al, expected_len in zip(res_als, gold_standard_al_lens):
        assert len(al) == expected_len

    best_al, _, _ = res_als
    best_al_rows = [tuple(row) for row in best_al.iter_rows()]
    assert (42, 291) in best_al_rows


def test_similar_different_structures(structures_dir, criterion):
    (stc1,) = get_structures_from_file(structures_dir / "ligands" / "1luf.cif")
    (stc2,) = get_structures_from_file(structures_dir / "ligands" / "2src.cif")
    cm1 = calculate_contact_map(stc1, criterion)
    cm2 = calculate_contact_map(stc2, criterion)

    desc1 = create_descriptor(stc1, stc1.pdb_ids["A:607"], cm1)
    desc2 = create_descriptor(stc2, stc2.pdb_ids["A:294"], cm2)

    fitter = CompDesc(desc1, desc2)
    res = fitter.compdesc()
    assert len(res) == 1
    rmsd, al, _ = res[0]
    assert pytest.approx(rmsd, 0.001) == 1.154
    assert len(al) == 27
    assert al.get_inds_aligned_with(stc1, [47]) == {stc2: [210]}


def test_CATrace_compdesc(structures_dir):
    mer_factory = BioPythonAtomSetFactory(classes=[CATrace])
    loader = StructureLoader(atom_set_factory=mer_factory)
    with open(structures_dir / "PorCA_only" / "1KAN.pdb") as fh:
        stc, = loader.load_structures([fh])
    criterion = get_gc_distance_criterion()
    cm = calculate_contact_map(stc, criterion)

    descs = tuple(filter(bool, create_descriptors(stc, cm)))
    desc1 = create_descriptor(stc, stc[19], cm)
    desc2 = create_descriptor(stc, stc[134], cm)

    fitter = CompDesc(desc1, desc2)
    res = fitter.compdesc()
    best_res, = res
    best_rmsd, _, _ = best_res
    assert best_rmsd == 1.007647156715393
    # confirmed by visual inspection
