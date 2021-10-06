import pytest

from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file
from pydesc.cydesc.overfit import Multifit
from pydesc.cydesc.overfit import Overfit
from pydesc.geometry import TRTMatrix


@pytest.fixture
def trivial_structures(nmr_structure_of_each_kind):
    stc1, stc2, *_ = get_structures_from_file(nmr_structure_of_each_kind)
    return stc1, stc2


@pytest.fixture
def trivial_structures3(nmr_structure_of_each_kind):
    stc1, stc2, stc3, *_ = get_structures_from_file(nmr_structure_of_each_kind)
    return stc1, stc2, stc3


@pytest.fixture
def sars_pair(alignments_dir):
    stc_path = alignments_dir / "structures" / "sars_pair.pdb"
    structures = get_structures_from_file(stc_path)
    path = alignments_dir / "pal" / "sars_pair.pal"
    with open(path) as fh:
        loader = PALLoader(fh)
    alignment = loader.load_alignment(structures)
    return alignment


class TestOverfit:
    def test_trivial(self, trivial_structures):
        fitter = Overfit()
        fitter.add_structures(*trivial_structures)
        rmsd, trt_mtx = fitter.overfit()

        assert 8.0 > rmsd >= 0.0
        assert isinstance(trt_mtx, TRTMatrix)

    def test_mers(self, trivial_structures):
        fitter = Overfit()
        for ind in range(len(trivial_structures[0]))[:3]:
            mer1, mer2 = [stc[ind] for stc in trivial_structures]
            fitter.add_mers(mer1, mer2)
        rmsd, trt_mtx = fitter.overfit()

        assert 8.0 > rmsd >= 0.0
        assert isinstance(trt_mtx, TRTMatrix)

    def test_points(self, trivial_structures):
        fitter = Overfit()
        for ind in range(len(trivial_structures[0]))[:3]:
            at1, at2 = [tuple(stc[ind])[0] for stc in trivial_structures]
            fitter.add_points(at1, at2)
        rmsd, trt_mtx = fitter.overfit()

        assert 8.0 > rmsd >= 0.0
        assert isinstance(trt_mtx, TRTMatrix)

    def test_points_different_lens(self,):
        fitter = Overfit()
        with pytest.raises(TypeError):
            fitter.add_points([0, 0, 0], [0, 0, 0, 0, 0])

    def test_artificial_points(self,):
        fitter = Overfit()
        fitter.add_points((1, 2, 3), (1, 2, 3))
        rmsd, trt = fitter.overfit()
        assert rmsd == 0.0

    def test_alignment(self, sars_pair):
        fitter = Overfit()
        fitter.add_alignment(sars_pair)
        rmsd, trt = fitter.overfit()
        assert rmsd >= 0.0

    def test_different_representation(self, protein_file, dna_file):
        (stc1,) = get_structures_from_file(protein_file)
        (stc2,) = get_structures_from_file(dna_file)
        fitter = Overfit()
        with pytest.raises(TypeError):
            fitter.add_mers(stc1[0], stc2[0])

    def test_different_structure_lens(self, protein_file, dna_file):
        (stc1,) = get_structures_from_file(protein_file)
        (stc2,) = get_structures_from_file(dna_file)
        fitter = Overfit()
        with pytest.raises(TypeError):
            fitter.add_structures(stc1, stc2)

    def test_sums(self):
        fitter = Overfit()
        sums = fitter.get_sums()
        fitter.add_sums(sums)
        sums_sq = sums + sums
        assert sums_sq
        # TODO: not sure what to test here


class TestMultifit:
    def test_points(self, trivial_structures3):
        fitter = Multifit(3)
        p1, p2, p3 = [stc[0].representation[0] for stc in trivial_structures3]
        fitter.add_points(p1, p2, p3)
        result = fitter.multifit_once()
        rmsds, matrices = zip(*result)
        for i in rmsds:
            assert i == 0.0
        for matrix in matrices:
            assert isinstance(matrix, TRTMatrix)

    def test_points_wrong_number(self):
        fitter = Multifit(4)
        with pytest.raises(TypeError):
            fitter.add_points([], [], [])

    def test_mers(self, trivial_structures3):
        fitter = Multifit(3)
        mers = [stc[0] for stc in trivial_structures3]
        fitter.add_mers(*mers)
        result = fitter.multifit_once()
        rmsds, matrices = zip(*result)
        for i in rmsds:
            assert i >= 0.0
        for matrix in matrices:
            assert isinstance(matrix, TRTMatrix)

    def test_structures(self, trivial_structures3):
        fitter = Multifit(3)
        fitter.add_structures(*trivial_structures3)
        result = fitter.multifit_once()
        rmsds, matrices = zip(*result)
        for i in rmsds:
            assert 0.0 < i < 8.0
        for matrix in matrices:
            assert isinstance(matrix, TRTMatrix)

    def test_alignment(self, sars_pair):
        fitter = Multifit(2)
        fitter.add_alignment(sars_pair)
        result = fitter.multifit()
        rmsds, matrices = zip(*result)
        for i in rmsds:
            assert 0.0 < i < 8.0
        for matrix in matrices:
            assert isinstance(matrix, TRTMatrix)
