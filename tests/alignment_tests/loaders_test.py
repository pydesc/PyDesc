from unittest.mock import MagicMock

import numpy
import pytest

from pydesc.alignment.base import Alignment
from pydesc.alignment.loaders import DASH
from pydesc.alignment.loaders import CSVLoader
from pydesc.alignment.loaders import FASTALoader
from pydesc.alignment.loaders import PALLoader
from pydesc.alignment.loaders import XMLLoader
from pydesc.alignment.processors import Close
from pydesc.api.structure import get_structures_from_file


@pytest.fixture
def structure_1l0q(alignments_dir):
    stc_path = alignments_dir / "structures" / "1l0q.pdb"
    return get_structures_from_file(stc_path)[0]


class TestCSVLoader:
    @pytest.fixture(scope="session")
    def csv_path(self, alignments_dir):
        return alignments_dir / "csv"

    def test_multi(self, csv_path):
        path = csv_path / "artificial_multi.csv"
        structures = [MagicMock() for _ in range(4)]
        for mocked_structure in structures:
            mocked_structure.converter.get_ind.side_effect = [i for i in range(1, 7)]

        with open(path) as fh:
            loader = CSVLoader(fh)

        alignment = loader.load_alignment(structures)

        assert alignment.get_inds_table().shape == (6, 4)
        expected_not_nans = [5, 4, 3, 5]
        real_not_nans = numpy.count_nonzero(alignment.get_inds_table() != DASH, axis=0)
        numpy.testing.assert_equal(real_not_nans, expected_not_nans)

        assert len(alignment.get_structures()) == 4

        stc_w_icodes = structures[2]
        calls = stc_w_icodes.converter.get_ind.call_args_list
        for call in calls[1:]:
            arg = call.args[0]
            assert arg.icode

        stc_wo_icodes = structures[1]
        calls = stc_wo_icodes.converter.get_ind.call_args_list
        for call in calls[1:]:
            arg = call.args[0]
            assert arg.icode is None

        stc_w_2letter_chain = structures[3]
        calls = stc_w_2letter_chain.converter.get_ind.call_args_list
        for call in calls[1:]:
            arg = call.args[0]
            assert len(arg.chain) == 2

    def test_pair(self, csv_path):
        path = csv_path / "artificial_pair.csv"
        structures = [MagicMock() for _ in range(2)]
        with open(path) as file_:
            loader = CSVLoader(file_)
            alignment = loader.load_alignment(structures)
        assert alignment.get_inds_table().shape == (6, 2)

    def test_no_content(self, csv_path):
        path = csv_path / "empty_pair.csv"
        with open(path) as fh:
            loader = CSVLoader(fh)
        structures = [MagicMock(name=i) for i in range(2)]
        alignment = loader.load_alignment(structures)
        assert len(alignment) == 0
        assert len(alignment.get_structures()) == 2

    def test_load_chosen_structures(self, csv_path):
        path = csv_path / "artificial_multi.csv"
        with open(path) as fh:
            loader = CSVLoader(fh)
        stc_map = {k: MagicMock(name=k) for k in ("stc1", "stc3A")}
        alignment = loader.load_alignment_mapping(stc_map)
        assert len(alignment.get_structures()) == 2
        for stc in alignment.get_structures():
            assert stc in stc_map.values()

    def test_self_aligned(self, alignments_dir, structure_1l0q):
        path_self_aligned = alignments_dir / "csv" / "1l0q_self.csv"
        with open(path_self_aligned) as fh:
            loader = CSVLoader(fh)
        alignment = loader.load_alignment_mapping({"1l0q": structure_1l0q})
        assert alignment.get_inds_table().shape == (41, 7)


class TestPALLoader:
    @pytest.fixture(scope="function")
    def artificial_multi_structures(self, mocked_structures3):
        mocked_chain = MagicMock(chain_name="K")
        mocked_structures3[0].chains = [mocked_chain]
        mocked_mers = [MagicMock(ind=i) for i in range(5)]
        # there should be 5 mers in each structure
        for structure in mocked_structures3:
            structure.converter.get_ind.return_value = 0
            structure.__getitem__.return_value = mocked_mers
        return mocked_structures3

    @pytest.fixture(scope="session")
    def artificial_multi_path(self, alignments_dir):
        return alignments_dir / "pal" / "artificial_multi.pal"

    def test_artificial_multi(self, artificial_multi_path, artificial_multi_structures):
        structures = artificial_multi_structures
        with open(artificial_multi_path) as fh:
            loader = PALLoader(fh)
        metadata = loader.read_metadata()
        for mol in "MOL1", "MOL2", "MOL3":
            assert mol in metadata["labels"]
        joined_pairs = loader.load_joined_pairs(structures)

        assert len(joined_pairs) == 3

        alignment = loader.load_alignment(structures)

        def get_stc_calls(index):
            return structures[index].converter.get_ind.call_args_list

        stc1_call_args_set = {call.args[0] for call in get_stc_calls(0)}
        for i in stc1_call_args_set:
            assert i.chain == "K"
        assert ("K", 31, "i") in stc1_call_args_set

        stc2_call_args_set = {call.args[0] for call in get_stc_calls(1)}
        for i in stc2_call_args_set:
            assert i.chain == "AB"

        stc3_call_args_set = {call.args[0] for call in get_stc_calls(2)}
        stc3_chains = {i.chain for i in stc3_call_args_set}
        assert stc3_chains == {"A", "B"}

    @pytest.mark.system
    def test_pair(self, alignments_dir):
        stc_path = alignments_dir / "structures" / "sars_pair.pdb"
        structures = get_structures_from_file(stc_path)
        path = alignments_dir / "pal" / "sars_pair.pal"

        with open(path) as fh:
            loader = PALLoader(fh)
        alignment = loader.load_alignment(structures)

        assert len(alignment.get_structures()) == 2
        assert alignment.get_inds_table().shape == (148, 2)

        expected_first = numpy.array([0, 2])
        numpy.testing.assert_equal(alignment.get_inds_table()[0], expected_first)

        expected_7_8 = numpy.array(
            [[8, 10], [9, 12]]
        )  # there is 1 residue gap in right structure
        numpy.testing.assert_equal(alignment.get_inds_table()[(8, 9), :], expected_7_8)

    def test_v2_v1(self, alignments_dir, artificial_multi_structures):
        structures = artificial_multi_structures
        path1 = alignments_dir / "pal" / "artificial_multi.pal"
        path2 = alignments_dir / "pal" / "artificial_multi_v2.pal"
        with open(path1) as fh1, open(path2) as fh2:
            loader1 = PALLoader(fh1)
            loader2 = PALLoader(fh2)

        al1 = loader1.load_alignment(structures)
        al2 = loader2.load_alignment(structures)

        numpy.testing.assert_array_equal(al1.get_inds_table(), al2.get_inds_table())
        assert al1.get_structures() == al2.get_structures()

    def test_load_chosen_structures(
        self, artificial_multi_path, artificial_multi_structures
    ):
        with open(artificial_multi_path) as fh:
            loader = PALLoader(fh)
        structures = artificial_multi_structures[:2]
        alignment = loader.load_alignment(structures)
        assert alignment.get_inds_table().shape == (5, 2)

    def test_load_lacking_chain(self, alignments_dir):
        path_wrong = alignments_dir / "pal" / "no_chain_when_needed.pal"
        with open(path_wrong) as fh:
            loader = PALLoader(fh)
        stc_path = alignments_dir / "structures" / "3g67.pdb"
        stc_3g67 = get_structures_from_file(stc_path)[0]
        mocked_structure = MagicMock()
        mocked_structure.chains = [MagicMock(chain_name="K")]
        mocked_structure.__getitem__.return_value = [MagicMock() for _ in range(5)]
        structures = [stc_3g67, mocked_structure]
        with pytest.raises(ValueError) as err_info:
            loader.load_alignment(structures)

        assert "matches none or more than one chain" in str(err_info.value)

        correct_path = alignments_dir / "pal" / "no_chain_when_needed_fixed.pal"
        with open(correct_path) as fh:
            loader = PALLoader(fh)
        alignment = loader.load_alignment(structures)
        assert len(alignment.get_structures()) == 2
        assert alignment.get_inds_table().shape == (5, 2)

    def test_self_aligned(self, alignments_dir, structure_1l0q):
        path_self_aligned = alignments_dir / "pal" / "1l0q_self.pal"
        with open(path_self_aligned) as fh:
            loader = PALLoader(fh)
        alignment = loader.load_alignment_mapping({"1l0q": structure_1l0q})
        assert alignment.get_inds_table().shape == (246, 2)
        alignment = Close(alignment).process()
        assert alignment.get_inds_table().shape == (41, 6)


class TestFASTALoader:
    @pytest.fixture(scope="session")
    def artificial_uneven_path(self, alignments_dir):
        path = alignments_dir / "fasta" / "uneven_multi.fasta"
        return path

    @pytest.fixture(scope="session")
    def artificial_multi_path(self, alignments_dir):
        path = alignments_dir / "fasta" / "artificial_multi.fasta"
        return path

    @pytest.fixture(scope="session")
    def artificial_multi_structures(self):
        structures = [MagicMock(name=f"mol{i}") for i in range(6)]
        for stc in structures[:4]:
            mer_list = [MagicMock(ind=i) for i in range(10, 15)]
            stc.__getitem__.return_value = mer_list
            stc.__iter__.return_value = iter(mer_list)
        structures[-2].__getitem__.return_value = [MagicMock(ind=1), MagicMock(ind=2)]
        structures[-1].__getitem__.return_value = [
            MagicMock(ind=i) for i in range(-8, -3)
        ]
        return structures

    def test_artificial_multi(self, artificial_multi_path, artificial_multi_structures):
        with open(artificial_multi_path) as fh:
            loader = FASTALoader(fh)
        alignment = loader.load_alignment(artificial_multi_structures)

        ranges = loader.read_metadata()["ranges"]
        assert ranges["mol1.A"] == "[A:2-4]"

        for label in ("mol1.A", "mol2", "mol3_chainAB", "mol4", "mol5", "mol6"):
            assert label in loader._structure_labels

        assert alignment.get_inds_table().shape == (4, 6)

        expected = numpy.array(
            [
                [DASH, DASH, DASH, 10, 1, -8],
                [10, 11, DASH, 11, 2, -7],
                [11, 12, 10, 12, 1, -6],
                [12, 13, 11, 13, DASH, DASH],
            ]
        )
        numpy.testing.assert_equal(alignment.get_inds_table(), expected)

    @pytest.mark.system
    def test_from_dama(self, alignments_dir):
        file_name = "mers_sars2.%s"
        fasta_path = alignments_dir / "fasta" / (file_name % "fasta")
        pal_path = alignments_dir / "pal" / (file_name % "pal")
        mers_path = alignments_dir / "structures" / "mers.pdb"
        (mers_stc,) = get_structures_from_file(mers_path)
        sars_path = alignments_dir / "structures" / "sars2.pdb"
        (sars_stc,) = get_structures_from_file(sars_path)
        structures_map = {
            "MOL1": mers_stc,
            "MOL2": sars_stc,
        }

        with open(fasta_path) as fh:
            fasta_loader = FASTALoader(fh)
        fasta = fasta_loader.load_alignment_mapping(structures_map)
        with open(pal_path) as fh:
            pal_loader = PALLoader(fh)
        pal = pal_loader.load_alignment_mapping(structures_map)

        numpy.testing.assert_equal(fasta.get_inds_table(), pal.get_inds_table())
        assert fasta.get_structures() == pal.get_structures()

    def test_load_chosen_structures(
        self, artificial_multi_path, artificial_multi_structures
    ):
        with open(artificial_multi_path) as fh:
            loader = FASTALoader(fh)
        alignment = loader.load_alignment(artificial_multi_structures[:2])
        assert alignment.get_inds_table().shape == (3, 2)

    def test_uneven(self, artificial_uneven_path, artificial_multi_structures):
        with open(artificial_uneven_path) as fh:
            loader = FASTALoader(fh)
        stc1, stc2, stc3, _, _, _ = artificial_multi_structures
        stc_map = {
            "mol1": stc1,
            "mol2": stc2,
            "mol3": stc3,
        }

        with pytest.raises(ValueError) as err_info:
            loader.load_alignment(artificial_multi_structures)

        assert "uneven" in str(err_info.value)

        stc_map.pop("mol2")
        with open(artificial_uneven_path) as fh:
            loader = FASTALoader(fh)
        alignment = loader.load_alignment_mapping(stc_map)
        assert alignment.get_inds_table().shape == (4, 2)


class TestXMLLoader:
    @pytest.fixture(scope="session")
    def sisy_sample(self, alignments_dir):
        path = alignments_dir / "xml" / "AL00051392.xml"
        return path

    @pytest.fixture(scope="session")
    def sisy_structures(self, alignments_dir):
        stcs = []
        for stc in ("1xi3A", "2tpsA", "1yadA"):
            pth = alignments_dir / "structures" / f"{stc}.pdb"
            (stc_obj,) = get_structures_from_file(pth)
            stcs.append(stc_obj)
        return stcs

    def test_sisy_sample(self, sisy_sample, sisy_structures):
        with open(sisy_sample) as fh:
            loader = XMLLoader(fh)

        alignment = loader.load_alignment(sisy_structures)

        assert len(alignment) == 168
        numpy.testing.assert_array_equal(alignment.get_inds_table()[0], [6, 17, 1])
