import os.path
from unittest.mock import MagicMock

import numpy
import pytest

from pydesc.alignment.base import MultipleAlignment
from pydesc.alignment.base import PairAlignment
from pydesc.alignment.loaders import CSVLoader
from pydesc.alignment.loaders import DASH
from pydesc.alignment.loaders import FASTALoader
from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file


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

        assert alignment.inds.shape == (6, 4)
        expected_not_nans = [5, 4, 3, 5]
        real_not_nans = numpy.count_nonzero(alignment.inds != DASH, axis=0)
        numpy.testing.assert_equal(real_not_nans, expected_not_nans)

        assert len(alignment.structures) == 4

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

        assert isinstance(alignment, MultipleAlignment)

    def test_pair(self, csv_path):
        path = csv_path / "artificial_pair.csv"
        structures = [MagicMock() for _ in range(2)]
        with open(path) as file_:
            loader = CSVLoader(file_)
            alignment = loader.load_alignment(structures)
        assert isinstance(alignment, PairAlignment)

    def test_no_content(self, csv_path):
        path = csv_path / "empty_pair.csv"
        with open(path) as fh:
            loader = CSVLoader(fh)
        structures = [MagicMock(name=i) for i in range(2)]
        alignment = loader.load_alignment(structures)
        assert len(alignment) == 0
        assert len(alignment.structures) == 2

    def test_single_col(self, csv_path):
        path = csv_path / "single_col.csv"
        with open(path) as fh:
            loader = CSVLoader(fh)
        with pytest.raises(ValueError):
            loader.load_alignment([MagicMock()])

    def test_load_some(self, csv_path):
        path = csv_path / "artificial_multi.csv"
        with open(path) as fh:
            loader = CSVLoader(fh)
        stc_map = {k: MagicMock(name=k) for k in ("stc1", "stc3A")}
        alignment = loader.load_alignment_mapping(stc_map)
        assert isinstance(alignment, PairAlignment)


class TestPALLoader:
    @pytest.fixture(scope="session")
    def artificial_multi_structures(self):
        structures = [MagicMock() for _ in range(3)]
        mocked_chain = MagicMock(chain_name="K")
        structures[0].chains = [mocked_chain]
        mocked_mers = [MagicMock(ind=i) for i in range(5)]
        for structure in structures:
            structure.converter.get_ind.return_value = 0
            structure.__getitem__.return_value = mocked_mers
        return structures

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

        assert len(joined_pairs.pair_alignments) == 3

        alignment = loader.load_alignment(structures)
        alignment_inds = alignment.inds
        joined_pair_as_columns = joined_pairs.to_columns()
        joined_pairs_inds = joined_pair_as_columns.limit_to_structures(structures).inds
        numpy.testing.assert_array_equal(alignment_inds, joined_pairs_inds)

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

        assert isinstance(alignment, PairAlignment)
        assert alignment.inds.shape == (148, 2)

        expected_first = numpy.array([0, 2])
        numpy.testing.assert_equal(alignment.inds[0], expected_first)

        expected_7_8 = numpy.array(
            [[8, 10], [9, 12]]
        )  # there is 1 residue gap in right structure
        numpy.testing.assert_equal(alignment.inds[(8, 9), :], expected_7_8)

    def test_v2_v1(self, alignments_dir, artificial_multi_structures):
        structures = artificial_multi_structures
        path1 = alignments_dir / "pal" / "artificial_multi.pal"
        path2 = alignments_dir / "pal" / "artificial_multi_v2.pal"
        with open(path1) as fh1, open(path2) as fh2:
            loader1 = PALLoader(fh1)
            loader2 = PALLoader(fh2)

        al1 = loader1.load_alignment(structures)
        al2 = loader2.load_alignment(structures)

        numpy.testing.assert_array_equal(al1.inds, al2.inds)
        assert al1.structures == al2.structures

    def test_load_some(self, artificial_multi_path, artificial_multi_structures):
        with open(artificial_multi_path) as fh:
            loader = PALLoader(fh)
        structures = artificial_multi_structures[:2]
        alignment = loader.load_alignment(structures)
        assert isinstance(alignment, PairAlignment)

    def test_load_lacking_chain(self, alignments_dir):
        path_wrong = alignments_dir / "pal" / "no_chain_when_needed.pal"
        with open(path_wrong) as fh:
            loader = PALLoader(fh)
        stc_path = alignments_dir / "structures" / "3g67.pdb"
        stc_3g67 = get_structures_from_file(stc_path)[0]
        mocked_structure = MagicMock()
        mocked_structure.chains = [MagicMock(chain_name="K")]
        structures = [stc_3g67, mocked_structure]
        with pytest.raises(ValueError) as err_info:
            loader.load_alignment(structures)

        assert "matches none or more than one chain" in str(err_info.value)

        correct_path = alignments_dir / "pal" / "no_chain_when_needed_fixed.pal"
        with open(correct_path) as fh:
            loader = PALLoader(fh)
        alignment = loader.load_alignment(structures)
        assert isinstance(alignment, PairAlignment)


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
        structures = [MagicMock(name=f"mol{i}") for i in range(5)]
        for stc in structures[:4]:
            stc.__getitem__.return_value = [MagicMock(ind=i) for i in range(10, 15)]
        structures[-1].__getitem__.return_value = [MagicMock(ind=1), MagicMock(ind=2)]
        return structures

    def test_artificial_multi(self, artificial_multi_path, artificial_multi_structures):
        with open(artificial_multi_path) as fh:
            loader = FASTALoader(fh)
        alignment = loader.load_alignment(artificial_multi_structures)

        ranges = loader.read_metadata()["ranges"]
        assert ranges["mol1.A"] == "[A:2-4]"

        for label in ("mol1.A", "mol2", "mol3_chainAB", "mol4", "mol5"):
            assert label in loader._structure_labels

        assert alignment.inds.shape == (4, 5)

        expected = numpy.array(
            [
                [DASH, DASH, DASH, 10.0, 1.0],
                [10.0, 11.0, DASH, 11.0, 2.0],
                [11.0, 12.0, 10.0, 12.0, 1.0],
                [12.0, 13.0, 11.0, 13.0, DASH],
            ]
        )
        numpy.testing.assert_equal(alignment.inds, expected)

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

        numpy.testing.assert_equal(fasta.inds, pal.inds)
        assert fasta.structures == pal.structures

    def test_load_some(self, artificial_multi_path, artificial_multi_structures):
        with open(artificial_multi_path) as fh:
            loader = FASTALoader(fh)
        alignment = loader.load_alignment(artificial_multi_structures[:2])
        assert isinstance(alignment, PairAlignment)

    def test_uneven(self, artificial_uneven_path, artificial_multi_structures):
        with open(artificial_uneven_path) as fh:
            loader = FASTALoader(fh)
        stc1, stc2, stc3, _, _ = artificial_multi_structures
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
        assert isinstance(alignment, PairAlignment)
