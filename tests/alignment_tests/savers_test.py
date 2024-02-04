from io import StringIO
from pathlib import Path

import numpy
import pytest

from pydesc.alignment.base import DASH
from pydesc.alignment.base import Alignment
from pydesc.alignment.loaders import PALLoader
from pydesc.alignment.savers import CSVSaver
from pydesc.alignment.savers import FASTASaver
from pydesc.alignment.savers import PALSaver
from pydesc.api.structure import get_structures_from_file


@pytest.fixture(scope="session")
def structure1(structures_dir):
    path = Path(structures_dir) / "prots_only" / "2BLL.pdb"
    (structure1,) = get_structures_from_file(str(path))
    return structure1


@pytest.fixture(scope="session")
def structure2(structures_dir):
    path = Path(structures_dir) / "prots_only" / "2BLL.pdb"
    (structure1,) = get_structures_from_file(str(path))
    return structure1


@pytest.fixture(scope="session")
def structure3(structures_dir):
    pth = Path(structures_dir) / "rna_only" / "1KIS.pdb"
    structure = get_structures_from_file(str(pth))
    return structure[0]


@pytest.fixture(scope="session")
def structure4(structures_dir):
    pth = Path(structures_dir) / "rna_only" / "1KIS.pdb"
    structure = get_structures_from_file(str(pth))
    return structure[0]


def make_trivial_pair_w_dashes(stc1, stc2):
    arr = numpy.array(
        [
            [1, 1],
            [DASH, 2],
            [3, DASH],
            [4, 4],
        ]
    )
    alignment = Alignment((stc1, stc2), arr)
    return alignment


def make_full_trivial_alignment(stc1, stc2):
    arr = numpy.array(
        [
            [i.ind, j.ind]
            for i, j in zip(stc1, stc2)
            if i.is_chainable() and j.is_chainable()
        ]
    )
    alignment = Alignment((stc1, stc2), arr)
    return alignment


class TestCSVSaver:
    @pytest.mark.parametrize(
        "names, expected_header", [(None, "2BLL\t2BLL"), (["A", "B"], "A\tB")]
    )
    def test_trivial(self, names, expected_header, structure1, structure2):
        saver = CSVSaver()
        alignment = make_trivial_pair_w_dashes(structure1, structure2)

        stream = StringIO()
        if names is not None:
            names = {structure1: names[0], structure2: names[1]}
        saver.save(stream, alignment, names=names)

        stream.seek(0)
        header = stream.readline()
        assert header.strip() == expected_header
        lines = [i.strip() for i in stream.readlines()]

        assert "-" in lines[1]
        assert "-" in lines[2]

        line0_mer = "A:317:R"
        assert lines[0].startswith(line0_mer)
        assert lines[0].endswith(line0_mer)

        line3_mer = "A:320:I"
        assert lines[3].startswith(line3_mer)
        assert lines[3].endswith(line3_mer)


class TestFASTASaver:
    def test__get_segments_simple(self):
        inds = [1, 2, 3, 10, 11, 12, 13]
        result = FASTASaver._get_segments(inds)
        assert result[1] == 3
        assert result[10] == 13

    def test_get_segments_empty(self):
        result = FASTASaver._get_segments([])
        assert not result

    def test_get_segments_trailing(self):
        inds = [1, 2, 3, 42]
        result = FASTASaver._get_segments(inds)
        assert result[1] == 3
        assert result[42] == 42
        inds = [1, 42, 43, 44]
        result = FASTASaver._get_segments(inds)
        assert result[1] == 1
        assert result[42] == 44

    def test_trivial(self, structure1, structure2):
        saver = FASTASaver()
        alignment = make_trivial_pair_w_dashes(structure1, structure2)

        stream = StringIO()
        saver.save(stream, alignment)

        expected_sequences = {
            structure1: "R-lI",
            structure2: "Rv-I",
        }
        expected_segments = {
            structure1: "[A:317-317, A:319-320]",
            structure2: "[A:317-318, A:320-320]",
        }

        stream.seek(0)
        for stc in alignment.get_structures():
            header = stream.readline()
            sequence = stream.readline()
            expected_sequence = expected_sequences[stc]
            assert sequence.strip() == expected_sequence
            assert header.startswith(">")
            assert stc.name in header
            expected_segment = expected_segments[stc]
            assert expected_segment in header

        assert not stream.readline()

    def test_full_trivial(self, structure3, structure4):
        saver = FASTASaver()
        alignment = make_full_trivial_alignment(structure3, structure4)

        stream = StringIO()
        saver.save(stream, alignment)

        stream.seek(0)
        header1 = stream.readline()
        seq1 = stream.readline()
        header2 = stream.readline()
        seq2 = stream.readline()

        expected_range = "[A:2-16, B:18-32]"
        assert expected_range in header1
        assert expected_range in header2

        assert len(seq1.strip()) == len(seq2.strip()) == 30

    def test_wrap_sequence(self, structure3, structure4):
        saver = FASTASaver(wrap_at=15)
        alignment = make_full_trivial_alignment(structure3, structure4)

        stream = StringIO()
        saver.save(stream, alignment)

        stream.seek(0)
        stream.readline()
        seq1_line1 = stream.readline()
        seq1_line2 = stream.readline()

        assert len(seq1_line1.strip()) == 15
        assert len(seq1_line2.strip()) == 15


class TestPALSaver:
    def test_trivial_pair(self, structure1, structure2):
        saver = PALSaver()
        alignment = make_full_trivial_alignment(structure1, structure2)

        tmp_file = StringIO()
        saver.save(tmp_file, alignment)

        tmp_file.seek(0)
        string = tmp_file.read()
        assert string != ""
        expected_header = "v:2.0\n2\n2BLL\n2BLL\n"
        assert string.startswith(expected_header)

        assert ">2BLL 2BLL" in string

        expected_line = "A:316 -- A:657 <--> A:316 -- A:657  (0,1)"
        assert expected_line in string

    def test_dashes_pair(self, structure1, structure2):
        saver = PALSaver()
        alignment = make_trivial_pair_w_dashes(structure1, structure2)

        stream = StringIO()
        saver.save(stream, alignment)

        stream.seek(0)
        string = stream.read()

        expected_line1 = "A:317 -- A:317 <--> A:317 -- A:317  (0,1)"
        expected_line2 = "A:320 -- A:320 <--> A:320 -- A:320  (0,1)"
        assert expected_line1 in string
        assert expected_line2 in string
        assert "@2" in string

    def test_names_mapping(self, structure1, structure2, structure3):
        names = {
            structure1: "AA",
            structure2: "BB",
        }
        saver = PALSaver()
        arr, _ = numpy.indices((10, 3))
        structures = (structure1, structure2, structure3)
        alignment = Alignment(structures, arr)
        stream = StringIO()

        saver.save(stream, alignment, names=names)

        stream.seek(0)
        string = stream.read()

        expected_header = "v:2.0\n3\nAA\nBB\n1KIS"
        assert string.startswith(expected_header)

        assert ">AA BB" in string
        assert ">AA 1KIS" in string
        assert ">BB 1KIS" in string

    def test_trivial_multiple(self, structure1, structure2, structure3):
        saver = PALSaver()
        arr, _ = numpy.indices((10, 3))
        structures = (structure1, structure2, structure3)
        alignment = Alignment(structures, arr)
        stream = StringIO()

        saver.save(stream, alignment)

        stream.seek(0)
        string = stream.read()

        expected_header = "v:2.0\n3\n2BLL\n2BLL\n1KIS"
        assert string.startswith(expected_header)

        expected_line1 = "A:316 -- A:325 <--> A:316 -- A:325  (0,1)"
        expected_line2 = "A:316 -- A:325 <--> A:1 -- A:10  (0,1)"

        assert expected_line1 in string
        assert expected_line2 in string

    def test_multiple_chains(self, structure1, structure2, structure3):
        saver = PALSaver()
        arr, _ = numpy.indices((30, 3))
        structures = (structure1, structure2, structure3)
        alignment = Alignment(structures, arr)
        stream = StringIO()

        saver.save(stream, alignment)

        stream.seek(0)
        string = stream.read()

        expected_header = "v:2.0\n3\n2BLL\n2BLL\n1KIS"
        assert string.startswith(expected_header)

        expected_stc1_stc2 = ">2BLL 2BLL\n@1\nA:316 -- A:345 <--> A:316 -- A:345  (0,1)"
        assert expected_stc1_stc2 in string

        expected_stc1_stc3 = (
            ">2BLL 1KIS\n@2\n"
            "A:316 -- A:331 <--> A:1 -- A:16  (0,1)\n"
            "A:332 -- A:345 <--> B:17 -- B:30  (0,1)\n"
        )
        assert expected_stc1_stc3 in string

    @pytest.mark.system
    def test_read_written(self, structure1, structure4, tmp_path):
        alignment = make_full_trivial_alignment(structure1, structure4)
        saver = PALSaver()

        save_stream = StringIO()
        saver.save(save_stream, alignment)

        save_stream.seek(0)
        tmp_file = tmp_path / "test.pal"
        tmp_file.write_text(save_stream.read())
        with open(tmp_file) as fh:
            loader = PALLoader(fh)

        new_alignment = loader.load_alignment((structure1, structure4))

        numpy.testing.assert_array_equal(
            alignment.get_inds_table(), new_alignment.get_inds_table()
        )
        assert alignment.get_structures() == new_alignment.get_structures()
