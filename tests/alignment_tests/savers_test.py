from io import StringIO
from pathlib import Path

import numpy
import pytest

from pydesc.alignment.base import DASH
from pydesc.alignment.base import PairAlignment
from pydesc.alignment.savers import CSVSaver
from pydesc.api.structure import get_structures_from_file


class TestCSVSaver:
    @pytest.mark.parametrize(
        "names, expected_header", [(None, "2BLL\t2BLL"), (["A", "B"], "A\tB")]
    )
    def test_trivial(self, names, expected_header, structures_dir):
        path = Path(structures_dir) / "prots_only" / "2BLL.pdb"
        saver = CSVSaver()
        (structure1,) = get_structures_from_file(str(path))
        (structure2,) = get_structures_from_file(str(path))

        arr = numpy.array([[1, 1], [DASH, 2], [3, DASH], [4, 4],])
        trivial_alignment = PairAlignment((structure1, structure2), arr)

        stream = StringIO()
        if names is not None:
            names = {structure1: names[0], structure2: names[1]}
        saver.save(stream, trivial_alignment, names=names)

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
