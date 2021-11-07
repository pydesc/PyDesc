import pytest

from pydesc.api.structure import get_structures_from_file
from pydesc.cydesc.fitdesc import FitDesc


def test_fit_protein_gold_standard(structures_dir):
    (stc,) = get_structures_from_file(structures_dir / "prots_only" / "4NJ6.pdb")
    segment = stc[0:7]
    fitter = FitDesc(segment, stc)
    res = fitter.fitdesc(2.0, 5)

    # gold standard comes from last working version of PyDesc in Python 2.7
    expected_rmsds = [
        1.381876870709675e-07,
        0.1570712774991989,
        0.184758260846138,
        0.1918950080871582,
        0.2038000226020813,
    ]
    rmsds, alignments, _ = zip(*res)

    for rmsd, expected in zip(rmsds, expected_rmsds):
        assert rmsd == expected

    for al in alignments:
        assert len(al) == 8
