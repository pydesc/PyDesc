import os.path

from pydesc.structure import StructureLoader
from tests.conftest import TEST_STRUCTURES_DIR
from pydesc.contacts.geometrical import PointsDistanceCriterion


def test_point_distance():
    sl = StructureLoader()
    path_str = os.path.join(TEST_STRUCTURES_DIR, "rna_only", "1KIS.pdb")
    stc, = sl.load_structures(path=path_str)

    mer1 = stc[18]  # B:19
    mer2 = stc[31]  # B:32
    mer3 = stc[4]   # A:5
    # B:19.rc is 10.5 or so far from B:32.rc
    # B:19.rc to A:5 is way more

    crt = PointsDistanceCriterion("rc", 12, 0)

    res=crt.calculate_contacts(stc)

    assert res[mer1.ind, mer2.ind] == 2
    assert res[mer1.ind, mer3.ind] == 0

