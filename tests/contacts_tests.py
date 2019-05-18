import os.path
import pytest
from tests.conftest import TEST_STRUCTURES_DIR

from pydesc.contacts import (
    RcContact,
    NxContact,
    RingCenterContact,
    PrcContact,
)
from pydesc.structure import StructureLoader


class TestNucleotides:

    @pytest.mark.parametrize('crit', [RcContact, NxContact, PrcContact, RingCenterContact])
    def test_point_distance(self, crit):
        sl = StructureLoader()
        stc, = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, 'rna_only', '1KIS.pdb'))

        mer1 = stc[18]
        mer2 = stc[31]
        mer3 = stc[4]

        crt = crit(distance_threshold=10.)

        val_close = crt.is_in_contact(mer1, mer2)
        val_far = crt.is_in_contact(mer1, mer3)

        dist_close = crt.calculate_distance(mer1, mer2)
        dist_far = crt.calculate_distance(mer1, mer3)

        assert val_close == 2
        assert val_far == 0

        assert dist_far > 10
        assert dist_close <= 10
