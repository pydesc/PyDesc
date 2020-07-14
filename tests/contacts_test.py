import os.path

import pytest

from pydesc.contacts import NxContact
from pydesc.contacts import PrcContact
from pydesc.contacts import RcContact
from pydesc.contacts import RingCenterContact
from pydesc.structure import StructureLoader

CRITERIA = [RcContact, NxContact, PrcContact, RingCenterContact]


class TestNucleotides:
    @pytest.mark.parametrize("criterion", CRITERIA)
    def test_point_distance(self, structures_dir, criterion):
        sl = StructureLoader()
        path_str = os.path.join(structures_dir, "rna_only", "1KIS.pdb")
        (stc,) = sl.load_structures(path=path_str)

        mer1 = stc[18]  # B:19
        mer2 = stc[31]  # B:32
        mer3 = stc[4]  # A:5

        crt = criterion(distance_threshold=12.5)

        val_close = crt.is_in_contact(mer1, mer2)
        val_far = crt.is_in_contact(mer1, mer3)

        dist_close = crt.calculate_distance(mer1, mer2)
        dist_far = crt.calculate_distance(mer1, mer3)

        assert val_close == 2
        assert val_far == 0

        assert dist_far > 12.5
        assert dist_close <= 12.5
