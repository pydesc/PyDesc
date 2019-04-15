import pytest
import os.path

from conftest import TEST_STRUCTURES_DIR
from conftest import PDB_FILES_WITH_TYPE

from pydesc.structure import StructureLoader
from pydesc.contacts import ContactMapCalculator
from pydesc.contacts.contacts import RcContact


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE)
def test_rc_contact_map(type_, struc_file):

    sl = StructureLoader()

    structures = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, type_, struc_file))

    for structure in structures:
        cm_calc = ContactMapCalculator(structure, contact_criterion_obj=RcContact())
        cm = cm_calc.calculate_contact_map()

        assert len(cm) > 0, 'No contacts in structure %s' % str(structure)
        for (k1, k2), v in cm:
            assert (structure[k1].rc - structure[k2].rc).calculate_length() < 10