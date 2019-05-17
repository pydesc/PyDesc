import pytest
import os.path
import pickle

from tests.conftest import TEST_STRUCTURES_DIR
from tests.conftest import TEST_CMAPS_DIR
from tests.conftest import PDB_FILES_WITH_TYPE
from tests.conftest import PDB_FILES_DICT

from pydesc.structure import StructureLoader
from pydesc.contacts import ContactMapCalculator
from pydesc.contacts.contacts import RcContact
from pydesc.config import ConfigManager

ConfigManager.warnings_and_exceptions.set("quiet", True)

PROTS_DIR = 'prots_only'


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


@pytest.mark.parametrize('structure_file', PDB_FILES_DICT[PROTS_DIR])
def test_golden_standard_pydesc_criterion_protein(structure_file):
    structure_name = os.path.splitext(structure_file)[0]
    with open(os.path.join(TEST_CMAPS_DIR, "%s_default.cmp" % structure_name), 'rb') as fh:
        golden_cmap_dict = pickle.load(fh)

    sl = StructureLoader()

    structure = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, PROTS_DIR, structure_file))[0]

    cm_calc = ContactMapCalculator(structure)
    cm = cm_calc.calculate_contact_map()
    res = {frozenset(k): v for k, v in cm._contacts.items()}

    assert golden_cmap_dict == res


@pytest.mark.parametrize('structure_file', PDB_FILES_DICT[PROTS_DIR])
def test_golden_standard_rc_protein(structure_file):
    structure_name = os.path.splitext(structure_file)[0]
    with open(os.path.join(TEST_CMAPS_DIR, "%s_rc.cmp" % structure_name), 'rb') as fh:
        golden_cmap_dict = pickle.load(fh)

    sl = StructureLoader()

    structure = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, PROTS_DIR, structure_file))[0]

    cm_calc = ContactMapCalculator(structure, contact_criterion_obj=RcContact())
    cm = cm_calc.calculate_contact_map()

    res = {frozenset(k): v for k, v in cm._contacts.items()}

    assert golden_cmap_dict == res
