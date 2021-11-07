import pytest

from pydesc.api.criteria import get_rc_distance_criterion
from pydesc.api.descriptor import create_descriptor
from pydesc.api.descriptor import create_descriptors
from pydesc.api.structure import get_structures_from_file
from pydesc.contacts.maps import ContactMapCalculator
from pydesc.selection import AtomSetSubclass
from pydesc.chemistry.base import Mer


@pytest.fixture(scope="function")
def protein_cmap(protein_file):
    (s,) = get_structures_from_file(protein_file)
    cmc = ContactMapCalculator(s, get_rc_distance_criterion())
    cm = cmc.calculate_contact_map()
    return s, cm


def test_build_descriptor(protein_cmap):
    s, cm = protein_cmap
    if len(s) < 15:
        pytest.xfail("Not enough mers to test descriptors")
    mer = tuple(AtomSetSubclass(Mer).create_structure(s))[4]
    descriptor = create_descriptor(s, mer, cm)
    assert mer in descriptor
    assert mer in descriptor.central_element
    assert len(descriptor.central_element) == 5
    assert len(descriptor.segments) >= 1
    assert len(descriptor.contacts) >= 1


def test_create_descriptors(protein_cmap):
    s, cm = protein_cmap
    descs = create_descriptors(s, cm)
    for mer, desc in zip(s, descs):
        if desc is None:
            continue
        assert mer in desc
        assert mer == desc.central_element.central_mer
