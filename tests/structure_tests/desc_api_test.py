import pytest

from pydesc.api.criteria import get_rc_distance_criterion
from pydesc.api.descriptor import create_descriptor
from pydesc.api.descriptor import create_descriptors
from pydesc.api.structure import get_structures_from_file
from pydesc.contacts.maps import ContactMapCalculator


@pytest.fixture(scope="function")
def protein_cmap(protein_file):
    (s,) = get_structures_from_file(protein_file)
    cmc = ContactMapCalculator(s, get_rc_distance_criterion())
    cm = cmc.calculate_contact_map()
    return s, cm


def test_build_descriptor(protein_cmap):
    s, cm = protein_cmap
    descriptor = create_descriptor(s, s[22], cm)
    assert s[22] in descriptor
    assert len(descriptor.segments) == 3
    assert len(descriptor.contacts) == 11


def test_create_descriptors(protein_cmap):
    s, cm = protein_cmap
    descs = create_descriptors(s, cm)
    for mer, desc in zip(s, descs):
        if desc is None:
            continue
        assert mer in desc
        assert mer == desc.central_element.central_mer
