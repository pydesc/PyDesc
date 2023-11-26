import pytest

from pydesc.api.criteria import get_rc_distance_criterion
from pydesc.api.descriptor import create_descriptor
from pydesc.api.descriptor import create_descriptors
from pydesc.api.structure import get_structures_from_file
from pydesc.chemistry.base import Mer
from pydesc.contacts.maps import ContactMapCalculator
from pydesc.selection import AtomSetSubclass


@pytest.fixture(scope="function")
def protein_cmap(protein_file):
    (s,) = get_structures_from_file(protein_file)
    cmc = ContactMapCalculator(s, get_rc_distance_criterion())
    cm = cmc.calculate_contact_map()
    return s, cm


def test_build_descriptor(protein_cmap):
    s, cm = protein_cmap
    most_contacted_mer, max_contacts = max(
        [(mer, len(cm[mer.ind])) for mer in s], key=lambda x: x[1]
    )
    descriptor = create_descriptor(s, most_contacted_mer, cm)
    assert most_contacted_mer in descriptor
    assert most_contacted_mer in descriptor.central_element
    assert len(descriptor.central_element) == 5
    assert len(descriptor.segments) > 1
    assert 1 < len(descriptor.contacts) < max_contacts


def test_create_descriptors(protein_cmap):
    s, cm = protein_cmap
    descs = create_descriptors(s, cm)
    for mer, desc in zip(s, descs):
        if desc is None:
            continue
        assert mer in desc
        assert mer == desc.central_element.central_mer
