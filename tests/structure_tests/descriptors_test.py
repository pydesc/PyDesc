import pytest

from pydesc.api.criteria import get_rc_distance_criterion
from pydesc.api.structure import get_structures_from_file
from pydesc.contacts.maps import ContactMapCalculator
from pydesc.structure.descriptors import DescriptorBuilder
from pydesc.structure.descriptors import DescriptorBuilderDriver
from pydesc.structure.descriptors import ElementFactory
from pydesc.structure.topology import ElementChainable
from pydesc.warnexcept import ElementCreationFail


class TestElementBuilder:
    def test_build(self, any_structure_file):
        structure, *_ = get_structures_from_file(any_structure_file)
        for chain in structure.chains:
            chainable = [i for i in chain if i.is_chainable()]
            for mer in chainable[2:-2]:
                try:
                    elem = ElementFactory.build(structure, mer)
                except ElementCreationFail:
                    next_mer = mer.next_mer
                    prev_mer = mer.prev_mer
                    neighbours = [prev_mer, next_mer]
                    if None not in neighbours:
                        new = [prev_mer.prev_mer, next_mer.next_mer]
                        neighbours.extend(new)
                    assert None in neighbours
                else:
                    assert type(elem) is ElementChainable

            for mer in chainable[:2] + chainable[-2:]:
                with pytest.raises(ElementCreationFail):
                    ElementFactory.build(structure, mer)


class TestProteinDescriptor:
    def test_build_rc_criterion(self, protein_file):
        """Test building descriptors with default side chain geometrical
        center distance criterion."""
        (s,) = get_structures_from_file(protein_file)

        if max(list(map(len, s.chains))) < 10:
            pytest.skip("Structure %s has chains below 10 mers long, " "thus skipping.")

        cmc = ContactMapCalculator(s, get_rc_distance_criterion())
        cm = cmc.calculate_contact_map()

        desc_builder = DescriptorBuilder()
        dbd = DescriptorBuilderDriver(desc_builder)

        descs = []
        descriptable_mers = []
        failed = 0
        for mer in s:
            try:
                element = ElementFactory().build(s, mer)
            except ElementCreationFail:
                failed += 1
                continue
            desc = dbd.build(s, element, cm)
            descs.append(desc)
            descriptable_mers.append(mer)

        assert len(descs) > 0, "No descriptors for structure %s" % protein_file

        for mer, desc in zip(descriptable_mers, descs):
            if desc is None:
                continue
            central_element = desc.central_element
            central_mer = tuple(central_element)[2]
            assert len(desc) > 0, "failed for %s" % mer
            assert mer in desc
            assert len(desc.contacts) > 0
            # keep it independent from structures and mers comparison
            for contact in desc.contacts:
                contact_c_inds = [tuple(i)[2].ind for i in contact.elements]
                assert central_mer.ind in contact_c_inds
