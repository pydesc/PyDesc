import pytest

from pydesc.api.criteria import get_rc_distance_criterion
from pydesc.api.structure import get_structures_from_file
from pydesc.contacts.maps import ContactMapCalculator
from pydesc.structure.descriptors import DescriptorBuilderDriver
from pydesc.structure.descriptors import ElementFactory
from pydesc.structure.topology import ElementChainable


class TestElementBuilder:
    def test_build(self, any_structure_file):
        structure = get_structures_from_file(any_structure_file)[0]
        for chain in structure.chains:
            chainable = [i for i in chain if i.is_chainable()]
            failed = 0
            for mer in chainable[2:-2]:
                try:
                    elem = ElementFactory.build(mer, structure)
                except ValueError:
                    failed += 1
                else:
                    assert type(elem) is ElementChainable
            assert failed <= 0.25 * len(chainable)

            for mer in chainable[:2] + chainable[-2:]:
                with pytest.raises(ValueError):
                    ElementFactory.build(mer, structure)


class TestProteinDescriptor:
    def test_build_rc_criterion(self, protein_file):
        """Test building descriptors with default side chain geometrical
        center distance criterion."""
        (s,) = get_structures_from_file(protein_file)

        if max(list(map(len, s.chains))) < 10:
            pytest.skip("Structure %s has chains below 10 mers long, " "thus skipping.")

        cmc = ContactMapCalculator(s, get_rc_distance_criterion())
        cm = cmc.calculate_contact_map()

        dbd = DescriptorBuilderDriver()

        descs = dbd.create_descriptors(s, cm)

        assert len([i for i in descs if i is not None]) > 0, (
            "No descriptors for structure %s" % protein_file
        )

        for mer, desc in zip(s, descs):
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
