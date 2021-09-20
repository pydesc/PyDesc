"""Convenience functions for creating descriptors."""

from pydesc.structure.descriptors import DescriptorBuilder
from pydesc.structure.descriptors import DescriptorBuilderDriver
from pydesc.structure.descriptors import ElementFactory
from pydesc.warnexcept import ElementCreationFail


def create_descriptor(structure, mer, contact_map):
    element_factory = ElementFactory()
    descriptor_builder = DescriptorBuilder(element_factory)
    builder_driver = DescriptorBuilderDriver(descriptor_builder)

    element = element_factory.build(structure, mer)
    descriptor = builder_driver.build(structure, element, contact_map)
    return descriptor


def create_descriptors(structure, contact_map):
    def mk_desc(mer):
        try:
            central_element = element_factory.build(structure, mer)
        except (ElementCreationFail):
            return
        descriptor = builder_driver.build(structure, central_element, contact_map)
        return descriptor

    element_factory = ElementFactory()
    descriptor_builder = DescriptorBuilder(element_factory)
    builder_driver = DescriptorBuilderDriver(descriptor_builder)

    return [mk_desc(mer) for mer in structure]
