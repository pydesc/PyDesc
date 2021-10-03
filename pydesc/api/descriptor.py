"""Convenience functions for creating descriptors."""
from pydesc.structure.descriptors import DescriptorBuilder
from pydesc.structure.descriptors import DescriptorBuilderDriver
from pydesc.structure.descriptors import ElementFactory
from pydesc.structure.topology import ElementOther
from pydesc.warnexcept import ElementCreationFail


def create_descriptor(
    structure, mer, contact_map, element_factory=None, builder_driver=None,
):
    element_factory, builder_driver = _create_factory_and_driver(
        element_factory, builder_driver
    )

    element = element_factory.build(structure, mer)
    descriptor = builder_driver.build(structure, element, contact_map)
    return descriptor


def create_descriptors(
    structure, contact_map, element_factory=None, builder_driver=None,
):
    def mk_desc(mer):
        try:
            central_element = element_factory.build(structure, mer)
        except (ElementCreationFail):
            return
        if isinstance(central_element, ElementOther):
            return
        descriptor = builder_driver.build(structure, central_element, contact_map)
        return descriptor

    element_factory, builder_driver = _create_factory_and_driver(
        element_factory, builder_driver
    )
    return [mk_desc(mer) for mer in structure]


def _create_factory_and_driver(
    element_factory=None, builder_driver=None,
):
    if element_factory is None:
        element_factory = ElementFactory()
    if builder_driver is None:
        descriptor_builder = DescriptorBuilder(element_factory)
        builder_driver = DescriptorBuilderDriver(descriptor_builder)
    return element_factory, builder_driver
