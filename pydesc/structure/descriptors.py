import operator
from abc import ABCMeta
from abc import abstractmethod
from functools import reduce

from pydesc.config import ConfigManager
from pydesc.mers import Ion
from pydesc.mers import MerChainable
from pydesc.mers import Nucleotide
from pydesc.mers import Residue
from pydesc.selection import SelectionsUnion
from pydesc.selection import Set
from . import AbstractStructure
from . import Contact
from . import ElementChainable
from . import ElementOther
from . import Segment


class ElementFactory:

    @staticmethod
    def build(mer):
        """Static builder.

        Returns appropriate Element subclass instance.

        Argument:
        mer -- instance of pydesc.monomer.Monomer.
        """
        if isinstance(mer, MerChainable):
            return ElementChainable(mer)
        return ElementOther(mer)


class DescriptorBuilderDriver:

    def __init__(self, class_builder_map=None):
        if class_builder_map is None:
            class_builder_map = {
                Nucleotide: NucleotideDescriptorBuilder(),
                Residue: ProteinDescriptorBuilder(),
            }
        self.class_map = class_builder_map

    def get_builder(self, mer):
        try:
            return self.class_map[type(mer)]
        except KeyError:
            for klass, builder in self.class_map.items():
                if isinstance(mer, klass):
                    return self.class_map[klass]
            raise ValueError('No descriptor builder for given mer: %s' % mer)

    def build(self, central_element, contact_map):
        """Static builder.

        Arguments:
        central_element -- ElementChainable instance.
        contact_map -- instance of contact map. Set on None by default. If so,
         contact map of structure from which given element was derived from
         is taken.

        Returns appropriate descriptor.
        """
        builder = self.get_builder(central_element.central_monomer)
        builder.set_central_element(central_element)
        builder.create_contacts(central_element, contact_map)
        builder.set_mers()
        builder.set_segments()
        builder.set_elements()
        return builder.build()

    def create_descriptors(self, structure_obj, contact_map):
        """Static method that creates all possible Descriptor instances for
        a given (sub)structure.

        Arguments:
        structure_obj -- any AbstractStructure subclass instance
        contact_map -- ContactMap instance that will be submitted to
        AbstractDescriptor.__init__ method, initially set to None.
        """

        def mk_desc(mer):
            try:
                return self.build(ElementFactory.build(mer), contact_map)
            except (TypeError, ValueError, AttributeError):
                return

        return [mk_desc(mer) for mer in structure_obj]


class DescriptorBuilder:
    __metaclass__ = ABCMeta

    """Builder-like class preparing data to create descriptor."""

    def __init__(self):
        """Set initial attributes."""
        self.central_element = None
        self.mers = []
        self.contacts = []
        self.segments = []
        self.elements = []

    def build(self):
        """Build protein descriptor using set attributes."""
        return Descriptor(
            self.central_element,
            self.mers,
            self.elements,
            self.segments,
            self.contacts
        )
        """Build appropriate descriptor."""
        pass

    def set_mers(self):
        """Set mers on basis of calculated contacts."""
        mers = set([])
        for cnt in self.contacts:
            mers = mers.union(cnt)
        self.mers = sorted(mers, key=lambda x: x.ind)

    def set_central_element(self, element):
        """Assign given element as central element of descriptor to be
        built."""
        self.central_element = element

    @abstractmethod
    def create_contacts(self, central_element, contact_map):
        """Create contacts for given central element using given contact
        map. Abstract method, should be overwritten in subclasses."""
        pass

    def set_segments(self):  # TODO: nightmare code -- fix it please
        """Fills initial attr segments."""

        def add_segment(start, end):
            """Adds ElementChainable to segments list, or Segment if start and
            end are distant more then proper segment length.
            """
            element_length = ConfigManager.element.element_chainable_length
            segment = Segment(start, end)
            if len(segment) == element_length:
                central_mer = tuple(segment)[element_length // 2]
                segment = ElementChainable(central_mer)
            segments.append(segment)

        def reduce_tuple(presegment_1, presegment_2):
            """Reduces two tuple to one, if they contain subsequent mers."""
            if presegment_1[1] == presegment_2[0]:
                return presegment_1[0], presegment_2[1]
            add_segment(*presegment_1)
            return presegment_2

        def switch(mer):
            """Return False if mer is last in current structure, else return
            tuple (this mer, next mer)."""
            return (mer, mer.next_mer) if mer.next_mer in self.mers else False

        segments = []
        switched_mers = [switch(mer) for mer in self.mers]
        last_segment = reduce(
            reduce_tuple,
            [i for i in switched_mers if bool(i)])
        add_segment(*last_segment)
        self.segments = segments

    def set_elements(self):
        """Fills initial attrs elements and elements_values."""
        neighbours = {self.central_element.central_monomer: []}
        elements = {}
        for con in self.contacts:
            element1, element2 = con.elements
            for el1, el2 in zip((element1, element2), (element2, element1)):
                elements[el1.central_monomer] = el1
                neighbours.setdefault(el1.central_monomer, []).append(
                    el2.central_monomer)
        self.elements = sorted(elements.values(),
                               key=lambda elm: elm.central_monomer.ind)


class ProteinDescriptorBuilder(DescriptorBuilder):
    """Protein descriptor builder-like object. Implements create_contacts
    method creating contacts for central mer alone."""

    def create_contacts(self, element_obj, contact_map):
        """Create star topology contacts for central mer of given element.
        
        Arguments:
        element_obj -- ElementChainable instance.
        contact_map -- instance of contact map. Set on None by default. If 
        so, contact map of structure from which given element was derived 
        from is taken.
        """  # TODO: fix
        stc = element_obj.derived_from
        central_mer = element_obj.central_monomer
        contacts_ = sorted(contact_map.get_mer_contacts(central_mer.ind))

        def create_contact(input_):
            """Returns Contact or None if failed to create one.

            Argument:
            input_ -- tuple containing: ind_1, ind_2 -- mers inds, value --
            contact value.
            """
            ((ind_1, ind_2), value) = input_
            try:
                return Contact(
                    ElementFactory.build(central_mer),
                    ElementFactory.build(stc[ind_2]))
            except (ValueError,):
                # TypeError is raised during creation of Contacts based on
                # pairs of mers that cannot be used to create Elements
                # ValueError -- non chainable or too close to terminus
                return

        self.contacts = [create_contact(cnt) for cnt in contacts_]


class NucleotideDescriptorBuilder(DescriptorBuilder):

    # TODO: does that make any sense since we have proof that ion contacts
    #  are not influential?
    def create_contacts(self, element_obj, contact_map):
        """Create nucleotide descriptor contacts.

        That means star topology contacts between nucleotides and
        additionally all nucleotides that are in contact with any ion
        central nucleotide is in contact with (so contacts via ion).
        """

        def bld_con(der, *args):
            """Builds a contact."""
            return Contact(
                ElementFactory.build(der[args[0]]),
                ElementFactory.build(der[args[1]]))

        derf = element_obj.derived_from
        cmer = element_obj.central_monomer

        if isinstance(element_obj.central_monomer, Ion):
            raise TypeError

        if not contact_map:
            contact_map = element_obj.derived_from.contact_map

        contacts_objs = []
        contacts_inds = [cmer.ind]
        ions = []
        for i in contact_map.get_monomer_contacts(cmer.ind):
            if type(derf[i[1]]) != Ion:
                try:
                    contacts_objs.append(bld_con(derf, *i))
                    contacts_inds.append(i[1])
                except ValueError:
                    continue
            else:
                ions.append(derf[i[1]])
        ion_contacts = {}
        for ion in ions:
            cdist = (cmer.ring_center - ion.rc).calculate_length()
            for c_ind in contact_map.contacts[ion.ind]:
                if c_ind in contacts_inds:
                    continue
                diff_vec = derf[c_ind].ring_center - ion.rc
                dist = cdist + diff_vec.calculate_length()
                try:
                    cnto = bld_con(derf, ion.ind, c_ind)
                except ValueError:
                    continue
                try:
                    shortest_cnt = min(ion_contacts[c_ind],
                                       (dist, cnto, ion.ind),
                                       key=lambda x: x[0])
                except KeyError:
                    shortest_cnt = (dist, cnto, ion.ind)
                ion_contacts[c_ind] = shortest_cnt

        qualif_ion_cons = [bld_con(derf, i, cmer.ind) for i in
                           set([i[2] for i in ion_contacts.values()])]
        contacts_objs.extend(
            [i[1] for i in ion_contacts.values()] + qualif_ion_cons)

        return filter(bool, contacts_objs)


class Descriptor(AbstractStructure):
    """Representation of PyDesc spatial unit.

    It consists of the central Element and all Elements it is in contact with,
    according to a given ContactMap ContactCriterion.
    """

    __metaclass__ = ABCMeta

    def __init__(self, central_element, mers, elements, segments, contacts):
        """Descriptor constructor.

        Arguments:
        element_obj -- Element instance that becomes central_element.
        list_of_contact_obj -- list of Contact instances that are present in
        descriptor to be built.
        """
        AbstractStructure.__init__(self, central_element.derived_from)
        self.central_element = central_element
        self.contacts = contacts
        self._mers = mers
        self.segments = segments
        self.elements = elements

    @property
    def cm_pid(self):
        """Returns PDB id of central element central monomer as string."""
        return str(self.central_element.central_monomer.pid)

    def __repr__(self):
        return '<%s of %s:%s>' % (
            self.__class__.__name__, str(self.derived_from), self.cm_pid)

    def select(self):
        """Overridden select method.

        Returns union of range selections (for self segments) and set of all 
         remaining mers.
        """  # TODO: should be easier
        components = list(map(Segment.select, self.segments))
        other_mers = [
            mer for mer in self if
            mer not in reduce(operator.add, self.segments)]
        components.append(
            Set(map(operator.methodcaller('get_pdb_id'), other_mers)))
        return SelectionsUnion(components)
