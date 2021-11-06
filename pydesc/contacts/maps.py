# Copyright 2017 Tymoteusz Oleniecki
#
# This file is part of PyDesc.
#
# PyDesc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyDesc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.

"""Contact maps and auxiliary classes providing method to create maps."""

from functools import wraps

import numpy
from scipy.sparse import dok_matrix

from pydesc.contacts.base import create_empty_contact_map


def same_structure_only(method):
    @wraps(method)
    def method_wrapper(self, other, *args, **kwargs):
        if self.converter != other.converter:
            msg = (
                "Sum of maps of contacts calculated for different structures is not"
                " supported."
            )
            raise ValueError(msg)
        return method(self, other, *args, **kwargs)

    return method_wrapper


class AbstractContactMapCalculator:
    """Abstract superclass responsible for calculating contact (frequency) maps.

    Args:
        structure: structure instance.
        contact_criterion: instance of contact criterion.

    """

    def __init__(self, structure, contact_criterion, selections=(None, None)):
        self.contact_criterion = contact_criterion
        self.structure = structure
        self.selection1, self.selection2 = selections

    def _calculate_contact_matrix(self):
        if self.selection1 is not None:
            structure1 = self.selection1.create_structure(self.structure)
            structure2 = self.selection2.create_structure(self.structure)
            contacts_mtx = self.contact_criterion.calculate_inter_contacts(
                structure1, structure2
            )
        else:
            contacts_mtx = self.contact_criterion.calculate_contacts(self.structure)
        contacts_mtx = dok_matrix(contacts_mtx)
        contacts_mtx.setdiag(0)

        return contacts_mtx


class ContactMapCalculator(AbstractContactMapCalculator):
    """Class responsible for calculating contact maps."""

    def calculate_contact_map(self):
        """Perform calculation of contact map for structure passed to initialization.

        Returns:
            ContactMap: object storing values of all contacts defined by contact
            criterion.

        """
        contacts_mtx = self._calculate_contact_matrix()
        converter = self.structure.derived_from.converter
        return ContactMap(contacts_mtx, converter)


class FrequencyMapCalculator(AbstractContactMapCalculator):

    # TODO: add hook for vectorized calculations

    def calculate_frequency_map(self):
        original_frame = self.structure.derived_from.get_frame()
        self.structure.derived_from.set_frame(0)
        end_frame = self.structure.derived_from.get_n_frames()
        contact_mtx = self._calculate_contact_matrix().astype(numpy.float64)
        for frame_ind in range(1, end_frame):
            self.structure.derived_from.set_frame(frame_ind)
            contact_mtx += self._calculate_contact_matrix()
        contact_mtx /= 2
        frequency_map = FrequencyMap(contact_mtx, self.structure, end_frame)
        self.structure.derived_from.set_frame(original_frame)
        return frequency_map


class DescriptorMapFactory:
    def __init__(self, descriptor):
        self.descriptor = descriptor

    def create_contact_map(self):
        contact_map = create_empty_contact_map(self.descriptor)
        for structural_contact in self.descriptor.contacts:
            ind1, ind2 = [
                element.central_mer.ind for element in structural_contact.elements
            ]
            contact_map[ind1, ind2] = structural_contact.value
        cmap = ContactMap(contact_map, self.descriptor.derived_from.converter)
        return cmap


class AbstractContactMap:
    """Abstract superclass for maps of contacts."""

    def __init__(self, contacts_mtx, converter):
        self.converter = converter
        self._contacts = contacts_mtx

    def __iter__(self):
        """Returns iterator that runs over all non-zero contacts in contact map."""
        return iter(self._contacts.items())

    def __len__(self):
        """Return number of non-zero contacts in this contact map."""
        return len(self._contacts)

    def get_dok_matrix(self):
        """Return contacts as scipy dok matrix."""
        return self._contacts

    def _get_atom_set_contacts(self, ind):
        """Get values of non-zero contacts of single atom set."""
        contacts = [(ind, value) for (_, ind), value in self._contacts[ind].items()]
        return contacts

    def _get_contact_value(self, ind1, ind2):
        """Return value of contact between atom sets of two given inds."""
        return self._contacts[ind1, ind2]

    def to_string(self, stream_out):
        """Dumps pairs of atom sets in contacts to a given file-like object in CSV
        format.

        Args:
          stream_out: file-like object (opened file or StringIO etc.).

        """
        with stream_out:
            for (ind1, ind2), value in self:
                pdb_id1 = self.converter.get_pdb_id(ind1)
                pdb_id2 = self.converter.get_pdb_id(ind2)
                line = f"{pdb_id1}\t{pdb_id2}\t{value}\n"
                stream_out.write(line)


class ContactMap(AbstractContactMap):
    """Map of contacts present in a given (sub)structure.

    Args:
        contacts_mtx(scipy.sparse.dok_matrix): matrix storing contact values.
        converter: instance of structure or substructure.

    """

    def __init__(self, contacts_mtx, converter):
        mtx = dok_matrix(contacts_mtx, dtype=numpy.uint8)
        super().__init__(mtx, converter)

    def __getitem__(self, item):
        try:
            return self.get_contact_value(*item)
        except TypeError:
            return self.get_atom_set_contacts(item)

    @same_structure_only
    def combine(self, other):
        """Return new contact map storing higher contact values from two maps.

        Args:
            other(ContactMap): contact map to combine with.

        Returns:
            ContactMap: map of maximal values of two maps.

        """
        matrix = other.get_dok_matrix().maximum(self.get_dok_matrix())
        return ContactMap(matrix, self.converter)

    @same_structure_only
    def get_relative_compliment_map(self, other):
        """Return relative compliment of other in self (self / other; contacts from
        self not occurring in other).

        Args:
            other(ContactMap): contact map to combine with.

        Returns:
            ContactMap: map storing only relative compliment of other in self.

        """
        matrix = self.get_dok_matrix().copy()
        other_matrix = other.get_dok_matrix()
        matrix[other_matrix == 2] = 0
        matrix[(other_matrix == 1).toarray() & (matrix != 0).toarray()] = 1
        return ContactMap(matrix, self.converter)

    def get_atom_set_contacts(self, ind):
        """Get values of non-zero contacts of single atom set.

        Args:
             ind(int): atom set ind.

        Returns:
             : sequence of tuples storing inds and contact values. Only non-zero
             values are present (1 or 2).

        """
        return self._get_atom_set_contacts(ind)

    def get_contact_value(self, ind1, ind2):
        """Return value of contact between atom sets of two given inds.

        Args:
            ind1(int): ind of first atom set.
            ind2(int): ind of second atom set.

        Returns:
            int: contact values (0, 1 or 2).

        """
        return self._get_contact_value(ind1, ind2)


class FrequencyMap(AbstractContactMap):
    """Map storing frequency of contacts in given (sub)structure or trajectory.

    Args:
        contacts_mtx(scipy.sparse.dok_matrix): matrix storing contact occurs.
        converter: (sub)structure or (sub)trajectory.
        n_frames(int): number of frames. Only necessary if second argument is
        structure, not trajectory. In such case structure is actually only needed for
        its converter (and that's typeist).

    """

    def __init__(self, contacts_mtx, converter, n_frames=None):
        mtx = dok_matrix(contacts_mtx, dtype=numpy.float64)
        super().__init__(mtx, converter)
        if n_frames is None:
            n_frames = converter.get_n_frames()
        self.n_frames = n_frames

    def __iter__(self):
        for pair, value in super().__iter__():
            yield pair, value / self.n_frames

    def get_contacts_occurs(self, ind):
        """Get non-zero occurs of contacts of single atom set.

        Args:
             ind(int): atom set ind.

        Returns:
             : sequence of tuples storing inds and contact occurs.

        """
        return self._get_atom_set_contacts(ind)

    def get_contacts_frequencies(self, ind):
        """Get non-zero frequencies of contacts of single atom set.

        Args:
             ind(int): atom set ind.

        Returns:
             : sequence of tuples storing inds and contact frequencies.

        """
        items = self._contacts[ind].items()
        contacts = [(ind, value / self.n_frames) for (_, ind), value in items]
        return contacts

    def get_contact_occurs(self, ind1, ind2):
        """Return number of occurs of contact between atom sets of two given inds.

        Args:
            ind1(int): ind of first atom set.
            ind2(int): ind of second atom set.

        Returns:
            float: number of occurs. Uncertain contacts are counted as 0.5.

        """
        return self._get_contact_value(ind1, ind2)

    def get_contact_frequency(self, ind1, ind2):
        """Return frequency of contact between atom sets of two given inds.

        Args:
            ind1(int): ind of first atom set.
            ind2(int): ind of second atom set.

        Returns:
            float: contact frequency (0.0-1.0).

        """
        value = self._get_contact_value(ind1, ind2)
        return value / self.n_frames
