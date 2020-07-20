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

from scipy.sparse import dok_matrix


class ContactMapCalculator:
    """Class responsible for calculating contact maps.

    Args:
        structure: structure instance.
        contact_criterion: instance of contact criterion.

    """

    def __init__(self, structure, contact_criterion, selections=(None, None)):
        self.contact_criterion = contact_criterion
        self.structure = structure
        self.selection1, self.selection2 = selections

    def calculate_contact_map(self):
        """Perform calculation of contact map for structure passed to initialization.

        Returns:
            ContactMap: object storing values of all contacts defined by contact
            criterion.

        """
        if self.selection1 is not None:
            structure1 = self.selection1.create_structure(self.structure)
            structure2 = self.selection2.create_structure(self.structure)
            contacts_mtx = self.contact_criterion.calculate_inter_contacts(
                structure1, structure2)
        else:
            contacts_mtx = self.contact_criterion.calculate_contacts(self.structure)
        contacts_mtx = dok_matrix(contacts_mtx)
        contacts_mtx.setdiag(0)

        return ContactMap(contacts_mtx, self.structure)


class ContactMap:
    """Map of contacts present in a given (sub)structure.

    Args:
        contacts_mtx(scipy.sparse.dok_matrix): matrix storing contact values.
        structure: instance of structure or substructure.

    """

    def __init__(self, contacts_mtx, structure):
        self.structure = structure
        self._contacts = contacts_mtx

    def __iter__(self):
        """Returns iterator that runs over all non-zero contacts in contact map."""
        return iter(list(self._contacts.items()))

    def __len__(self):
        """Return number of non-zero contacts in this contact map."""
        return len(self._contacts)

    def __getitem__(self, item):
        try:
            return self.get_contact_value(*item)
        except TypeError:
            return self.get_mer_contacts(item)

    def get_mer_contacts(self, mer_id):
        """Get values of non-zero contacts of mer with given ind.

        Args:
             mer_id(int): mer ind.

        Returns:
             : sequence of tuples storing mer inds and contact values. Only non-zero
             values are present (1 or 2).

        """
        contacts = [(ind, value) for (_, ind), value in self._contacts[mer_id].items()]
        return contacts

    def get_contact_value(self, mer_id1, mer_id2):
        """Return value of contact between mers of two given inds.

        Args:
            mer_id1(int): ind of first mer.
            mer_id2(int): ind of second mer.

        Returns:
            int: contact values (0, 1 or 2).

        """
        return self._contacts[mer_id1, mer_id2]

    def to_string(self, stream_out):
        """Dumps pairs of mers in contacts to a given file-like object in CSV format.

        Args:
          stream_out: file-like object (opened file or StringIO etc.).

        """
        converter = self.structure.converter
        with stream_out:
            for (ind1, ind2), value in self:
                pdb_id1 = converter.get_pdb_id(ind1)
                pdb_id2 = converter.get_pdb_id(ind2)
                line = f"{pdb_id1}\t{pdb_id2}\t{value}\n"
                stream_out.write(line)
