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

"""
Classes that deal with contacts among mers present in PyDesc (sub)structures.

created: 28.04.2014 - , Tymoteusz 'hert' Oleniecki
"""

from scipy.sparse import dok_matrix

from pydesc.selection import Everything


class ContactMapCalculator:
    """Class responsible for calculating contact maps."""

    def __init__(self, structure_obj, contact_criterion_obj):
        self.contact_criterion = contact_criterion_obj
        self.structure = structure_obj

    def calculate_contact_map(self):
        """Return ContactMap for structure set during initialization.

        Contact is a tuple containing first and second Monomer, distance(s)
        between them and a contact value under the ContactMap criterion.
        """
        contacts_mtx = self.contact_criterion.calculate_contacts(self.structure)
        contacts_mtx = dok_matrix(contacts_mtx)
        contacts_mtx.setdiag(0)

        return ContactMap(contacts_mtx, self.structure)


class ContactMap:
    """Map of contacts present in a given (sub)structure."""

    def __init__(self, contacts_mtx, structure):
        """ContactMap constructor.

        Arguments:
        contacts_mtx -- sparse dict-like contact matrix representing
        contacts in biopolymer. 2 indicates contact,
        1 -- plausible contact.
        structure -- instance of any pydesc.structure.AbstractStructure
        subclass for which contact map is to be created.
        """
        self.structure = structure
        self.frame = None
        self._contacts = contacts_mtx

    def __iter__(self):
        """Returns iterator that runs over all contacts in contact map."""
        return iter(list(self._contacts.items()))

    def __len__(self):
        """Return length of self._contacts dok_matrix keys list."""
        return len(self._contacts)

    def __getitem__(self, item):
        try:
            return self.get_contact_value(*item)
        except TypeError:
            return self.get_mer_contacts(item)

    def get_mer_contacts(self, mer_id):
        """Returns list of given monomer contacts.

        Contact is a tuple containing given monomer ind, ind of monomer that
        stays in contacts with it and contact value in three-valued logic.
        Distance evaluation is based on settings in configuration manager.

        Arguments:
        mer_id -- monomer instance, PDB_id instance (or tuple containing
        proper values; see number converter docstring for more information),
        PyDesc ind or monomer index on a list of structure mers.
        """
        return list(self._contacts[mer_id].items())

    def get_contact_value(self, mer_id1, mer_id2):
        """Returns value of contact between two given mers according to
        three-valued logic.

        Arguments:
        mer_id1 -- reference to first monomer in contact which value is to
        be checked. Reference could be:
        monomer instance itself, its PyDesc ind, its PDB_id instance (or
        tuple containing proper values; see number converter docstring for
        more information) or monomer index on a list of structure mers.
        mer_id2 -- reference to second monomer corresponding to mer_id1.
        """
        return self._contacts[mer_id1, mer_id2]

    def to_string(self, stream_out):
        """Dumps pairs of mers in contacts to a given file-like object in
        CSV format."""
        converter = self.structure.converter
        with stream_out:
            for (ind1, ind2), value in self:
                pdb_id1 = converter.get_pdb_id(ind1)
                pdb_id2 = converter.get_pdb_id(ind2)
                line = f"{pdb_id1}\t{pdb_id2}\t{value}\n"
                stream_out.write(line)


class FrequencyContactMap:
    """Class representing maps of contact frequencies in trajectories or NMR
    structures.
    """

    def __init__(
        self,
        structures,
        contact_criterion_obj=None,
        ignore1=True,
        select1=Everything(),
        select2=Everything(),
    ):
        """Initialize ContactMap.

        Arguments:
        structures -- list of pydesc.structure.AbstractStructure subclass
        instances representing the same structure.
        contact_criterion_obj -- instance of
        pydesc.contacts.ContactCriterion determining how to calculate
        contacts. Initially set to None.
        If so, contact for residues is based on ca-cbx criterion,
        and contact for nucleotides is based on ion contact or ring center
        contact.
        ignore1 -- bool; determines if contact value 1 is to be treated as
        0 (True) or 2 (False).
        select1, select2 -- pydesc.selection.Selection subclass instances to be
        used in contact map calculation against each other.
        By default Everything selection is used.
        """
        self.contact_criterion = contact_criterion_obj
        self.selA = select1
        self.selB = select2
        self.stcA = select1.create_structure(structures[0])
        self.stcB = select2.create_structure(structures[0])
        self.contacts = None
        self.structures = structures
        self.ignore = ignore1

    def __iter__(self):
        return iter(
            [
                (list(self.stcA)[i].ind, list(self.stcB)[j].ind, v)
                for i, j, v in list(self.contacts.items())
            ]
        )

    @property
    def frames(self):
        """Returns number of trajectory frames."""
        return len(self.structures)

    # TODO: create FCM calculator and transfer this there
    def calculate_frequencies(self):
        def get_value(val):
            """Returns value depended on self.ignore."""
            if self.ignore:
                return int(val == 2)
            return int(bool(val))

        for stc in self.structures:
            cmap_calculator = ContactMapCalculator(
                stc, self.contact_criterion, self.selA, self.selB
            )
            cmap = cmap_calculator.calculate_contact_map()
            new_mtx = cmap_calculator.get_as_sparse_mtx(get_value)
            try:
                self.contacts += new_mtx
            except TypeError:
                self.contacts = new_mtx
