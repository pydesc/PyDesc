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

import numpy as np
import scipy.spatial
from scipy.sparse import dok_matrix

from pydesc.contacts import contacts
from pydesc.selection import Everything
from pydesc.warnexcept import WrongMerType


class ContactMapCalculator:
    """Class responsible for calculating contact maps."""

    def __init__(
        self,
        structure_obj,
        contact_criterion_obj=None,
        select1=Everything(),
        select2=None,
    ):
        """ContactMapCalculator initializer.

        Arguments:
        structure_obj -- instance of any pydesc.structure.AbstractStructure
        subclass for which contact map is to be created.
        contact_criterion_obj -- instance of
        pydesc.contacts.ContactCriterion determining how to calculate contacts.
        Initially set to None.
        If so, contact for residues is based on ca-cbx criterion,
        and contact for nucleotides is based on ion contact or ring center
        contact.
        select1, select2 -- pydesc.selection.Selection subclass instances to be
        used in contact map calculation against each other.
        By default Everything selection is used.
        """
        if contact_criterion_obj is None:
            contact_criterion_obj = contacts.ContactsAlternative(
                contacts.CaCbxContact(),
                contacts.ContactsAlternative(
                    contacts.NIContact(), contacts.RingCenterContact()
                ),
            )
        self.contact_criterion = contact_criterion_obj
        self.structure = structure_obj
        self._rc_distC = None
        self._rc_dist1C = None
        self._rc_dist2C = None
        self._rc_dist12 = None

        if select2 is None:
            self.sel1 = tuple(select1.create_structure(structure_obj))
            self.sel2 = self.sel1
            self.sel12 = self.sel1
            self.sel1uni = ()
            self.sel2uni = ()
            return
        self.sel1 = tuple(select1.create_structure(structure_obj))
        self.sel2 = tuple(select2.create_structure(structure_obj))
        sel12 = select1 * select2
        self.sel12 = tuple(sel12.create_structure(structure_obj))
        self.sel1uni = tuple((select1 - sel12).create_structure(structure_obj))
        self.sel2uni = tuple((select2 - sel12).create_structure(structure_obj))

    def calculate_rc_dist(self):
        items = (
            (self.sel12, self.sel12, "_rc_distC"),
            (self.sel1uni, self.sel12, "_rc_dist1C"),
            (self.sel2uni, self.sel12, "_rc_dist2C"),
            (self.sel1uni, self.sel2uni, "_rc_dist12"),
        )
        for sel1, sel2, attr_name in items:
            points1 = np.array([i.rc.vector for i in sel1])
            points2 = np.array([i.rc.vector for i in sel2])
            try:
                res = scipy.spatial.distance.cdist(points1, points2)
            except ValueError:  # TODO: determine what exceptions
                res = None
            setattr(self, attr_name, res)

    def _get_dict_like_matrix(self):
        structure_length = self.structure.derived_from[-1].ind + 1
        dimension = (structure_length, structure_length)
        return dok_matrix(dimension, dtype=int)

    def _fill_interaction_matrix(self, interaction_function):
        """Return interaction matrix filled with values returned by given
        interaction function.

        Argument:
        interaction_function -- function taking mer1, mer2 and interaction
        matrix and returning None.
        Purpose of this function is to fill mtx with appropriate value in
        appropriate way.
        """
        self.calculate_rc_dist()
        max_rc_dist = getattr(self.contact_criterion, "max_rc_dist", None)
        contacts_mtx = self._get_dict_like_matrix()

        if max_rc_dist:
            # common(C) vs C
            try:
                indexes = np.transpose(np.where(self._rc_distC <= max_rc_dist))
                indexes = indexes[indexes[:, 0] < indexes[:, 1]]

                for (i, j) in indexes:
                    mer1 = self.sel12[i]
                    mer2 = self.sel12[j]
                    interaction_function(mer1, mer2, contacts_mtx)

            except (ValueError, IndexError):
                pass
        else:
            for i, mer1 in enumerate(self.sel12):
                for mer2 in self.sel12[i + 1 :]:
                    interaction_function(mer1, mer2, contacts_mtx)

        items = (
            (self.sel1uni, self.sel12, "_rc_dist1C"),  # uniA vs C
            (self.sel2uni, self.sel12, "_rc_dist2C"),  # uniB vs C
            (self.sel1uni, self.sel2uni, "_rc_dist12"),  # uniA vs uniB
        )
        for mer_tuple_1, mer_tuple_2, rc_name in items:
            if max_rc_dist:
                try:
                    rc_mtx = getattr(self, rc_name)
                    indexes = np.transpose(np.where(rc_mtx <= max_rc_dist))
                    for (i, j) in indexes:
                        mer1 = mer_tuple_1[i]
                        mer2 = mer_tuple_2[j]
                        interaction_function(mer1, mer2, contacts_mtx)
                except (ValueError, TypeError):
                    pass
            else:
                for mer1 in mer_tuple_1:
                    for mer2 in mer_tuple_2:
                        interaction_function(mer1, mer2, contacts_mtx)

        return contacts_mtx

    def calculate_contact_map(self):
        """Return ContactMap for structure set during initialization.

        Contact is a tuple containing first and second Monomer, distance(s)
        between them and a contact value under the ContactMap criterion.
        """

        def compare_mers(mer_1, mer_2, mtx):
            try:
                value = self.contact_criterion.is_in_contact_no_pre_check(mer_1, mer_2)
            except WrongMerType:
                return
            if value > 0:
                mtx[mer_1.ind, mer_2.ind] = value
                mtx[mer_2.ind, mer_1.ind] = value

        contacts_mtx = self._fill_interaction_matrix(compare_mers)

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
        with stream_out:
            for (k1, k2), value in self._contacts:
                m1 = self.structure[k1]
                m2 = self.structure[k2]
                line = "%s\t%s\t%i\n" % (str(m1.pid), str(m2.pid), value)
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
