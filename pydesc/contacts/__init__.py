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

import pydesc.mers
import pydesc.contacts.contacts
from pydesc.warnexcept import WrongMerType
from pydesc.selection import Everything


class ContactMapCalculator(object):

    def __init__(self, structure_obj, contact_criterion_obj=None, select1=Everything(), select2=None):
        """ContactMapCalculator constructor.

        Arguments:
        structure_obj -- instance of any pydesc.structure.AbstractStructure subclass for which contact map is to be
        created.
        contact_criterion_obj -- instance of pydesc.contacts.ContactCriterion determinion how to calculate contacts.
        Initially set to None.
        If so, contact for residues is based on ca-cbx criterion, and contact for nucleotides is based on ion contact or
        ring center contact.
        select1, select2 -- pydesc.selection.Selection subclass instances to be used in contact map calculation against
        each other.
        By default Everything selection is used.
        """
        if contact_criterion_obj is None:
            contact_criterion_obj = pydesc.contacts.contacts.ContactsAlternative(
                pydesc.contacts.contacts.CaCbxContact(),
                pydesc.contacts.contacts.ContactsAlternative(
                    pydesc.contacts.contacts.NIContact(),
                    pydesc.contacts.contacts.RingCenterContact()
                )
            )
        self.contact_criterion = contact_criterion_obj
        self.structure = structure_obj
        structure_length = structure_obj.derived_from[-1].ind
        self._contacts = dok_matrix((structure_length + 1, structure_length + 1), dtype=int)
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

    def __iter__(self):
        """Returns iterator that runs over all contacts in contact map."""
        return iter([(i, j, v) for i, cs in self._contacts.items() for j, v in cs.items()])

    def calculate_rc_dist(self):
        items = (
            (self.sel12, self.sel12, '_rc_distC'),
            (self.sel1uni, self.sel12, '_rc_dist1C'),
            (self.sel2uni, self.sel12, '_rc_dist2C'),
            (self.sel1uni, self.sel2uni, '_rc_dist12'),
        )
        for sel1, sel2, attr_name in items:
            points1 = np.array([i.rc.vector for i in sel1])
            points2 = np.array([i.rc.vector for i in sel2])
            try:
                res = scipy.spatial.distance.cdist(points1, points2)
            except Exception as e:  # TODO: determine what exceptions
                res = None
            setattr(self, attr_name, res)

    def calculate_contact_map(self):
        """Return ContactMap for structure set during initialization.

        Contact is a tuple containing first and second Monomer, distance(s) between them and a contact value under the ContactMap criterion.
        """

        def compare_mers(mer_1, mer_2):
            try:
                value = self.contact_criterion.is_in_contact_no_pre_check(mer_1, mer_2)
            except WrongMerType:
                return
            if value > 0:
                self._contacts[mer_1.ind, mer_2.ind] = value
                self._contacts[mer_2.ind, mer_1.ind] = value

        self.calculate_rc_dist()
        max_rc_dist = getattr(self.contact_criterion, "max_rc_dist", None)

        if max_rc_dist:
            # common(C) vs C
            try:
                indexes = np.transpose(np.where(self._rc_distC <= max_rc_dist))
                indexes = indexes[indexes[:, 0] < indexes[:, 1]]

                for (i, j) in indexes:
                    mer1 = self.sel12[i]
                    mer2 = self.sel12[j]
                    compare_mers(mer1, mer2)

            except (ValueError, IndexError):
                pass
        else:
            for i, mer1 in enumerate(self.sel12):
                for mer2 in self.sel12[i + 1:]:
                    compare_mers(mer1, mer2)

        items = (
            (self.sel1uni, self.sel12, '_rc_dist1C'),  # uniA vs C
            (self.sel2uni, self.sel12, '_rc_dist2C'),  # uniB vs C
            (self.sel1uni, self.sel2uni, '_rc_dist12'),  # uniA vs uniB
        )
        for mer_tuple_1, mer_tuple_2, rc_name in items:
            if max_rc_dist:
                try:
                    rc_mtx = getattr(self, rc_name)
                    indexes = np.transpose(np.where(rc_mtx <= max_rc_dist))
                    for (i, j) in indexes:
                        mer1 = mer_tuple_1[i]
                        mer2 = mer_tuple_2[j]
                        compare_mers(mer1, mer2)
                except ValueError:
                    pass
            else:
                for mer1 in mer_tuple_1:
                    for mer2 in mer_tuple_2:
                        compare_mers(mer1, mer2)

        return ContactMap(self._contacts, self.structure)


class ContactMap(object):
    """Map of contacts present in a given (sub)structure."""

    def __init__(self, contacts_mtx, structure):
        """ContactMap constructor.

        Arguments:
        contacts_mtx -- sparse dict-like contact matrix representing contacts in biopolymer. 2 indicates contact,
        1 -- plausible contact.
        structure -- instance of any pydesc.structure.AbstractStructure subclass for which contact map is to be created.
        """
        self.structure = structure
        self.frame = None
        self._contacts = contacts_mtx

    def __iter__(self):
        """Returns iterator that runs over all contacts in contact map."""
        return iter(self._contacts.items())

    def __len__(self):
        """Return length of self._contacts dok_matrix keys list."""
        return len(self._contacts)

    def __getitem__(self, item):
        try:
            return self.get_contact_value(*item)
        except TypeError:
            return self.get_monomer_contacts(item)

    def get_monomer_contacts(self, monomer_id, raw_numbering=False):
        """Returns list of given monomer contacts.

        Contact is a tuple containing given monomer ind, ind of monomer that stays in contacts with it and
        contact value in three-valued logic. Distance evaluation is based on settings in configuration manager.

        Arguments:
        monomer_id -- monomer instance, PDB_id instance (or tuple containing proper values; see number converter
        docstring for more information), PyDesc ind or monomer index on a list of structure mers.
        raw_numbering -- True or False. Idicates if monomer_id is a PyDesc ind (False) or a monomer list index (True).
        """
        ind = self._convert_to_ind(monomer_id, raw_numbering)
        return self._contacts[ind].items()

    def get_contact_value(self, monomer_id_1, monomer_id_2, raw_numbering=False):
        """Returns value of contact between two given mers according to three-valued logic.

        Arguments:
        monomer_id_1 -- reference to first monomer in concatc which value is to be checked. Reference could be:
        monomer instance itself, its PyDesc ind, its PDB_id instance (or tuple containing proper values; see number
        converter docstring for more information) or monomer index on a list of structure mers.
        monomer_id_2 -- reference to second monomer corresponding to monomer_id_1.
        raw_numbering -- True or False. Idicates if given ids are PyDesc inds (False) or a monomer list indexes (True).
        """
        ind1 = self._convert_to_ind(monomer_id_1, raw_numbering=raw_numbering)
        ind2 = self._convert_to_ind(monomer_id_2, raw_numbering=raw_numbering)
        return self._contacts[ind1, ind2]

    def _convert_to_ind(self, monomer_id, raw_numbering=False):
        """Returns PyDesc ind based on any given reference to monomer.

        Arguments:
        monomer_id -- reference that could be monomer instance itself, its PyDesc ind, its PDB_id instance (or tuple
        containing proper values; see numberconverter docstring for more information) or monomer index on a list
        of structure mers.
        raw_numbering -- True or False. Indicates if given ids are PyDesc inds (False) or a monomer list indexes (True).
        """
        if raw_numbering:
            monomer_id = tuple(self.structure)[monomer_id]

        try:
            return monomer_id.ind
        except AttributeError:
            try:
                return self.structure.derived_from.converter.get_ind(monomer_id)
            except Exception as e:  # TODO: what error is raised?
                return monomer_id  # should be int in that case

    def to_string(self, stream_out):
        """Dumps pairs of mers in contacts to a given file-like object in CSV format."""
        with stream_out:
            for (k1, k2), value in self._contacts:
                m1 = self.structure[k1]
                m2 = self.structure[k2]
                line = "%s\t%s\t%i\n" % (str(m1.pid), str(m2.pid), value)
                stream_out.write(line)


class FrequencyContactMap(object):
    """Class representing maps of contact frequencies in trajectories or NMR structures."""

    def __init__(self, structures, contact_criterion_obj=None, ignore1=True, select1=Everything(),
                 select2=Everything()):
        """ContactMap costructor.

        Arguments:
        structures -- list of pydesc.structure.AbstractStructure subclass instances representing the same structure.
        contact_criterion_obj -- instance of pydesc.contacts.ContactCriterion determinion how to calculate contacts. Initailly set to None.
        If so, contact for residues is based on ca-cbx criterion, and contact for nucleotides is based on ion contact or ring center contact.
        ignore1 -- bool; determines if contact value 1 is to be treted as 0 (True) or 2 (False).
        select1, select2 -- pydesc.selection.Selection subclass instances to be used in contact map calculation against each other.
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
        return iter([(list(self.stcA)[i].ind, list(self.stcB)[j].ind, v) for i, j, v in self.contacts.items()])

    @property
    def frames(self):
        """Returns number of trajectory frames."""
        return len(self.structures)

    def calculate_frequencies(self):
        def get_value(val):
            """Returns value depended on self.ignore."""
            if self.ignore:
                return int(val == 2)
            return int(bool(val))

        for stc in self.structures:
            cmap = ContactMap(stc, self.contact_criterion, self.selA, self.selB)
            cmap.calculate_contacts()
            new_mtx = cmap.get_as_sparse_mtx(get_value)
            try:
                self.contacts += new_mtx
            except TypeError:
                self.contacts = new_mtx
