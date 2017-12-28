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

import pydesc.monomer
import pydesc.contacts.contacts
from pydesc.warnexcept import WrongMerType
from pydesc.selection import Everything

import numpy as np
import scipy.spatial


class ContactMapCalculator(object):

    def __init__(self, structure_obj, contact_criterion_obj=None, select1=Everything(), select2=Everything()):
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
        self.derived_from = structure_obj.derived_from
        self._contact_criterion = contact_criterion_obj
        self.structure = structure_obj
        self.frame = None
        self._contacts = None
        self.sel1 = select1.create_structure(structure_obj)
        self.sel2 = select2.create_structure(structure_obj)
        sel12 = select1 * select2
        self.sel12 = sel12.create_structure(structure_obj)
        self.sel1uni = (select1 - sel12).create_structure(structure_obj)
        self.sel2uni = (select2 - sel12).create_structure(structure_obj)
        self._rc_dist = None

    def __iter__(self):
        """Returns iterator thet runs over all contacts in contact map."""
        return iter([(i, j, v) for i, cs in self._contacts.items() for j, v in cs.items()])

    def calculate_rc_dist(self):
        items = (
            (self.sel12, self.sel12, '_rcdistC'),
            (self.sel1uni, self.sel12, '_rcdist1C'),
            (self.sel2uni, self.sel12, '_rcdist2C'),
            (self.sel1uni, self.sel2uni, '_rcdist12'),
            )
        for sel1, sel2, attrn in items:
            points1 = [i.rc.vector for i in sel1]
            points2 = [i.rc.vector for i in sel2]
            try:
                res = scipy.spatial.distance.cdist(points1, points2)
            except:
                res = None
            setattr(self, attrn, res)

    def calculate_contact_map(self):
        """Return ContactMap for structure set during initialization.

        Contact is a tuple containing first and second Monomer, distance(s) between them and a contact value under the ContactMap criterion.
        """

        def cmp(mer1, mer2):
            try:
                value = iic(mer1, mer2)
            except WrongMerType:
                return
            if value > 0:
                self._contacts[mer1.ind][mer2.ind] = value
                self._contacts[mer2.ind][mer1.ind] = value

        self.calculate_rc_dist()
        iic = self._contact_criterion.is_in_contact
        self._contacts = dict((monomer_obj.ind, {}) for monomer_obj in self.structure)
        max_rc_dist = getattr(self._contact_criterion, "max_rc_dist", None)

        if max_rc_dist is not None:
            #common(C) vs C
            indexes = np.transpose(np.where(self._rcdistC <= max_rc_dist))
            indexes = indexes[indexes[:, 0] < indexes[:, 1]]

            for (i, j) in indexes:
                mer1 = self.sel12._monomers[i]
                mer2 = self.sel12._monomers[j]
                cmp(mer1, mer2)

            items = (
                (self.sel1uni, self.sel12, '_rcdist1C'),    #uniA vs C
                (self.sel2uni, self.sel12, '_rcdist2C'),    #uniB vs C
                (self.sel1uni, self.sel2uni, '_rcdist12'),  #uniA vs uniB
                )
            for sel1, sel2, rcname in items:
                rcmtx = getattr(self, rcname)
                if not rcmtx:
                    continue
                indexes = np.transpose(np.where(rcmtx <= max_rc_dist))
                for (i, j) in indexes:
                    mer1 = sel1._monomers[i]
                    mer2 = sel2._monomers[j]
                    cmp(mer1, mer2)
        else:
            #common(C) vs C
            for i, mer1 in enumerate(self.sel12):
                for mer2 in self.sel12._monomers[i + 1:]:
                    cmp(mer1, mer2)
            items = (
                (self.sel1uni, self.sel12),     #uniA vs C
                (self.sel2uni, self.sel12),     #uniB vs C
                (self.sel1uni, self.sel2uni),   #uniA vs uniB
                )
            for sel1, sel2 in items:
                for mer1 in sel1:
                    for mer2 in sel2:
                        cmp(mer1, mer2)



class ContactMap(object):

    """Map of contacts present in a given (sub)structure."""

    def __init__(self, contacts, structure):
        """ContactMap constructor.

        Arguments:
        structure_obj -- instance of any pydesc.structure.AbstractStructure subclass for which contact map is to be created.
        contact_criterion_obj -- instance of pydesc.contacts.ContactCriterion determinion how to calculate contacts. Initailly set to None.
        If so, contact for residues is based on ca-cbx criterion, and contact for nucleotides is based on ion contact or ring center contact.
        select1, select2 -- pydesc.selection.Selection subclass instances to be used in contact map calculation against each other.
        By default Everything selection is used.
        """
        self.structure = structure
        self.frame = None
        self._contacts = contacts

    def __iter__(self):
        """Returns iterator that runs over all contacts in contact map."""
        return iter([(i, j, v) for i, cs in self._contacts.items() for j, v in cs.items()])

    def get_monomer_contacts(self, monomer_id, raw_numbering=False, internal=False):
        """Returns list of given monomer contacts.

        Contact is a tuple containing given monomer ind, ind of monomer that stays in contacts with it and
        contact value in three-valued logic. Distance evaluation is based on settings in configuration manager.

        Arguments:
        monomer_id -- monomer instance, PDB_id instance (or tuple containing proper values; see number converter
        docstring for more information), PyDesc ind or monomer index on a list of structure monomers.
        raw_numbering -- True or False. Idicates if monomer_id is a PyDesc ind (False) or a monomer list index (True).
        internal -- True or False. Indicates if method was called by other method.
        """
        ind = self._convert_to_ind(monomer_id, raw_numbering)
        if internal:
            return self.contacts[ind]
        return [(ind,) + i for i in self.contacts[ind].items()]

    def get_contact_value(self, monomer_id_1, monomer_id_2, raw_numbering=False):
        """Returns value of contact between two given mers according to three-valued logic.

        Arguments:
        monomer_id_1 -- reference to first monomer in concatc which value is to be checked. Reference could be:
        monomer instance itself, its PyDesc ind, its PDB_id instance (or tuple containing proper values; see number
        converter docstring for more information) or monomer index on a list of structure monomers.
        monomer_id_2 -- reference to second monomer corresponding to monomer_id_1.
        raw_numbering -- True or False. Idicates if given ids are PyDesc inds (False) or a monomer list indexes (True).
        """
        contacts_1 = self.get_monomer_contacts(monomer_id_1, raw_numbering, True)
        ind_2 = self._convert_to_ind(monomer_id_2, raw_numbering)
        try:
            return contacts_1[ind_2]
        except KeyError:
            return 0

    def get_contact_criterion(self, monomer_id_1, monomer_id_2, raw_numbering=False):
        """Returns criterion of contact between two given mers.

        If contact map criterion is an alternative of contacts - returns first matching subcriterion.
        If subcriterions are ContactsAlternative instance - their subcriterion are considered subcriterion of top level alternative.
        Otherwise contact map criterion is returned if given mers are in contact.
        None is returned if given mers are not in contact.

        Arguments:
        monomer_id_1 -- reference to first monomer in concatc which value is to be checked. Reference could be:
        monomer instance itself, its PyDesc ind, its PDB_id instance (or tuple containing proper values; see number
        converter docstring for more information) or monomer index on a list of structure monomers.
        monomer_id_2 -- reference to second monomer corresponding to monomer_id_1.
        raw_numbering -- True or False. Idicates if given ids are PyDesc inds (False) or a monomer list indexes (True).
        """
        if not self.get_contact_value(monomer_id_1, monomer_id_2, raw_numbering):
            return None
        mers = [self.structure[self._convert_to_ind(mer_id, raw_numbering)] for mer_id in (monomer_id_1, monomer_id_2)]
        try:
            return self.contact_criterion.get_validating_subcriterion(*mers)
        except AttributeError:
            return self.contact_criterion

    def _convert_to_ind(self, monomer_id, raw_numbering=False):
        """Returns PyDesc ind based on any given reference to monomer.

        Arguments:
        monomer_id -- reference that could be monomer instance itself, its PyDesc ind, its PDB_id instance (or tuple
        containing proper values; see numberconverter docstring for more information) or monomer index on a list
        of structure monomers.
        raw_numbering -- True or False. Indicates if given ids are PyDesc inds (False) or a monomer list indexes (True).
        """
        if raw_numbering:
            monomer_id = self.structure._mers[monomer_id]    # pylint: disable=protected-access

        if type(monomer_id) is int:
            return monomer_id
        try:
            return monomer_id.ind
        except AttributeError:
            return self.structure.derived_from.converter.get_ind(monomer_id)

    def dump(self, fobj):
        """Dumps pairs of mers in contacts to a given file-like object in CSV format."""
        with fobj:
            for k1 in self.contacts:
                m1 = self.structure[k1]
                for k2 in self.contacts[k1]:
                    m2 = self.structure[k2]
                    line = "%s\t%s\t%i\n" % (str(m1.pid), str(m2.pid), self.contacts[k1][k2])
                    fobj.write(line)

    @property
    def contact_criterion(self):
        """Property returning contact map contact criterion."""
        return self._contact_criterion

    @property
    def contacts(self):
        """Property returning contacts.

        Contacts are stored in nested dicts containing indexes of all mers as keys of 1st level dict,
        and indexes of all mers in contact with given mer as values of 2nd level dict.
        Contact values are values of 2nd level dicts. Mers that are not with contact with given mer
        have no entry in 2nd level dicts.
        """

        str_frame = self.substructure.frame

        if self.frame != str_frame:
            self.calculate_contacts()
            self.frame = str_frame

        return self._contacts
