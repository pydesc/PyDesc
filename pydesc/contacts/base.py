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
Base class for contact criteria, combined criteria and utility functions
that deal helps in dealing with contacts among mers present in PyDesc
(sub)structures.

created: 20.05.2019 - , Tymoteusz 'hert' Oleniecki
"""

import re
from abc import ABCMeta
from abc import abstractmethod

import numpy

from pydesc.warnexcept import WrongMerType
from pydesc.selection import Everything


class ContactCriterion(metaclass=ABCMeta):
    """Abstract class, criteria instances."""

    def __init__(self):
        self.selection1 = Everything()
        self.selection2 = Everything()

    def set_selections(self, selection1, selection2):
        self.selection1 = selection1
        self.selection2 = selection2

    def calculate_contacts(self, structure_obj):

        mers1 = self.selection1.create_structure(structure_obj)
        mers2 = self.selection2.create_structure(structure_obj)

        contact_map = numpy.zeros((len(mers1), len(mers2)), dtype=numpy.uint8)
        contact_map = self._fill_contact_matrix(mers1, mers2, contact_map)

        return contact_map

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        for mer1 in mers1:
            for mer2 in mers2:
                value = self._calculate_contact(mer1, mer2)
                matrix[mer1.ind, mer2.ind] = value
        return matrix

    def _calculate_contact(self, mer1, mer2):
        raise NotImplemented

    def __or__(self, other):
        """Returns ContactsAlternative of self and other contact criterion"""
        return ContactsAlternative(self, other)

    def __and__(self, other):
        """Returns ContactsConjunction of self and other contact criterion"""
        return ContactsConjunction(self, other)

    def __xor__(self, other):
        """Returns ContactsExclusiveDisjunction of self and other contact
        criterion"""
        return ContactsExclusiveDisjunction(self, other)

    @abstractmethod
    def __str__(self):
        return "contact criterion"


class NotCriterion(ContactCriterion):
    def __str__(self):
        criterion_repr = str(self.criterion)
        return f"not ({criterion_repr})"

    def __init__(self, criterion):
        self.criterion = criterion
        super().__init__()
        self.set_selections(criterion.selection1, criterion.selection2)

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        method = self.criterion._fill_contact_matrix
        method(mers1, mers2, matrix)
        return (matrix - 1) * -1 + 1  # revert values


class CombinedContact(ContactCriterion, metaclass=ABCMeta):
    """Abstract class, criteria obtained via logical operations on contact
    criteria."""

    def __init__(self, *criteria_objs):
        """Combined criteria constructor.

        Arguments:
        criteria_objs -- basic criteria objects.
        """
        n_criteria = len(criteria_objs)
        if n_criteria < 2:
            raise AttributeError("Need at least two criteria to combine.")
        self.n_criteria = n_criteria
        self.criteria = criteria_objs
        super().__init__()

    def _repr_operation(self):
        """Returns regular expression operation to be used in __repr__ and
        __str__."""
        pattern = re.compile("[A-Z][a-z]*")
        name = self.__class__.__name__.replace("Contacts", "")
        operation = " ".join(re.findall(pattern, name)).lower()
        return operation

    def __repr__(self):
        self_repr = self._repr_operation().capitalize()
        sub_criteria_repr = " and ".join(map(str, self.criteria))
        return f"<{self_repr} of criteria based on {sub_criteria_repr}>"

    def __str__(self):
        self_repr = self._repr_operation()
        sub_criteria_repr = " and ".join(map(str, self.criteria))
        return f"{self_repr} of criteria based on {sub_criteria_repr}"

    @abstractmethod
    def _fill_contact_matrix(self, mers1, mers2, matrix):
        method = self.criteria[0]._fill_contact_matrix
        matrix = method(mers1, mers2, matrix)
        for criterion in self.criteria[1:]:
            method = criterion._fill_contact_matrix
            new_layer = method(mers1, mers2, matrix)
            matrix = numpy.stack((matrix, new_layer), axis=2)
        return matrix


class ContactsConjunction(CombinedContact):
    def _fill_contact_matrix(self, mers1, mers2, matrix):
        contact_map = super()._fill_contact_matrix(mers1, mers2, matrix)
        return numpy.min(contact_map, axis=2)


class ContactsAlternative(CombinedContact):
    def _fill_contact_matrix(self, mers1, mers2, matrix):
        contact_map = super()._fill_contact_matrix(mers1, mers2, matrix)
        return numpy.max(contact_map, axis=2)

    def get_validating_sub_criterion(self, mer_1, mer_2):
        for contact_criterion in self.criteria:
            try:
                return contact_criterion.get_validating_sub_criterion(mer_1, mer_2)
            except AttributeError:
                try:
                    if contact_criterion.is_in_contact(mer_1, mer_2):
                        return contact_criterion
                except WrongMerType:
                    continue
            except ValueError:
                continue
        raise ValueError("Given mers are not in contact.")


class ContactsExclusiveDisjunction(CombinedContact):
    def _fill_contact_matrix(self, mers1, mers2, matrix):
        contact_map = super()._fill_contact_matrix(mers1, mers2, matrix)
        # fields where only one value is sure are 2
        twos = numpy.count_nonzero(contact_map == 2, axis=2)
        result = (twos == 1) * 2
        # only fields where 2 occurred less than once are taken into account further
        mask = twos == 0
        undecided = numpy.max(contact_map == 1, axis=2) * mask
        result = numpy.stack((result, undecided.astype(bool)), axis=2)
        # 2 is picked over 1, 1 over 0.
        return numpy.max(result, axis=2)
