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

"""Base and auxiliary classes for contact criteria.

Contact criteria are used to calculate contact maps for structures.
Pre-defined, ready to use criteria are stored in `criteria.py` submodule.

"""

import re

import numpy
from scipy.sparse import dok_matrix

from pydesc.selection import Everything


class ContactCriterion:
    """Abstract criterion, base for all other criteria."""

    def __init__(self):
        self.selection1 = Everything()
        self.selection2 = Everything()
        self.asymmetric = False

    def set_selections(self, selection1, selection2):
        """Set selections determining for which mers contacts will be calculated.
        
        E.g. for criteria that only make sense for residues, but does not for
        nucleotides, one would like to add selections picking residues only.
        
        Note that this allows to define contacts between different types of mers,
        e.g. between nucleotides and residues.

        Args:
            selection1: selection instance. By default Everything selection is set.
            selection2: selection instance. By default Everything selection is set.

        """
        self.selection1 = selection1
        self.selection2 = selection2
        self.asymmetric = True

    def set_selection(self, selection):
        """Set single selection determining for which mers contacts will be
        calculated.

        Calls set_selection with *selection* argument twice.

        """
        self.set_selections(selection, selection)
        self.asymmetric = False

    def calculate_contacts(self, structure_obj):
        """Calculate all contacts in given structure.
        
        Contacts are usually sparse, so this method by default returns scipy dok
        matrix. Note that shape of returned matrix takes into account all mers from
        structure given structure was derived from. Size of returned matrix depends
        on number of inds recognised in structure's converter.
        
        In subclasses in most cases it should be sufficient to leave that method
        untouched and overwrite _fill_contact_matrix instead.

        Args:
            structure_obj: any pydesc

        Returns:
            scipy.sparse.dok_matrix: square matrix filled with contact values. Each
            row and column corresponds with mer ind, so for substructures returned
            matrix will have probably greater dimension than expected.

        """
        mers1 = self.selection1.create_structure(structure_obj)
        mers2 = self.selection2.create_structure(structure_obj)

        total_len = structure_obj.derived_from.converter.get_max_ind()
        contact_map = dok_matrix((total_len, total_len), dtype=numpy.uint8)
        contact_map = self._fill_contact_matrix(mers1, mers2, contact_map)

        return contact_map

    def calculate_inter_contacts(self, structure1, structure2):
        """

        Args:
            structure1:
            structure2:

        Returns:

        """
        mers1 = self.selection1.create_structure(structure1)
        mers2 = self.selection2.create_structure(structure2)

        if structure1.derived_from is not structure2.derived_from:
            raise ValueError(
                "Both given sub structures have to be part of the same " "structure"
            )
        total_len = structure1.derived_from.converter.get_max_ind()
        contact_map = dok_matrix((total_len, total_len), dtype=numpy.uint8)
        contact_map = self._fill_contact_matrix(mers1, mers2, contact_map)

        if self.asymmetric:
            mers1 = self.selection1.create_structure(structure2)
            mers2 = self.selection2.create_structure(structure1)

            contact_map = self._fill_contact_matrix(mers1, mers2, contact_map)

        return contact_map

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        """Fill given contact matrix for given two sets of mers.
        
        Meant to be overwritten in subclasses. ContactCriterion gives rather
        suboptimal, generic implementation, which requires implementation of
        _calculate_contact method. So there are two ways of writing subclasses:
        either overwrite _fill_contact_matrix, or leave it as is and implement
        _calculate_contact.
        
        This method suppose to fill given matrix with values so that element i,
        j store value of contact between i-th mer from mers1 and j-th mer from mers2,
        where i and j are corresponding mers inds (values of ind attribute).

        Args:
            mers1: sequence of mers.
            mers2: other sequence of mers (possibly the same).
            matrix: matrix of appropriate shape to be filled with calculated values.

        Returns:
            given matrix, filled with calculated values.

        """
        for mer1 in mers1:
            for mer2 in mers2:
                value = self._calculate_contact(mer1, mer2)
                matrix[mer1.ind, mer2.ind] = value
        return matrix

    def _calculate_contact(self, mer1, mer2):
        """Calculate value of single contact between given mers.
        
        Meant to be overwritten in subclasses that do not have optimal vectorised
        implementation (in which case they overwrite _fill_contact_matrix method).
        For all others -- raise NotImplemented error.

        Args:
            mer1: first mer instance.
            mer2: second mer instance.

        Returns:
            int: contact value (0 if there is no contact, 1 for uncertain contact
            or 2 for sure one).

        """
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

    def __str__(self):
        return "contact criterion"


class NotCriterion(ContactCriterion):
    """Class reverting criterion calculations.

    Args:
        criterion: contact criterion instance.

    """

    def __str__(self):
        criterion_repr = str(self.criterion)
        return f"not ({criterion_repr})"

    def __init__(self, criterion):
        self.criterion = criterion
        super().__init__()
        self.set_selections(criterion.selection1, criterion.selection2)

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        matrix[:, :] = 2
        contact_map = self.criterion.calculate_inter_contacts(mers1, mers2)
        matrix = dok_matrix(matrix + (-1 * contact_map), dtype=numpy.uint8)
        return matrix


class CombinedContact(ContactCriterion):
    """Abstract class, criteria obtained via logical operations on contact
    criteria.

    Args:
        : any number of other criteria.

    """

    def __init__(self, *criteria_objs):
        n_criteria = len(criteria_objs)
        if n_criteria < 2:
            raise AttributeError("Need at least two criteria to combine.")
        self.n_criteria = n_criteria
        self.criteria = criteria_objs
        super().__init__()

    def _repr_operation(self):
        """Returns regular expression operation to be used in __repr__ and __str__."""
        pattern = re.compile("[A-Z][a-z]*")
        name = self.__class__.__name__.replace("Contacts", "")
        operation = " ".join(re.findall(pattern, name)).lower()
        return operation

    def __repr__(self):
        self_repr = self._repr_operation().capitalize()
        return f"<{self_repr} of {self.n_criteria} contact criteria>"

    def __str__(self):
        self_repr = self._repr_operation()
        sub_criteria_repr = " and ".join(map(str, self.criteria))
        return f"{self_repr} of criteria based on {sub_criteria_repr}"


class ContactsConjunction(CombinedContact):
    """Conjunction of contact criteria."""

    def __init__(self, *criteria_objs):
        super().__init__(*criteria_objs)

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        contact_map = self.criteria[0].calculate_inter_contacts(mers1, mers2)
        for criterion in self.criteria[1:]:
            new_map = criterion.calculate_inter_contacts(mers1, mers2)
            contact_map = contact_map.minimum(new_map)
        return contact_map


class ContactsAlternative(CombinedContact):
    """Alternative of contact criteria."""

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        contact_map = self.criteria[0].calculate_inter_contacts(mers1, mers2)
        for criterion in self.criteria[1:]:
            new_map = criterion.calculate_inter_contacts(mers1, mers2)
            contact_map = contact_map.maximum(new_map)
        return contact_map


class ContactsExclusiveDisjunction(CombinedContact):
    """Exclusive disjunction (xor) of contact criteria."""

    def _fill_contact_matrix(self, mers1, mers2, matrix):
        new_map = self.criteria[0].calculate_inter_contacts(mers1, mers2)
        count_ones = (new_map == 1).astype(numpy.uint8)
        count_twos = (new_map == 2).astype(numpy.uint8)
        for criterion in self.criteria[1:]:
            new_map = criterion.calculate_inter_contacts(mers1, mers2)
            count_ones += (new_map == 1).astype(numpy.uint8)
            count_twos += (new_map == 2).astype(numpy.uint8)

        contact_map = dok_matrix(count_twos.shape, dtype=numpy.uint8)
        contact_map[(count_twos == 1)] = 2
        contact_map[count_ones != 0] = 1
        return contact_map
        # TODO: is that right?
