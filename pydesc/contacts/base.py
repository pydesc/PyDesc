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
from copy import deepcopy
from functools import wraps

import pydesc.mers
from pydesc.warnexcept import WrongMerType


def for_monomer_type_only(type_1, type_2=None):
    """Class decorator used to assert correct type of mers for subsequent
    contact evaluation.

    This decorator sets type_1 and type_2 class attributes to provided
    values. It is implemented mostly for backward
    compatibility.

    Arguments:
    type_1 -- class of monomer for first mer.
    type_2 -- class of monomer for second mer; initially set to None,
    if so type_1 is taken as type_2.
    """

    def proper_type_decorator(criterion_class):
        """Decorator returning a criterion class with type_1 and type_2
        class attributes set to provided values."""

        criterion_class.set_types_cls(type_1, type_2)
        criterion_class.is_in_contact_no_pre_check = check_type(
            criterion_class.is_in_contact_no_pre_check
        )
        try:
            criterion_class.calculate_distance = check_type(
                criterion_class.calculate_distance
            )
        except AttributeError:
            pass

        return criterion_class

    return proper_type_decorator


def ignore_exceptions(*exceptions):
    """Class decorator that allows to ignore given errors raised by criteria.

    Argument:
    exceptions -- type of errors that are to be ignored.

    Decorator wraps _is_in_contact method that raises CannotCalculateContact
    error and returns 0 instead.
    """

    def decorator(criterion_class):
        original_is_in_contact = criterion_class._is_in_contact

        @wraps(original_is_in_contact)
        def wrapped_is_in_contact(self, *args, **kwargs):
            """Wrapped is_in_contact method that returns 0 instead of
            raising CannotCalculateContact error."""
            try:
                return original_is_in_contact(self, *args, **kwargs)
            except exceptions:
                return 0

        criterion_class._is_in_contact = wrapped_is_in_contact
        return criterion_class

    return decorator


def not_decorator(criterion_class):
    """Class decorator that changes _is_in_contact method.

    Argument:
    criterion_class -- instance of ContactCriterion class.

    Wrapped _is_in_contact method returns opposite values than original
    method: 0 for 2, 1 for 1 and 2 for 0.
    """

    original_is_in_contact = criterion_class._is_in_contact

    @wraps(original_is_in_contact)
    def wrapped_is_in_contact(self, *args, **kwargs):
        """Wrapped _is_in_contact method that returns 0 instead of raising
        CannotCalculateContact error."""
        res = original_is_in_contact(self, *args, **kwargs)
        if res == 0:
            return 2
        elif res == 2:
            return 0
        # when res == 1
        return res

    criterion_class._is_in_contact = wrapped_is_in_contact

    return criterion_class


def check_type(ory_mth):
    """
    """

    @wraps(ory_mth)
    def new_mth(self, monomer_1, monomer_2, *args, **kwargs):
        if self.type_2 is None or self.type_1 == self.type_2:
            if self.type_1 is None:
                return ory_mth(self, monomer_1, monomer_2, *args, **kwargs)
            else:
                if self._test_type(monomer_1, self.type_1) and self._test_type(
                    monomer_2, self.type_1
                ):
                    return ory_mth(self, monomer_1, monomer_2, *args, **kwargs)
        else:
            if self._test_type(monomer_1, self.type_1) and self._test_type(
                monomer_2, self.type_2
            ):
                return ory_mth(self, monomer_1, monomer_2, *args, **kwargs)
            elif self._test_type(monomer_1, self.type_2) and self._test_type(
                monomer_2, self.type_1
            ):
                return ory_mth(self, monomer_2, monomer_1, *args, **kwargs)

        msg_tup = (monomer_1, monomer_2, self.type_1, self.type_2, self.__class__)
        raise WrongMerType(*msg_tup)

    addition = "\n\nWorks only for certain type of mers."
    new_mth.__doc__ = ory_mth.__doc__ + addition

    return new_mth


class ContactCriterion(metaclass=ABCMeta):
    """Abstract class, criteria instances."""

    type_1 = None
    type_2 = None

    _test_type_cache = {}

    @property
    def criteria(self):
        """Returns list containing current object."""
        return [self]

    def set_types(self, type_1, type_2=None):
        """Sets types of mers for which criterion is to be applied.

        Arguments:
        type_1 -- pydesc.monomer.Monomer subclass.
        type_2 -- if not given, type_1 assumed for both.
        """
        if type_1 == type_2:
            type_2 = None

        if type_1 == pydesc.mers.Mer and type_2 is None:
            type_1 = None

        self.type_1 = type_1
        self.type_2 = type_2

    def get_types(self):
        """Returns types of mer for which criterion is created."""
        type_1 = pydesc.mers.Mer if self.type_1 is None else self.type_1
        type_2 = type_1 if self.type_2 is None else self.type_2

        return type_1, type_2

    set_types_cls = classmethod(set_types)

    def is_in_contact(self, monomer_1, monomer_2, *args, **kwargs):
        """Returns three-valued logic contact value.

        This method checks monomer types and calls
        is_in_contact_no_pre_check which does actual job.

        Arguments:
        monomer_1 -- first monomer instance.
        monomer_2 -- second mers instance.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more
        information.
        """

        if not self._pre_check(monomer_1, monomer_2, **kwargs):
            return 0

        return self.is_in_contact_no_pre_check(monomer_1, monomer_2, *args, **kwargs)

    def is_in_contact_no_pre_check(self, monomer_1, monomer_2, *args, **kwargs):
        """Like is_in_contact, but without pre-check."""
        return self._is_in_contact(monomer_1, monomer_2, *args, **kwargs)

    @staticmethod
    def _test_type(monomer, mtype):
        """
        Checks if monomer matches a given type.

        This is a wrapper for the isinstance built-in which caches most
        frequent queries.

        Arguments:
            monomer -- monomer
            mtype -- type
        """

        if type(monomer) == mtype:
            return True

        try:
            (good, bad) = ContactCriterion._test_type_cache[mtype]

            if type(monomer) in good:
                return True
            elif type(monomer) in bad:
                return False
        except KeyError:
            pass

        if isinstance(monomer, mtype):
            try:
                ContactCriterion._test_type_cache[mtype][0].append(type(monomer))
            except KeyError:
                ContactCriterion._test_type_cache[mtype] = ([type(monomer)], [])
            return True

        try:
            ContactCriterion._test_type_cache[mtype][1].append(type(monomer))
        except KeyError:
            ContactCriterion._test_type_cache[mtype] = ([], [type(monomer)])

        return False

    def _pre_check(self, monomer_1, monomer_2, rc_dist=None, **kwargs):
        """A method for checking a quick and easy precondition of a contact.

        This method should be overridden whenever possible to provide a
        quick and dirty checking, before computing actual criteria.

        This method accepts arguments of any type and returns a boolean.

        Ideally a condition tested here should be simpler even than type
        checking.
        """

        try:
            if rc_dist is not None:
                return rc_dist <= self.max_rc_dist
            else:
                return abs(monomer_1.rc - monomer_2.rc) <= self.max_rc_dist
        except:  # If anything goes wrong, just disregard the whole precheck.
            pass

        return True

    @abstractmethod
    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, lazy=True):
        """Abstract method overridden in subclasses.

        Returns three-valued logic contact value.

        This method is called by is_in_contact, which is supposed to check
        monomer types.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more
        information.
        """
        pass

    def __eq__(self, criterion_obj):
        """Checks if given objects are equal.

        Argument:
        criterion_obj -- object to be compared with current object.
        """
        if type(self) != type(criterion_obj):
            return False
        if self.__dict__ != criterion_obj.__dict__:
            return False
        return True

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

    def __invert__(self):
        new_class = not_decorator(type(self))
        copy = deepcopy(self)
        copy.__class__ = new_class
        return copy


class CombinedContact(ContactCriterion, metaclass=ABCMeta):
    """Abstract class, criteria obtained via logical operations on contact
    criteria."""

    def __init__(self, *criteria_objs):
        """Combined criteria constructor.

        Arguments:
        criteria_objs -- basic criteria objects.
        """
        self._criteria = criteria_objs
        if len(self.criteria) < 2:
            raise AttributeError("Need at least two criteria to combine.")

    def _repr_operation(self):
        """Returns regular expression operation to be used in __repr__ and
        __str__."""
        pattern = re.compile("[A-Z][a-z]*")
        name = self.__class__.__name__.replace("Contacts", "")
        operation = " ".join(re.findall(pattern, name)).lower()
        return operation

    def __repr__(self):
        return "<%s of criteria based on %s>" % (
            self._repr_operation().capitalize(),
            " and ".join(map(str, self.criteria)),
        )

    def __str__(self):
        return "%s of %s criteria" % (
            self._repr_operation(),
            " and ".join(map(str, self.criteria)),
        )

    def is_in_contact(self, monomer_1, monomer_2, *args, **kwargs):
        return self._is_in_contact(monomer_1, monomer_2, *args, no_pre_check=True)

    def is_in_contact_no_pre_check(self, monomer_1, monomer_2, *args, **kwargs):
        return self._is_in_contact(monomer_1, monomer_2, *args)

    @abstractmethod
    def _is_in_contact(
        self, monomer_1_obj, monomer_2_obj, lazy=True, no_pre_check=False, **kwargs
    ):
        """Returns value of combined contact criterion.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- True or False, initially set to True. Determines if
        sub-criteria values are to be calculated lazy or not.
        Lazy calculation means that if program is able to assume that
        criterion is not satisfied during calculation sub-criterion -
        further sub-criteria are not calculated.
        """
        pass

    def _calculate_distance(self, mer1, mer2, *args, **kwargs):
        raise NotImplementedError("Combined contacts cannot calc distances.")

    @property
    def criteria(self):
        """Returns sequence of combined criterion criteria."""
        return self._criteria


class ContactsConjunction(CombinedContact):
    """Conjunction of criteria given in a list.

    Computes criteria type as an intersection of types accepted by
    sub-criteria. Algorithm used to resolve types may fail if multiple
    inheritance is used. Given criteria could be a CombinedContact instance.
    """

    def __init__(self, *criteria_objs):
        """Conjunction of criteria constructor.

        Arguments:
        criteria_objs -- any number of basic or combined criteria objects.
        """
        CombinedContact.__init__(self, *criteria_objs)

        type_1 = pydesc.mers.Mer
        type_2 = pydesc.mers.Mer

        for i in criteria_objs:
            (t1, t2) = i.get_types()

            if issubclass(t1, type_1):
                type_1 = t1
            elif not issubclass(type_1, t1):
                raise AttributeError("Given criteria require incompatible mer types.")

            if issubclass(t2, type_2):
                type_2 = t2
            elif not issubclass(type_2, t2):
                raise AttributeError("Given criteria require incompatible mer types.")

            try:
                d = i.max_rc_dist
                self.max_rc_dist = min(getattr(self, "max_rc_dist", d), d)
            except AttributeError:
                pass

        self.set_types(type_1, type_2)

    def _is_in_contact(
        self, monomer_1_obj, monomer_2_obj, lazy=True, no_pre_check=False, **kwargs
    ):
        """Returns contact value under conjunction of the given criteria.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- True or False, initially set to True. Determines if
        sub-criteria values are to be calculated lazy or not.
        Lazy calculation means that if program is able to assume that
        criterion is not satisfied during calculation sub-criterion - further
        sub-criteria are not calculated.
        """
        values = []
        for contact_criterion in self.criteria:
            try:
                value = contact_criterion.is_in_contact(
                    monomer_1_obj, monomer_2_obj, lazy=lazy, **kwargs
                )
            except WrongMerType:
                value = 0
            values.append(value)
            if lazy and value == 0:
                break
        if all(value == 2 for value in values):
            return 2
        elif all(value >= 1 for value in values):
            return 1
        else:
            return 0


class ContactsDisjunction(CombinedContact):
    """
    Abstract class grouping combined contacts which require only some
    criteria to be satisfied.

    Computes criteria type to meet any type of sub-criteria.
    """

    def __init__(self, *criteria_objs):
        """Conjunction of criteria constructor.

        Arguments:
        criteria_objs -- basic criteria objects.
        """
        CombinedContact.__init__(self, *criteria_objs)

        def lcs(types_):
            """ Lowest common superclass"""
            if len(types_) == 0:
                return None

            mros = [x.mro() for x in types_]
            for x in mros[0]:
                if all(x in mro for mro in mros):
                    return x

        types = [crit.get_types() for crit in criteria_objs]

        type_1 = lcs([t[0] for t in types])
        type_2 = lcs([t[1] for t in types])

        self.set_types(type_1, type_2)

        self.max_rc_dist = 0
        for i in self.criteria:
            try:
                self.max_rc_dist = max(self.max_rc_dist, i.max_rc_dist)
            except AttributeError:
                del self.max_rc_dist
                break


class ContactsAlternative(ContactsDisjunction):
    """Alternative of criteria given in the list.

    Given criteria could be a CombinedContact instance.
    """

    def _is_in_contact(
        self, monomer_1_obj, monomer_2_obj, lazy=True, no_pre_check=False, **kwargs
    ):
        """Returns contact value under an alternative of given criteria.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- True or False, initially set to True. Determines if
        sub-criteria values are to be calculated lazy or not.
        Lazy calculation means that if program is able to assume that
        criterion is not satisfied during calculation sub-criterion -
        further sub-criteria are not calculated.
        """
        attr = "is_in_contact_no_pre_check" if no_pre_check else "is_in_contact"
        values = []
        for contact_criterion in self.criteria:
            try:
                method = getattr(contact_criterion, attr)
                value = method(monomer_1_obj, monomer_2_obj, lazy=lazy, **kwargs)
            except WrongMerType:
                value = 0
            values.append(value)
            if lazy and value == 2:
                break
        if any(value == 2 for value in values):
            return 2
        elif any(value == 1 for value in values):
            return 1
        else:
            return 0

    def get_validating_sub_criterion(self, mer_1, mer_2):
        """Returns sub-criterion for which given mers are in contact.

        Arguments:
        mer_1, mer_2 -- pydesc.monomer.Monomer instances.

        Raises ValueError in given mers are not in contact.
        """
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


class ContactsExclusiveDisjunction(ContactsDisjunction):
    """Exclusive Disjunction of criteria given in the list.

    Given criteria could be a CombinedContact instance.
    """

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, lazy=True, **kwargs):
        """Returns contact value under an exclusive disjunction of given
        criteria.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- True or False, initially set to True. Determines if
        sub-criteria values are to be calculated lazy or not.
        Lazy calculation means that if program is able to assume that
        criterion is not satisfied during calculation sub-criterion -
        further sub-criteria are not calculated.
        """
        values = []
        for contact_criterion in self.criteria:
            try:
                value = contact_criterion.is_in_contact(
                    monomer_1_obj, monomer_2_obj, lazy=lazy, **kwargs
                )
            except WrongMerType:
                value = 0
            values.append(value)
            if lazy and values.count(2) > 1:
                break
        if values.count(2) == 1 and not any(value == 1 for value in values):
            return 2
        elif values.count(2) in (0, 1) and any(value == 1 for value in values):
            return 1
        else:
            return 0


class DescriptorCriterion(ContactCriterion):
    """Contacts present in a given descriptor.

    This is a helper class useful to extract topology of a descriptor.
    """

    def __init__(self, descriptor_obj):
        """Descriptor criteria constructor.

        Argument:
        descriptor_obj -- instance of PyDesc descriptor.
        """
        self.desc = descriptor_obj

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns value of contact if given mers have their Contact
        instance in criterion desc attr.

        Returns contact value if both given mers are central mers for
        elements in any pydesc.structure.Contact stored in current criterion
        desc.contacts. Otherwise returns zero.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- always set to None.
        """
        for contact in self.desc.contacts:
            if not max(contact.elements).central_monomer in [
                monomer_1_obj,
                monomer_2_obj,
            ]:
                continue

            if not min(contact.elements).central_monomer in [
                monomer_1_obj,
                monomer_2_obj,
            ]:
                continue

            return contact.value

        return 0
