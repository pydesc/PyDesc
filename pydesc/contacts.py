# Copyright 2017 Agnieszka Mykowiecka, Tymoteusz Oleniecki
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

created: 18.07.2013 - , Tymoteusz 'hert' Oleniecki, Agnieszka Mykowiecka
"""

import pydesc.monomer
from pydesc.config import ConfigManager
from pydesc.warnexcept import WrongMonomerType
from pydesc.warnexcept import CannotCalculateContact

import scipy.spatial
import numpy
import re
from abc import ABCMeta, abstractmethod

# pylint: disable=no-member
ConfigManager.new_branch("contacts")
ConfigManager.contacts.set_default("ca_contact_distance", 6.0)
ConfigManager.contacts.set_default("ca_contact_undecidable_range", 0.5)
ConfigManager.contacts.set_default("cbx_contact_distance", 6.5)
ConfigManager.contacts.set_default("cbx_contact_undecidable_range", 0.5)
ConfigManager.contacts.set_default("rc_contact_distance", 7.5)
ConfigManager.contacts.set_default("rc_contact_undecidable_range", 0.5)
ConfigManager.contacts.set_default("ring_center_contact_distance", 6.25)
ConfigManager.contacts.set_default("ring_center_contact_undecidable_range", 0.0)
ConfigManager.contacts.set_default("at_contact_distance", 5.0)
ConfigManager.contacts.set_default("at_contact_undecidable_range", 0)
ConfigManager.contacts.set_default("rpa_contact_distance", 0.9)
ConfigManager.contacts.set_default("rpa_contact_undecidable_range", 0)
ConfigManager.contacts.set_default("rcb_pairing_contact_distance", 7)
ConfigManager.contacts.set_default("rcb_pairing_contact_undecidable_range", 0)
ConfigManager.contacts.set_default("rcb_stacking_contact_distance", 2.4)
ConfigManager.contacts.set_default("rcb_stacking_contact_undecidable_range", 0)
ConfigManager.contacts.set_default("rc_pairing_contact_distance", 1.5)
ConfigManager.contacts.set_default("rc_pairing_contact_undecidable_range", 0)
ConfigManager.contacts.set_default("rc_stacking_contact_distance", 4.9)
ConfigManager.contacts.set_default("rc_stacking_contact_undecidable_range", 0)
ConfigManager.contacts.set_default("ni_contact_distance", 7.4)
ConfigManager.contacts.set_default("ni_contact_undecidable_range", 0.0)
ConfigManager.contacts.set_default("cacbx_contact_distance", 0.75)
ConfigManager.contacts.set_default("cacbx_undecidable_range", 0.05)
# pylint: enable=no-member


def CaCbxContact():     # pylint: disable=invalid-name
    # class-like name is required here
    """Function producing ca-cbx contact criterion.

    Returned criterion is an alternative of basic CaContact and conjunction of basic CbxCriterion and CaCbxSubtractionCriterion.
    """
    return ContactsAlternative(CaContact(),
                               ContactsConjunction(CbxContact(),
                                                   CaCbxSubtractionCriterion()))


def for_monomer_type_only(type_1, type_2=None):
    """Class decorator used to assert correct type of monomers for subsequent contact evaluation.

    This decorator sets type_1 and type_2 class attributes to provided values. It is implemented mostly for backward
    compatibility.

    Arguments:
    type_1 -- class of monomer for first mer.
    type_2 -- class of monomer for second mer; initially set to None, if so type_1 is taken as type_2.
    """

    def proper_type_decorator(criterion_class):
        """Decorator returning a criterion class with type_1 and type_2 class attributes set to provided values."""

        criterion_class.set_types_cls(type_1, type_2)

        return criterion_class

    return proper_type_decorator


def skip_missing(criterion_class):
    """Class decorator that allows to ignore errors raised by criterions that needs attributes not present in all mers.

    Argument:
    criterion_class -- subclass of ContactCriterion to be wrapped.

    Decorator wraps _is_in_contact method that raises CannotCalculateContact error and returns 0 insted.
    """

    original_is_in_contact = criterion_class._is_in_contact

    def wrapped_is_in_contact(self, *args, **kwargs):
        """Wrapped is_in_contact method that returns 0 instead of raising CannotCalculateContact error."""
        try:
            return original_is_in_contact(self, *args, **kwargs)
        except:
            return 0

    criterion_class._is_in_contact = wrapped_is_in_contact

    return criterion_class


def Not(criterion_class):
    """Class decorator that changes _is_in_contact method.

    Argument:
    criterion_class -- instance of ContactCriterion class.

    Wrapped _is_in_contact method returns opposite values than oryginal method: 0 for 2, 1 for 1 and 2 for 0.
    """

    original_is_in_contact = criterion_class._is_in_contact

    def wrapped_is_in_contact(self, *args, **kwargs):
        """Wrapped _is_in_contact method that returns 0 instead of raising CannotCalculateContact error."""
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
    def new_mth(self, monomer_1, monomer_2, *args, **kwargs):
        if self.type_2 is None or self.type_1 == self.type_2:
            if self.type_1 is None:
                return ory_mth(self, monomer_1, monomer_2, *args, **kwargs)
            else:
                if self._test_type(monomer_1, self.type_1) and self._test_type(monomer_2, self.type_1):
                    return ory_mth(self, monomer_1, monomer_2, *args, **kwargs)
        else:
            if self._test_type(monomer_1, self.type_1) and self._test_type(monomer_2, self.type_2):
                return ory_mth(self, monomer_1, monomer_2, *args, **kwargs)
            elif self._test_type(monomer_1, self.type_2) and self._test_type(monomer_2, self.type_1):
                return ory_mth(self, monomer_2, monomer_1, *args, **kwargs)

        msg_tup = (monomer_1, monomer_2, self.type_1, self.type_2, self.__class__)
        raise WrongMonomerType(*msg_tup)

    new_mth.__doc__ = ory_mth.__doc__ + "\n\nThis method checks monomer types and calls proper method which does actual job."

    return new_mth


class ContactCriterion(object):

    """Abstract class, criteria instances."""

    __metaclass__ = ABCMeta

    type_1 = None
    type_2 = None

    _test_type_cache = {}

    @staticmethod
    def _test_type(monomer, mtype):
        """
        Checks if monomer matches a given type.

        This is a wrapper for the isinstance built-in which caches most frequent queries.

        Arguments:
            monomer -- monomer
            mtype -- type
        """

        if monomer.__class__ == mtype:
            return True

        try:
            (good, bad) = ContactCriterion._test_type_cache[mtype]

            if monomer.__class__ in good:
                return True
            elif monomer.__class__ in bad:
                return False
        except KeyError:
            pass

        if isinstance(monomer, mtype):
            try:
                ContactCriterion._test_type_cache[mtype][0].append(monomer.__class__)
            except KeyError:
                ContactCriterion._test_type_cache[mtype] = ([monomer.__class__], [])
            return True

        try:
            ContactCriterion._test_type_cache[mtype][1].append(monomer.__class__)
        except KeyError:
            ContactCriterion._test_type_cache[mtype] = ([], [monomer.__class__])

        return False


    @check_type
    def is_in_contact_nopre(self, monomer_1, monomer_2, *args, **kwargs):
        """Like is_in_contact, but withour precheck."""
        return self._is_in_contact(monomer_1, monomer_2, *args, **kwargs)

    def is_in_contact(self, monomer_1, monomer_2, *args, **kwargs):
        """Returns three-valued logic contact value.

        This method checks monomer types and calls _is_in_contact_nopre which does actual job.

        Arguments:
        monomer_1 -- first monomer instance.
        monomer_2 -- second monomers instance.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more information.
        """

        if not self._precheck(monomer_1, monomer_2, **kwargs):
            return 0

        return self.is_in_contact_nopre(monomer_1, monomer_2, *args, **kwargs)

    @check_type
    def calculate_distance(self, monomer_1, monomer_2, *args, **kwargs):
        """Calculates distance evaluated by current criterion.

        Arguments:
        monomer_1, monomer_2 -- pydesc.monomer.Monomer subclass instances, for which distance is to be calculated.
        """
        return self._calculate_distance(monomer_1, monomer_2, *args, **kwargs)

    def set_types(self, type_1, type_2=None):
        if type_1 == type_2:
            type_2 = None

        if type_1 == pydesc.monomer.Monomer and type_2 == None:
            type_1 = None

        self.type_1 = type_1
        self.type_2 = type_2


    def get_types(self):
        type_1 = pydesc.monomer.Monomer if self.type_1 is None else self.type_1
        type_2 = type_1 if self.type_2 is None else self.type_2

        return (type_1, type_2)

    set_types_cls = classmethod(set_types)

    @abstractmethod
    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, lazy=True):
        """Abstract method overridden in subclasses.

        Returns three-valued logic contact value.

        This method is called by is_in_contact, which is supposed to check monomer types.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more information.
        """
        pass

    def _precheck(self, monomer_1, monomer_2, rcdist=None, **kwargs):
        """A method for checking a quick and easy precondition of a contact.

        This method should be overriden whenever possible to provide a quick and dirty checking,
        before computing actual criteria.

        This method accepts arguments of any type and returns a boolean.

        Ideally a condition tested here should be simpler even than type checking.
        """

        try:
            if rcdist is not None:
                return rcdist <= self.max_rc_dist
            else:
                return abs(monomer_1.rc - monomer_2.rc) <= self.max_rc_dist
        except: # If anything goes wrong, just disregard the whole precheck.
            pass

        return True

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

    @property
    def criteria(self):
        """Returns list containing current object."""
        return [self]

    def __or__(self, othercc):
        """Returns ContactsAlternative of self and other contact criterion"""
        return ContactsAlternative(self, othercc)

    def __and__(self, othercc):
        """Returns ContactsConjunction of self and other contact criterion"""
        return ContactsConjunction(self, othercc)

    def __xor__(self, othercc):
        """Returns ContactsExclusiveDisjuntion of self and other contact criterion"""
        return ContactsExclusiveDisjunction(self, othercc)


class PointsDistanceCriterion(ContactCriterion):

    """Abstract class, criteria instances."""

    __metaclass__ = ABCMeta
    monomer_hallmark = None

    def __init__(self, criterion_distance=None, undecidable_range=None):
        """Contact criterion constructor.

        Arguments:
        criterion_distance -- radius of the sphere of a center at a given point, inside of which all points are in contact with the central point.
        undecidable_range -- distance from the surface of the sphere, at which all points assume contact value 1 according to the  three-valued logic. Innitially set to 0.

        See also config file docstring.
        """
        self._criterion_distance = criterion_distance
        self._undecidable_range = undecidable_range
        self.max_rc_dist = self.criterion_distance + self.undecidable_range + 10

    @property
    def criterion_distance(self):
        if self._criterion_distance:
            return self._criterion_distance
        return getattr(ConfigManager.contacts, self.monomer_hallmark + '_contact_distance')

    @property
    def undecidable_range(self):
        if self._undecidable_range:
            return self._undecidable_range
        return getattr(ConfigManager.contacts, self.monomer_hallmark + '_contact_undecidable_range')

    def __repr__(self):
        return '<Contact criterion based on %s distance>' % (self.monomer_hallmark,)

    def __str__(self):
        return '%s distance criterion' % (self.monomer_hallmark,)

    def _calculate_distance(self, monomer1obj, monomer2obj):
        """Calculates the distance between two given monomers.

        Arguments:
        monomer1obj -- first monomer instacne.
        monomer2obj -- second monomers instacne.

        Returns the distance between monomers' points in a given unit.
        Returns None if given Monomers do not have appropriate attribute.
        """
        try:
            return (getattr(monomer1obj, self.monomer_hallmark) - (getattr(monomer2obj, self.monomer_hallmark))).calculate_length()
        except AttributeError:
            raise CannotCalculateContact(monomer1obj, monomer2obj, self)

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        """
        distance = self._calculate_distance(monomer_1_obj, monomer_2_obj)
        min_value = self.criterion_distance - self.undecidable_range
        max_value = self.criterion_distance + self.undecidable_range
        if distance is not None:
            if distance <= min_value:
                return 2
            elif distance >= max_value:
                return 0
            else:
                return 1
        else:
            return 0


@for_monomer_type_only(pydesc.monomer.Residue)
class CaContact(PointsDistanceCriterion):

    """Carbon alfa distance criterion."""

    monomer_hallmark = "ca"


@for_monomer_type_only(pydesc.monomer.Residue)
class CbxContact(PointsDistanceCriterion):

    """C-beta extended points (carbon beta extended by 1 Angstrom) distance criterion."""

    monomer_hallmark = "cbx"


class RcContact(PointsDistanceCriterion):

    """Geometrical center distance criterion."""

    monomer_hallmark = "rc"


@for_monomer_type_only(pydesc.monomer.Nucleotide)
class RingCenterContact(PointsDistanceCriterion):

    """Nucleotide ring center distance criterion."""

    monomer_hallmark = "ring_center"

    def __init__(self, *args, **kwargs):
        """RingCenterContact constructor, extended PointsDistanceCriterion method."""
        PointsDistanceCriterion.__init__(self, *args, **kwargs)  # pylint: disable=no-member
        self.max_rc_dist = self.max_rc_dist - 2  # for the sake of compatibilty!


class VectorDistanceCriterion(ContactCriterion):    # pylint: disable=abstract-class-little-used
    # this class will be used more in future

    """Superclass for criteria based on operations made on mers represented by vectors."""

    __metaclass__ = ABCMeta

    def __init__(self, monomer_hallmark_1, monomer_hallmark_2):
        """Vector distance criterion constructor.

        Arguments:
        monomer_hallmark_1, monomer_hallmark_2 -- strings, names of monomer atoms, pseudoatoms or attributes that are
        pydesc.geometry.Coord instances and are to be, respectivelu, beggining and ending of vectore that represents
        monomer in current criterion.
        """
        self.monomer_hallmarks = (monomer_hallmark_1, monomer_hallmark_2)

    def __repr__(self):
        strings = tuple(self.monomer_hallmarks)
        return '<Contact criterion based on operation on vectors of %s and %s>' % strings

    def __str__(self):
        strings = tuple(self.monomer_hallmarks)
        return 'operation on vectors of %s and %s' % strings

    @abstractmethod
    def _is_in_contact(self, monomer_1, monomer_2, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more information.
        """
        pass


@for_monomer_type_only(pydesc.monomer.Residue)
class CaCbxSubtractionCriterion(VectorDistanceCriterion):

    """Criterion based on difference between carbon alpha and extended carbon beta distances.

    Checks if difference of distances between given residues cas' and cbxs' fits given range.
    Range is defined in configuration manager. See:
    pydesc.config.ConfigManager.contacts.cacbx_undecidable_range
    pydesc.config.ConfigManager.contacts.cacbx_contact_distance
    """

    def __init__(self, cirterion_distance=None, udecidable_range=None):
        """Ca-cbx subtraction criterion constructor.

        Arguments:
        criterion distance -- acceptable difference of distances between alpha carbons and cbx for two mers to be considered contacted.
        udecidable_range -- half of the length of undecidable interval spaning above and below criterion distance.

        Both values are initially set to None. If so - values from configuration manager are taken.
        """
        VectorDistanceCriterion.__init__(self, 'ca', 'cbx')
        if cirterion_distance is None:
            self.criterion_distance = ConfigManager.contacts.cacbx_contact_distance     # pylint: disable=no-member
        else:
            self.criterion_distance = cirterion_distance
        if udecidable_range is None:
            self.undecidable_range = ConfigManager.contacts.cacbx_undecidable_range     # pylint: disable=no-member
        else:
            self.undecidable_range = udecidable_range

    def _is_in_contact(self, monomer_1, monomer_2, **kwargs):
        """Returns three-valued logic contact value depending on given mers positions.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more information.
        """
        distances = [(getattr(monomer_1, hallmark) - getattr(monomer_2, hallmark)).calculate_length() for hallmark in self.monomer_hallmarks]
        difference = distances[0] - distances[1]
        upper_treshold = self.criterion_distance + \
            self.undecidable_range   # pylint: disable=no-member
        lower_treshold = self.criterion_distance - \
            self.undecidable_range   # pylint: disable=no-member
        # undecidable_range and criterion_distance are to be overridden in
        # subclasses (abstract attrs).
        if difference >= upper_treshold:
            return 2
        elif difference >= lower_treshold:
            return 1
        else:
            return 0


class SetDistanceCriterion(ContactCriterion):

    """Abstract class, criteria instances."""

    __metaclass__ = ABCMeta
    monomer_hallmark = None
    monomer_hallmark2 = None

    def __init__(self, criterion_distance, undecidable_range=0, num_of_checked_pairs=2):
        """
        Contact criterion constructor.

        Arguments:
        criterion_distance -- radius of the sphere of a center at a given point, inside of which all points are in contact with the central point.
        undecidable_range -- distance from the surface of the sphere, at which all points assume contact value 1 according to the  three-valued logic. Innitially set to 0.
        num_of_checked_pairs -- number of pairs of atoms for which the average distance should be lower than the given cutoff

        See also config file docstring.
        """
        self.criterion_distance = criterion_distance
        self.undecidable_range = undecidable_range
        self.num_of_checked_pairs = num_of_checked_pairs

    def __repr__(self):
        return '<Contact criterion based on %s distances>' % " and ".join(set(map(str, [self.monomer_hallmark, self.monomer_hallmark2])))

    def __str__(self):
        return '%s distances' % " and ".join(set(map(str, [self.monomer_hallmark, self.monomer_hallmark2])))

    # @profile
    def _calculate_distance(self, monomer_1_obj, monomer_2_obj):
        """Calculates the distance between two given monomers.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.

        Returns the distance between monomers' points in a given unit.
        Takes all atoms into consideration if monomer_hallmarks are none
        """


        def get_hallmark(mark_strnum, mer_obj):
            """Returns hallmarks atoms.

            Arguments:
            mark_strnum -- string, hallmark no.
            mer_obj -- pydesc.monomer.Monomer instance.
            """
            mark = getattr(self, 'monomer_hallmark' + mark_strnum)
            if mark is not None:
                atoms = (getattr(mer_obj, mark))
                if type(atoms) == dict:
                    atoms = atoms.values()
                else:
                    atoms = [atoms]
            else:
                atoms = [a for a in mer_obj]
            return atoms

        vecs1 = [a.vector for a in get_hallmark('', monomer_1_obj)]
        vecs2 = [a.vector for a in get_hallmark('2', monomer_2_obj)]

        dist_mat = scipy.spatial.distance.cdist(vecs1, vecs2)

        return dist_mat


    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more information.
        """
        if not self._precheck(monomer_1_obj, monomer_2_obj):
            return False

        dist_mat = self._calculate_distance(monomer_1_obj, monomer_2_obj)
        min_value = self.criterion_distance - self.undecidable_range
        max_value = self.criterion_distance + self.undecidable_range


        class Vertex(object):

            """Class of vertices for bipartite graphs"""

            def __init__(self):
                """Creates attributes edges and features for new vertex."""
                self.edges = []

            def add_edge(self, v2):
                """Adds given vertes to list of neighbours."""
                self.edges.append(v2)

        class BipartiteGraph(object):

            """Representation of bipartite graphs (not directed)."""

            def __init__(self, vs, us):
                """Creates attributes vertices1, vertices2 and edges.

                Arguments:
                vs -- first layer of vertices.
                us -- second layer of vertices.

                Edges are created automatically basing on list of neighbours in every vertex from 1st layer.
                """
                self.vertices1 = vs
                self.vertices2 = us
                self.edges = [frozenset([v, u]) for v in vs for u in v.edges]

            @staticmethod
            def make_from_adjacency_matrix(mtx):
                """Returns BiparireGraph object based on given adjacency matrix."""
                vs, us = map(range, mtx.shape)
                v_s = [Vertex() for i in vs]
                u_s = [Vertex() for i in us]
                c_dict = dict(zip(us, u_s))
                for r, v in zip(vs, v_s):
                    for i in numpy.where(mtx[r]==-1)[0]:
                        u = c_dict[i]
                        v.add_edge(u)
                        u.add_edge(v)
                return BipartiteGraph(v_s, u_s)

            def BFS(self, verts, match):
                """Performs BFS for Hopcroft-Karp algorithm.

                Arguments:
                verts -- list of free vertices for which augmentation paths are to be found.
                match -- list of already matched edges.

                Method attaches to each vertex its distance (dist attribute) to nearest free vertex from given list.
                Skips edges that are already matched for odd distances. Skips not matched edges for even distances.
                """
                queue = [i for i in verts]
                for i in self.vertices1 + self.vertices2:
                    i.visit = 0
                    i.dist = numpy.inf
                    i.parent = None
                for vert in verts:
                    vert.visit = 1
                    vert.dist = 0
                    vert.parent = None
                run = 1
                while queue:
                    odd = bool(run % 2)
                    thev = queue.pop(0)
                    for theu in thev.edges:
                        if (frozenset([theu, thev]) in match) == odd:
                            continue
                        if theu.visit != 0:
                            continue
                        theu.visit = 1
                        theu.dist = thev.dist + 1
                        theu.parent = theu
                        queue.append(theu)
                    thev.visit = 2
                    run += 1

            def DFS(self, vert):
                """Returns a path from given vertes to nearest free vertes marked by BFS.

                Argument:
                vert -- Vertex instance.

                Returns set of tuples representing edges (containing two vertices).
                Sets dist for each used vertes to infinity, so it could not be used by next DFS.
                """
                path = [vert]
                curr = vert
                while curr.dist != 0:
                    prev = min([n for n in curr.edges if n.dist + 1 == vert.dist])
                    path.append(prev)
                    curr.dist = numpy.inf
                    curr = prev
                curr.dist = numpy.inf
                return set([i for i in zip(path[1:], path[:-1])])

        for res, val in ((2, min_value), (1, max_value)):
            bool_mtx = numpy.sign(dist_mat - val)  # -1 indicates contact below threshold
            free_1, free_2 = map(list, map(set, numpy.where(bool_mtx == -1.)))
            bool_mtx = bool_mtx[free_1,].T[free_2,].T     # removing rows and columns that has no -1
            # Hopcroft-Karp algorithm
            graph = BipartiteGraph.make_from_adjacency_matrix(bool_mtx)
            mtchs = set([])
            free_1 = graph.vertices1
            free_2 = graph.vertices2
            mtchs_len = -1
            while True:
                graph.BFS(free_1, mtchs)
                # if there is no more augmention paths algorithm stops
                # that happens when: there is no more free vertices in any layer
                # OR
                # when length of match set has not been increased during last iteration
                if [] in (free_1, free_2) or len(mtchs) == mtchs_len:
                    break
                mtchs_len = len(mtchs)
                min_d = min([i.dist for i in free_2])
                for theu in [i for i in free_2 if i.dist == min_d]:   # one iteraion goes only for shortest paths
                    try:
                        path = graph.DFS(theu)
                        mtchs = mtchs.symmetric_difference(path)
                        if len(mtchs) >= self.num_of_checked_pairs:
                            # algorithm stops when number of matches is great enough, no need to finish it
                            return res
                    except ValueError:  # no paths for current vertex
                        continue
                usedv, usedu = zip(*mtchs)
                free_1 = [i for i in graph.vertices1 if i not in usedv]
                free_2 = [i for i in graph.vertices1 if i not in usedu]
            # algorithm stops here
        return 0

@for_monomer_type_only(pydesc.monomer.Nucleotide)
class RaContact(SetDistanceCriterion):

    """Nucleotide ring atoms distance criterion."""

    monomer_hallmark = "ring_atoms"
    monomer_hallmark2 = "ring_atoms"

    def __init__(self, range=None, und_range=None, no_pairs=None):
        """RingAtomsAverageDistContact constructor, extended SetDistanceCriterion method."""
        # pylint: disable=no-member
        rng = ConfigManager.contacts.at_contact_distance if range is None else range
        urng = ConfigManager.contacts.at_contact_undecidable_range if range is None else und_range
        nprs = 2 if range is None else no_pairs
        # pylint: enable=no-member
        SetDistanceCriterion.__init__(self, rng, urng, nprs)  # pylint: disable=no-member
        self.max_rc_dist = rng + urng + 18

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RaContact _is_in_contact, extended SetDistanceCriterion method."""
        return SetDistanceCriterion._is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs)


@for_monomer_type_only(pydesc.monomer.Nucleotide)
class AtContact(SetDistanceCriterion):

    """Nucleotide atoms distance criterion."""

    def __init__(self, range=None, und_range=None, no_pairs=None):
        """CaContact constructor, extended PointsDistanceCriterion method."""
        # pylint: disable=no-member
        rng = ConfigManager.contacts.at_contact_distance if range is None else range
        urng = ConfigManager.contacts.at_contact_undecidable_range if range is None else und_range
        nprs = 2 if range is None else no_pairs
        # pylint: enable=no-member
        SetDistanceCriterion.__init__(self, rng, urng, nprs)
        self.max_rc_dist =  rng + urng + 18

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """AtContact _is_in_contact, extended SetDistanceCriterion method."""
        return SetDistanceCriterion._is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs)


@for_monomer_type_only(pydesc.monomer.Nucleotide, pydesc.monomer.Ion)
class NIContact(SetDistanceCriterion):

    """Nucleotide-ion distance criterion."""

    def __init__(self):
        """."""
        SetDistanceCriterion.__init__(self, ConfigManager.contacts.ni_contact_distance, ConfigManager.contacts.ni_contact_undecidable_range, 1)   # pylint: disable=no-member
        self.monomer_hallmark = None
        self.monomer_hallmark2 = 'rc'
        self.max_rc_dist = ConfigManager.contacts.ni_contact_distance + ConfigManager.contacts.ni_contact_undecidable_range + 9

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """NIContact _is_in_contact, extended SetDistanceCriterion method."""

        contact_value = SetDistanceCriterion._is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs)
        return contact_value


class DihedralAngleCriterion(ContactCriterion):

    """Abstract class, criteria instances."""

    __metaclass__ = ABCMeta

    def __init__(self, criterion_distance, undecidable_range=0):
        """Contact criterion constructor.

        Arguments:
        criterion_distance -- radius of the sphere of a center at a given point, inside of which all points are in contact with the central point.
        undecidable_range -- distance from the surface of the sphere, at which all points assume contact value 1 according to the  three-valued logic. Innitially set to 0.

        See also config file docstring.
        """
        self.criterion_distance = criterion_distance
        self.undecidable_range = undecidable_range

    def __repr__(self):
        return '<Contact criterion based on dihedral angle between %ss>' % (self.monomer_hallmark,)     # pylint: disable=no-member
        # monomer_hallmark is an abstract attr

    def __str__(self):
        return 'dihedral angle between %ss' % (self.monomer_hallmark,)     # pylint: disable=no-member
        # monomer_hallmark is an abstract attr

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- ignored. See CombinedCriteria._is_in_contact to get more information.
        """
        # pylint: disable=no-member
        angle = abs((getattr(monomer_1_obj, self.monomer_hallmark)).dihedral_angle_cos(getattr(monomer_2_obj, self.monomer_hallmark)))
        # pylint: enable=no-member
        # monomer_hallmark is an abstract attr
        min_value = self.criterion_distance - self.undecidable_range
        max_value = self.criterion_distance + self.undecidable_range
        if angle <= min_value:
            return 0
        elif angle >= max_value:
            return 2
        else:
            return 1


@for_monomer_type_only(pydesc.monomer.Nucleotide)
class RpaContact(DihedralAngleCriterion):

    """Angle between nucleotide ring planes criterion."""

    def __init__(self):
        """RpaContact constructor, extended DihedralAngleCriterion method."""
        self.monomer_hallmark = "ring_plane"
        DihedralAngleCriterion.__init__(
            self, ConfigManager.contacts.rpa_contact_distance, ConfigManager.contacts.rpa_contact_undecidable_range)    # pylint: disable=no-member

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RpaContact _is_in_contact, extended DihedralAngleCriterion method."""
        return DihedralAngleCriterion._is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs)


class HorizontalBisectorDistanceCriterion(ContactCriterion):

    """Abstract class, criteria instances.

    Methods:
    _is_in_contact - returns contact value of a given contact distance.
    calculate_distance - calculates distance between two given monomers.
    """

    __metaclass__ = ABCMeta

    def __init__(self, criterion_distance, undecidable_range=0):
        """Contact criterion constructor.

        Arguments:
        criterion_distance -- radius of the sphere of a center at a given point, inside of which all points are in contact with the central point.
        undecidable_range -- distance from the surface of the sphere, at which all points assume contact value 1 according to the  three-valued logic. Innitially set to 0.

        See also config file docstring.
        """
        self.criterion_distance = criterion_distance
        self.undecidable_range = undecidable_range

    def __repr__(self):
        return '<Contact criterion based on bisector distance between %ss>' % (self.monomer_hallmark_plane,)        # pylint: disable=no-member

    def __str__(self):
        return ' bisector distance between %ss' % (self.monomer_hallmark_plane,)    # pylint: disable=no-member
        # monomer_hallmark_x is an abstract attr

    def _calculate_distance(self, monomer_1_obj, monomer_2_obj):
        """Calculates the distance between two given monomers.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        """
        # pylint: disable=no-member
        bisector = (getattr(monomer_1_obj, self.monomer_hallmark_plane)).bisection_plane(getattr(monomer_2_obj, self.monomer_hallmark_plane))
        ort_projection_point1 = bisector.ort_projection(getattr(monomer_1_obj, self.monomer_hallmark_point))
        ort_projection_point2 = bisector.ort_projection(getattr(monomer_2_obj, self.monomer_hallmark_point))
        # pylint: enable=no-member
        # monomer_hallmark_x is an abstract attr
        return (ort_projection_point1 - ort_projection_point2).calculate_length()

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- ignored. See CombinedCriteria._is_in_contact to get more information.
        """
        distance = self._calculate_distance(monomer_1_obj, monomer_2_obj)
        min_value = self.criterion_distance - self.undecidable_range
        max_value = self.criterion_distance + self.undecidable_range
        if distance <= min_value:
            return 2
        elif distance >= max_value:
            return 0
        else:
            return 1


@for_monomer_type_only(pydesc.monomer.Nucleotide)
class RcbpDistance(HorizontalBisectorDistanceCriterion):

    """Horizontal distance between ring centers contact criterion (for pairing)."""

    def __init__(self):
        """RcbpDistance constructor, extended HorizontalBisectorDistanceCriterion method for pairing criterion."""
        self.monomer_hallmark_point = "ring_center"
        self.monomer_hallmark_plane = "ring_plane"
        HorizontalBisectorDistanceCriterion.__init__(
            self, ConfigManager.contacts.rcb_pairing_contact_distance, ConfigManager.contacts.rcb_pairing_contact_undecidable_range)    # pylint: disable=no-member
        self.max_rc_dist = ConfigManager.contacts.rcb_pairing_contact_distance + ConfigManager.contacts.rcb_pairing_contact_undecidable_range + 6

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RcbpDistance _is_in_contact, extended HorizontalBisectorDistanceCriterion method."""
        return HorizontalBisectorDistanceCriterion._is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs)


@for_monomer_type_only(pydesc.monomer.Nucleotide)
class RcbsDistance(HorizontalBisectorDistanceCriterion):

    """Horizontal distance between ring centers contact criterion (for stacking)."""

    def __init__(self):
        """RcbsDistance constructor, extended HorizontalBisectorDistanceCriterion method for stacking criterion."""
        self.monomer_hallmark_point = "ring_center"
        self.monomer_hallmark_plane = "ring_plane"
        HorizontalBisectorDistanceCriterion.__init__(
            self, ConfigManager.contacts.rcb_stacking_contact_distance, ConfigManager.contacts.rcb_stacking_contact_undecidable_range)  # pylint: disable=no-member
        self.max_rc_dist = ConfigManager.contacts.rcb_stacking_contact_distance + ConfigManager.contacts.rcb_stacking_contact_undecidable_range + 6

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RcbsDistance _is_in_contact, extended HorizontalBisectorDistanceCriterion method."""
        return HorizontalBisectorDistanceCriterion._is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs)


class VerticalBisectorDistanceCriterion(ContactCriterion):

    """Abstract class, criteria instances.

    Methods:
    _is_in_contact - returns contact value of a given contact distance.
    calculate_distance - calculates distance between two given monomers.
    """

    __metaclass__ = ABCMeta

    def __init__(self, criterion_distance, undecidable_range=0):
        """Contact criterion constructor.

        Arguments:
        criterion_distance -- radius of the sphere of a center at a given point, inside of which all points are in contact with the central point.
        undecidable_range -- distance from the surface of the sphere, at which all points assume contact value 1 according to the  three-valued logic. Innitially set to 0.

        See also config file docstring.
        """
        self.criterion_distance = criterion_distance
        self.undecidable_range = undecidable_range

    def __repr__(self):
        return '<Vertical bisector of %ss distance criterion>' % (self.monomer_hallmark_plane,)     # pylint: disable=no-member

    def __str__(self):
        return 'bisector distance of %ss' % (self.monomer_hallmark_plane,)  # pylint: disable=no-member
        # monomer_hallmark_x is abstract attr

    def _calculate_distance(self, monomer_1_obj, monomer_2_obj):
        """Calculates the distance between two given monomers.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        """
        # pylint: disable=no-member
        bisector = (getattr(monomer_1_obj, self.monomer_hallmark_plane)).bisection_plane(getattr(monomer_2_obj, self.monomer_hallmark_plane))
        ort_projection_point1 = bisector.ort_projection(getattr(monomer_1_obj, self.monomer_hallmark_point))
        ort_projection_point2 = bisector.ort_projection(getattr(monomer_2_obj, self.monomer_hallmark_point))
        return (ort_projection_point1 - getattr(monomer_1_obj, self.monomer_hallmark_point)).calculate_length() + (ort_projection_point2 - getattr(monomer_2_obj, self.monomer_hallmark_point)).calculate_length()
        # pylint: enable=no-member
        # monomer_hallmark_x is abstract attr

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- ignored. See CombinedCriteria._is_in_contact to get more information.
        """
        distance = self._calculate_distance(monomer_1_obj, monomer_2_obj)
        min_value = self.criterion_distance - self.undecidable_range
        max_value = self.criterion_distance + self.undecidable_range
        if distance <= min_value:
            return 2
        elif distance >= max_value:
            return 0
        else:
            return 1


@for_monomer_type_only(pydesc.monomer.Nucleotide)
class RcpDistance(VerticalBisectorDistanceCriterion):

    """Vertical distance between ring centers contact criterion (for pairing)."""

    def __init__(self):
        """RcpDistance constructor, extended VerticalBisectorDistanceCriterion method for pairing criterion."""
        self.monomer_hallmark_point = "ring_center"
        self.monomer_hallmark_plane = "ring_plane"
        VerticalBisectorDistanceCriterion.__init__(
            self, ConfigManager.contacts.rc_pairing_contact_distance, ConfigManager.contacts.rc_pairing_contact_undecidable_range)  # pylint: disable=no-member
        self.max_rc_dist = ConfigManager.contacts.rc_pairing_contact_distance + ConfigManager.contacts.rc_pairing_contact_undecidable_range + 16

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RcpDistance _is_in_contact, extended VerticalBisectorDistanceCriterion method."""
        return VerticalBisectorDistanceCriterion._is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs)


@for_monomer_type_only(pydesc.monomer.Nucleotide)
class RcsDistance(VerticalBisectorDistanceCriterion):

    """Vertical distance between ring centers contact criterion (for stacking)."""

    def __init__(self):
        """RcsDistance constructor, extended VerticalBisectorDistanceCriterion method for stacking criterion."""
        self.monomer_hallmark_point = "ring_center"
        self.monomer_hallmark_plane = "ring_plane"
        VerticalBisectorDistanceCriterion.__init__(
            self, ConfigManager.contacts.rc_stacking_contact_distance, ConfigManager.contacts.rc_stacking_contact_undecidable_range)  # pylint: disable=no-member
        self.max_rc_dist = ConfigManager.contacts.rc_stacking_contact_distance + ConfigManager.contacts.rc_stacking_contact_undecidable_range + 16

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RcsDistance _is_in_contact, extended VerticalBisectorDistanceCriterion method."""
        return VerticalBisectorDistanceCriterion._is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs)


class CombinedContact(ContactCriterion):

    """Abstract class, criteria obtained via logical operations on contact criteria."""

    __metalclass__ = ABCMeta

    def __init__(self, *criteria_objs):
        """Combined criteria constructor.

        Arguments:
        criteria_objs -- basic criteria objects.
        """
        self._criteria = criteria_objs
        if len(self.criteria) < 2:
            raise AttributeError("Not enough criteria given to create combined contact criterion")

    def _repr_operation(self):
        """Returns regular expression operation to be used in __repr__ and __str__."""
        pattern = re.compile('[A-Z]{1}[a-z]*')
        operation = " ".join(
            re.findall(pattern, self.__class__.__name__.replace("Contacts", ""))).lower()
        return operation

    def __repr__(self):
        return '<%s of criteria based on %s>' % (self._repr_operation().capitalize(), " and ".join(map(str, self.criteria)))

    def __str__(self):
        return '%s of %s criteria' % (self._repr_operation(), " and ".join(map(str, self.criteria)))

    @abstractmethod
    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, lazy=True, **kwargs):
        """Returns value of combined contact criterion.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- True or False, initially set to True. Deremines if subcriteria values are to be calculated lazy or not.
        Lazy calculation means that if program is able to assume that criterion is not satisfied during calculation subcriterion - further subcriteria are not calculated.
        """
        pass

    @property
    def criteria(self):
        """Returns sequence of combined criterion criteria."""
        return self._criteria


class ContactsConjunction(CombinedContact):

    """Conjunction of critera given in a list.

    Computes criteria type as an intersection of types accepted by subcriteria.
    Algorithm used to resolve types may fail if multiple ineritance is used.
    Given criteria could be a CombinedContact instance.
    """

    def __init__(self, *criteria_objs):
        """Conjunction of criteria constructor.

        Arguments:
        criteria_objs -- basic criteria objects.
        """
        CombinedContact.__init__(self, *criteria_objs)

        type_1 = pydesc.monomer.Monomer
        type_2 = pydesc.monomer.Monomer


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




    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, lazy=True, **kwargs):
        """Returns contact value under conjunction of the given criteria.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- True or False, initially set to True. Deremines if subcriteria values are to be calculated lazy or not.
        Lazy calculation means that if program is able to assume that criterion is not satisfied during calculation subcriterion - further subcriteria are not calculated.
        """
        values = []
        for contact_criterion in self.criteria:
            try:
                value = contact_criterion._is_in_contact(
                    monomer_1_obj, monomer_2_obj, lazy=lazy, **kwargs)
            except WrongMonomerType:
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
    Abstract class grouping combined contacts which require only some criteria to be satisfied.

    Computes criteria type to meet any type of subcriterias.
    """

    def __init__(self, *criteria_objs):
        """Conjunction of criteria constructor.

        Arguments:
        criteria_objs -- basic criteria objects.
        """
        CombinedContact.__init__(self, *criteria_objs)

        def lcs(types):
            """ Lowest common superclass"""
            if len(types) == 0:
                return None

            mros = [x.mro() for x in types]
            for x in mros[0]:
                if all(x in mro for mro in mros):
                    return x

        types = [crit.get_types() for crit in criteria_objs]

        type_1 = lcs([t[0] for t in types])
        type_2 = lcs([t[1] for t in types])

        self.set_types(type_1, type_2)

        self.max_rc_dist=0
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

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, lazy=True, **kwargs):
        """Returns contact value under an alternative of given criteria.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- True or False, initially set to True. Deremines if subcriteria values are to be calculated lazy or not.
        Lazy calculation means that if program is able to assume that criterion is not satisfied during calculation subcriterion - further subcriteria are not calculated.
        """
        values = []
        for contact_criterion in self.criteria:
            try:
                value = contact_criterion.is_in_contact(monomer_1_obj, monomer_2_obj, lazy=lazy, **kwargs)
            except WrongMonomerType:
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

    def get_validating_subcriterion(self, mer_1, mer_2):
        """Returns subcriterion for which given mers are in contact.

        Arguments:
        mer_1, mer_2 -- pydesc.monomer.Monomer instances.

        Raises ValueError in given mers are not in contact.
        """
        for contact_criterion in self.criteria:
            try:
                return contact_criterion.get_validating_subcriterion(mer_1, mer_2)      # pylint: disable=no-member
                # some of criterions has this attribute
            except AttributeError:
                try:
                    if contact_criterion.is_in_contact(mer_1, mer_2):
                        return contact_criterion
                except WrongMonomerType:
                    continue
            except ValueError:
                continue
        raise ValueError('Given mers are not in contact.')


class ContactsExclusiveDisjunction(ContactsDisjunction):

    """Exclusive Disjunction of criteria given in the list.

    Given criteria could be a CombinedContact instance.
    """

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, lazy=True, **kwargs):
        """Returns contact value under an exclusive disjunction of given criteria.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- True or False, initially set to True. Deremines if subcriteria values are to be calculated lazy or not.
        Lazy calculation means that if program is able to assume that criterion is not satisfied during calculation subcriterion - further subcriteria are not calculated.
        """
        values = []
        for contact_criterion in self.criteria:
            try:
                value = contact_criterion.is_in_contact(monomer_1_obj, monomer_2_obj, lazy=lazy, **kwargs)
            except WrongMonomerType:
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
        """Returns value of contact if given mers have their Contact instance in criterion desc attr.

        Returns contact value if both given mers are central monomers for elements in
        any pydesc.structure.Contact stored in current criterion desc.contacts. Otherwise
        returns zero.

        Arguments:
        monomer_1_obj -- first monomer instacne.
        monomer_2_obj -- second monomers instacne.
        lazy -- always set to None.
        """
        for contact in self.desc.contacts:
            if not max(contact.elements).central_monomer in [monomer_1_obj, monomer_2_obj]:
                continue

            if not min(contact.elements).central_monomer in [monomer_1_obj, monomer_2_obj]:
                continue

            return contact.value

        return 0


# pylint: disable=missing-docstring, protected-access, attribute-defined-outside-init
