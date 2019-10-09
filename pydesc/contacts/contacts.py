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

from abc import ABCMeta
from abc import abstractmethod

import numpy
import scipy.spatial

import pydesc.mers
from pydesc.config import ConfigManager
from pydesc.warnexcept import CannotCalculateContact
from .base import check_type
from .base import ContactCriterion
from .base import ContactsAlternative
from .base import ContactsConjunction
from .base import for_monomer_type_only

# pylint: disable=no-member
ConfigManager.new_branch("contacts")
ConfigManager.contacts.set_default("ca_contact_distance", 6.0)
ConfigManager.contacts.set_default("ca_contact_undecidable_range", 0.5)
ConfigManager.contacts.set_default("cbx_contact_distance", 6.5)
ConfigManager.contacts.set_default("cbx_contact_undecidable_range", 0.5)
ConfigManager.contacts.set_default("rc_contact_distance", 7.5)
ConfigManager.contacts.set_default("rc_contact_undecidable_range", 0.5)
ConfigManager.contacts.set_default("ring_center_contact_distance", 6.25)
ConfigManager.contacts.set_default("ring_center_contact_undecidable_range",
                                   0.0)
ConfigManager.contacts.set_default("at_contact_distance", 5.0)
ConfigManager.contacts.set_default("at_contact_undecidable_range", 0.)
ConfigManager.contacts.set_default("rpa_contact_distance", 0.9)
ConfigManager.contacts.set_default("rpa_contact_undecidable_range", 0.)
ConfigManager.contacts.set_default("rcb_pairing_contact_distance", 7.)
ConfigManager.contacts.set_default("rcb_pairing_contact_undecidable_range", 0.)
ConfigManager.contacts.set_default("rcb_stacking_contact_distance", 2.4)
ConfigManager.contacts.set_default("rcb_stacking_contact_undecidable_range",
                                   0.)
ConfigManager.contacts.set_default("rc_pairing_contact_distance", 1.5)
ConfigManager.contacts.set_default("rc_pairing_contact_undecidable_range", 0.)
ConfigManager.contacts.set_default("rc_stacking_contact_distance", 4.9)
ConfigManager.contacts.set_default("rc_stacking_contact_undecidable_range", 0.)
ConfigManager.contacts.set_default("ni_contact_distance", 7.4)
ConfigManager.contacts.set_default("ni_contact_undecidable_range", 0.0)
ConfigManager.contacts.set_default("nx_contact_distance", 7.5)
ConfigManager.contacts.set_default("nx_contact_undecidable_range", 0.0)
ConfigManager.contacts.set_default("ncx_contact_distance", 7.5)
ConfigManager.contacts.set_default("ncx_contact_undecidable_range", 0.0)
ConfigManager.contacts.set_default("xnc_contact_distance", 5.75)
ConfigManager.contacts.set_default("xnc_contact_undecidable_range", 0.0)
ConfigManager.contacts.set_default("prc_contact_distance", 7.5)
ConfigManager.contacts.set_default("prc_contact_undecidable_range", 0.0)
ConfigManager.contacts.set_default("cacbx_contact_distance", 0.75)
ConfigManager.contacts.set_default("cacbx_undecidable_range", 0.05)


# pylint: enable=no-member


def CaCbxContact():  # pylint: disable=invalid-name
    # class-like name is required here
    """Function producing ca-cbx contact criterion.

    Returned criterion is an alternative of basic CaContact and conjunction
    of basic CbxCriterion and CaCbxSubtractionCriterion.
    """
    return ContactsAlternative(
        CaContact(),
        ContactsConjunction(
            CbxContact(),
            CaCbxSubtractionCriterion()))


class PointsDistanceCriterion(ContactCriterion, metaclass=ABCMeta):
    """Abstract class, criteria instances."""

    monomer_hallmark = None

    def __init__(self, distance_threshold=None, undecidable_range=None):
        """Contact criterion constructor.

        Arguments:
        distance_threshold -- radius of the sphere of a center at a given
        point, inside of which all points are in contact with the central
        point.
        undecidable_range -- distance from the surface of the sphere,
        at which all points assume contact value 1 according to the
        three-valued logic. Initially set to 0.

        See also config file docstring.
        """
        self._distance_threshold = distance_threshold
        self._undecidable_range = undecidable_range
        distance_threshold = self.distance_threshold
        undecidable_range = self.undecidable_range
        self.max_rc_dist = distance_threshold + undecidable_range + 10

    @property
    def distance_threshold(self):
        """Property returning criterion distance."""
        if self._distance_threshold is not None:
            return self._distance_threshold
        setting_name = '%s_contact_distance' % self.monomer_hallmark
        return getattr(ConfigManager.contacts, setting_name)

    @property
    def undecidable_range(self):
        """Property returning undecidable range."""
        if self._undecidable_range is not None:
            return self._undecidable_range
        setting_name = '%s_contact_undecidable_range' % self.monomer_hallmark
        return getattr(ConfigManager.contacts, setting_name)

    def calculate_distance(self, monomer_1, monomer_2, *args, **kwargs):
        """Calculates distance evaluated by current criterion.

        Arguments:
        monomer_1, monomer_2 -- pydesc.monomer.Monomer subclass instances,
        for which distance is to be calculated.
        """
        return self._calculate_distance(monomer_1, monomer_2, *args, **kwargs)

    def __repr__(self):
        text = '<Contact criterion based on %s distance>'
        return text % self.monomer_hallmark

    def __str__(self):
        return '%s distance criterion' % (self.monomer_hallmark,)

    @check_type
    def _calculate_distance(self, monomer1obj, monomer2obj, *args, **kwargs):
        """Calculates the distance between two given mers.

        Arguments:
        monomer1obj -- first monomer instance.
        monomer2obj -- second mers instance.

        Returns the distance between mers' points in a given unit.
        Returns None if given Monomers do not have appropriate attribute.
        """
        return (getattr(monomer1obj, self.monomer_hallmark) - (
            getattr(monomer2obj, self.monomer_hallmark))).calculate_length()

    def _is_in_contact(self, monomer1obj, monomer2obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer1obj -- first monomer instance.
        monomer2obj -- second mers instance.
        """
        distance = self._calculate_distance(monomer1obj, monomer2obj)
        min_value = self.distance_threshold - self.undecidable_range
        max_value = self.distance_threshold + self.undecidable_range
        if distance <= min_value:
            return 2
        elif distance >= max_value:
            return 0
        else:
            return 1


@for_monomer_type_only(pydesc.mers.Residue)
class CaContact(PointsDistanceCriterion):
    """Carbon alpha distance criterion."""

    monomer_hallmark = "ca"


@for_monomer_type_only(pydesc.mers.Residue)
class CbxContact(PointsDistanceCriterion):
    """C-beta extended points (carbon beta extended by 1 Angstrom) distance
    criterion."""

    monomer_hallmark = "cbx"


class RcContact(PointsDistanceCriterion):
    """Geometrical center distance criterion."""

    monomer_hallmark = "rc"


@for_monomer_type_only(pydesc.mers.Nucleotide)
class RingCenterContact(PointsDistanceCriterion):
    """Nucleotide ring center distance criterion."""

    monomer_hallmark = "ring_center"

    def __init__(self, *args, **kwargs):
        """RingCenterContact constructor, extended PointsDistanceCriterion
        method."""
        PointsDistanceCriterion.__init__(self, *args, **kwargs)
        self.max_rc_dist = self.max_rc_dist - 2  # for sake of compatibility!


@for_monomer_type_only(pydesc.mers.Nucleotide)
class PrcContact(PointsDistanceCriterion):
    """Nucleotide proximate ring center distance criterion."""

    monomer_hallmark = "prc"

    def __init__(self, *args, **kwargs):
        """RingCenterContact constructor, extended PointsDistanceCriterion
        method."""
        PointsDistanceCriterion.__init__(self, *args, **kwargs)
        self.max_rc_dist = self.max_rc_dist - 2


@for_monomer_type_only(pydesc.mers.Nucleotide)
class NxContact(PointsDistanceCriterion):
    """Nucleotide ring center distance criterion."""

    monomer_hallmark = "nx"

    def __init__(self, *args, **kwargs):
        """RingCenterContact constructor, extended PointsDistanceCriterion
        method."""
        PointsDistanceCriterion.__init__(self, *args, **kwargs)
        self.max_rc_dist = self.max_rc_dist - 2


class DifferentPointsDistanceCriterion(ContactCriterion):
    """Abstract class for criteria based on distances between two different
    points from two different mers.

    Methods in those classes should be symmetrical, i.e. _is_in_contact
    and _calculate_distance should return the same result when called for
    *(mer1, mer2) and *(mer2, mer1).
    """

    def __init__(self, criterion_distance=None, undecidable_range=None,
                 mer1_hallmark=None, mer2_hallmark=None):
        """Contact criterion constructor.

        Arguments:
        distance_threshold -- radius of the sphere of a center at a given
        point, inside of which all points are in contact with the central
        point.
        undecidable_range -- distance from the surface of the sphere,
        at which all points assume contact value 1 according to the
        three-valued logic. Initially set to 0.

        See also config file docstring.
        """
        self._criterion_distance = criterion_distance
        self._undecidable_range = undecidable_range
        self.mer1_hallmark = mer1_hallmark
        self.mer2_hallmark = mer2_hallmark
        distance = self.criterion_distance
        buffer = self.undecidable_range
        self.max_rc_dist = distance + buffer + 10
        self.min_value = self.criterion_distance - self.undecidable_range
        self.max_value = self.criterion_distance + self.undecidable_range

    @property
    def criterion_distance(self):
        if self._criterion_distance is not None:
            return self._criterion_distance
        setting_name = "_".join((
            self.mer1_hallmark,
            self.mer2_hallmark,
            'contact_distance'))
        return getattr(ConfigManager.contacts, setting_name)

    @property
    def undecidable_range(self):
        if self._undecidable_range is not None:
            return self._undecidable_range
        setting_name = "_".join((
            self.mer1_hallmark,
            self.mer2_hallmark,
            'contact_undecidable_range'))
        return getattr(ConfigManager.contacts, setting_name)

    def _calculate_distance(self, mer1, mer2, *args, **kwargs):
        res = set([])
        for m1, m2 in ((mer1, mer2), (mer2, mer1)):
            try:
                res.add((getattr(m1, self.mer1_hallmark) - (
                    getattr(m2, self.mer2_hallmark))).calculate_length())
            except AttributeError:
                pass
        try:
            return min(res)
        except ValueError:
            raise CannotCalculateContact(mer1, mer2, self)

    def _is_in_contact(self, mer1, mer2, *args, **kwargs):
        distance = self.calculate_distance(mer1, mer2, *args, **kwargs)
        if distance <= self.min_value:
            return 2
        elif distance >= self.max_value:
            return 0
        else:
            return 1


class VectorDistanceCriterion(ContactCriterion, metaclass=ABCMeta):
    """Superclass for criteria based on operations made on mers represented
    by vectors."""

    def __init__(self, monomer_hallmark_1, monomer_hallmark_2):
        """Vector distance criterion constructor.

        Arguments:
        monomer_hallmark_1, monomer_hallmark_2 -- strings, names of monomer
        atoms, pseudoatoms or attributes that are
        pydesc.geometry.Coord instances and are to be, respectively,
        beelining and ending of vector that represents
        monomer in current criterion.
        """
        self.monomer_hallmarks = (monomer_hallmark_1, monomer_hallmark_2)

    def __repr__(self):
        strings = tuple(self.monomer_hallmarks)
        msg = '<Contact criterion based on operation on vectors of %s and %s>'
        return msg % strings

    def __str__(self):
        strings = tuple(self.monomer_hallmarks)
        return 'operation on vectors of %s and %s' % strings

    @abstractmethod
    def _is_in_contact(self, monomer_1, monomer_2, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more
        information.
        """
        pass


@for_monomer_type_only(pydesc.mers.Residue)
class CaCbxSubtractionCriterion(VectorDistanceCriterion):
    """Criterion based on difference between carbon alpha and extended
    carbon beta distances.

    Checks if difference of distances between given residues cas' and cbxs'
    fits given range.
    Range is defined in configuration manager. See:
    pydesc.config.ConfigManager.contacts.cacbx_undecidable_range
    pydesc.config.ConfigManager.contacts.cacbx_contact_distance
    """

    def __init__(self, cirterion_distance=None, undecidable_range=None):
        """Ca-cbx subtraction criterion constructor.

        Arguments:
        criterion distance -- acceptable difference of distances between
        alpha carbons and cbx for two mers to be considered contacted.
        undecidable_range -- half of the length of undecidable interval
        spanning above and below criterion distance.

        Both values are initially set to None. If so - values from
        configuration manager are taken.
        """
        VectorDistanceCriterion.__init__(self, 'ca', 'cbx')
        if cirterion_distance is None:
            distance = ConfigManager.contacts.cacbx_contact_distance
            self.criterion_distance = distance
        else:
            self.criterion_distance = cirterion_distance
        if undecidable_range is None:
            distance = ConfigManager.contacts.cacbx_undecidable_range
            self.undecidable_range = distance
        else:
            self.undecidable_range = undecidable_range

    def _is_in_contact(self, monomer_1, monomer_2, **kwargs):
        """Returns three-valued logic contact value depending on given mers
        positions.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more
        information.
        """
        distances = []
        for hallmark in self.monomer_hallmarks:
            point1 = getattr(monomer_1, hallmark)
            point2 = getattr(monomer_2, hallmark)
            distance = (point1 - point2).calculate_length()
            distances.append(distance)
        difference = distances[0] - distances[1]
        upper_threshold = self.criterion_distance + self.undecidable_range
        lower_threshold = self.criterion_distance - self.undecidable_range
        if difference >= upper_threshold:
            return 2
        elif difference >= lower_threshold:
            return 1
        else:
            return 0


class SetDistanceCriterion(ContactCriterion, metaclass=ABCMeta):
    """Abstract class, criteria instances."""
    monomer_hallmark = None
    monomer_hallmark2 = None

    def __init__(self, criterion_distance, undecidable_range=0,
                 num_of_checked_pairs=2):
        """
        Contact criterion constructor.

        Arguments:
        distance_threshold -- radius of the sphere of a center at a given
        point, inside of which all points are in contact with the central
        point.
        undecidable_range -- distance from the surface of the sphere,
        at which all points assume contact value 1 according to the
        three-valued logic. Initially set to 0.
        num_of_checked_pairs -- number of pairs of atoms for which the
        average distance should be lower than the given cutoff.

        See also config file docstring.
        """
        self.criterion_distance = criterion_distance
        self.undecidable_range = undecidable_range
        self.num_of_checked_pairs = num_of_checked_pairs

    def __repr__(self):
        return '<Contact criterion based on %s distances>' % " and ".join(
            set(map(str, [self.monomer_hallmark, self.monomer_hallmark2])))

    def __str__(self):
        return '%s distances' % " and ".join(
            set(map(str, [self.monomer_hallmark, self.monomer_hallmark2])))

    def _calculate_distance(self, monomer_1_obj, monomer_2_obj):
        """Calculates the distance between two given mers.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.

        Returns the distance between mers' points in a given unit.
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
                    atoms = list(atoms.values())
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
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- ignored. See CombinedCriteria.is_in_contact to get more
        information.
        """
        if not self._pre_check(monomer_1_obj, monomer_2_obj):
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

                Edges are created automatically basing on list of neighbours in
                every vertex from 1st layer.
                """
                self.vertices1 = vs
                self.vertices2 = us
                self.edges = [frozenset([v, u]) for v in vs for u in v.edges]

            @staticmethod
            def make_from_adjacency_matrix(mtx):
                """Returns BiparireGraph object based on given adjacency
                matrix."""
                vs, us = list(map(range, mtx.shape))
                v_s = [Vertex() for i in vs]
                u_s = [Vertex() for i in us]
                c_dict = dict(list(zip(us, u_s)))
                for r, v in zip(vs, v_s):
                    for i in numpy.where(mtx[r] == -1)[0]:
                        u = c_dict[i]
                        v.add_edge(u)
                        u.add_edge(v)
                return BipartiteGraph(v_s, u_s)

            def BFS(self, verts, match):
                """Performs BFS for Hopcroft-Karp algorithm.

                Arguments:
                verts -- list of free vertices for which augmentation paths
                are to be found.
                match -- list of already matched edges.

                Method attaches to each vertex its distance (dist attribute) to
                nearest free vertex from given list.
                Skips edges that are already matched for odd distances.
                Skips not matched edges for even distances.
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
                """Returns a path from given vertes to nearest free vertes
                marked by BFS.

                Argument:
                vert -- Vertex instance.

                Returns set of tuples representing edges (containing two
                vertices).
                Sets dist for each used vertes to infinity, so it could not
                be used by next DFS.
                """
                path = [vert]
                curr = vert
                while curr.dist != 0:
                    prev = min(
                        [n for n in curr.edges if n.dist + 1 == vert.dist])
                    path.append(prev)
                    curr.dist = numpy.inf
                    curr = prev
                curr.dist = numpy.inf
                return set([i for i in zip(path[1:], path[:-1])])

        for res, val in ((2, min_value), (1, max_value)):
            bool_mtx = numpy.sign(
                dist_mat - val)  # -1 indicates contact below threshold
            free_1, free_2 = list(
                map(list, list(map(set, numpy.where(bool_mtx == -1.)))))
            bool_mtx = bool_mtx[free_1,].T[
                free_2,].T  # removing rows and columns that has no -1
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
                for theu in [i for i in free_2 if
                             i.dist == min_d]:  # one iteraion goes only for shortest paths
                    try:
                        path = graph.DFS(theu)
                        mtchs = mtchs.symmetric_difference(path)
                        if len(mtchs) >= self.num_of_checked_pairs:
                            # algorithm stops when number of matches is great enough, no need to finish it
                            return res
                    except ValueError:  # no paths for current vertex
                        continue
                usedv, usedu = list(zip(*mtchs))
                free_1 = [i for i in graph.vertices1 if i not in usedv]
                free_2 = [i for i in graph.vertices1 if i not in usedu]
            # algorithm stops here
        return 0


@for_monomer_type_only(pydesc.mers.Nucleotide)
class RaContact(SetDistanceCriterion):
    """Nucleotide ring atoms distance criterion."""

    monomer_hallmark = "ring_atoms"
    monomer_hallmark2 = "ring_atoms"

    def __init__(self, range=None, und_range=None, no_pairs=None):
        """RingAtomsAverageDistContact constructor, extended
        SetDistanceCriterion method."""
        # pylint: disable=no-member
        rng = ConfigManager.contacts.at_contact_distance if range is None else range
        urng = ConfigManager.contacts.at_contact_undecidable_range if range is None else und_range
        nprs = 2 if range is None else no_pairs
        # pylint: enable=no-member
        SetDistanceCriterion.__init__(self, rng, urng,
                                      nprs)  # pylint: disable=no-member
        self.max_rc_dist = rng + urng + 18

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RaContact _is_in_contact, extended SetDistanceCriterion method."""
        return SetDistanceCriterion._is_in_contact(self, monomer_1_obj,
                                                   monomer_2_obj, **kwargs)


@for_monomer_type_only(pydesc.mers.Nucleotide)
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
        self.max_rc_dist = rng + urng + 18

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """AtContact _is_in_contact, extended SetDistanceCriterion method."""
        return SetDistanceCriterion._is_in_contact(self, monomer_1_obj,
                                                   monomer_2_obj, **kwargs)


@for_monomer_type_only(pydesc.mers.Nucleotide, pydesc.mers.Ion)
class NIContact(SetDistanceCriterion):
    """Nucleotide-ion distance criterion."""

    def __init__(self):
        """."""
        SetDistanceCriterion.__init__(self,
                                      ConfigManager.contacts.ni_contact_distance,
                                      ConfigManager.contacts.ni_contact_undecidable_range,
                                      1)  # pylint: disable=no-member
        self.monomer_hallmark = None
        self.monomer_hallmark2 = 'rc'
        self.max_rc_dist = ConfigManager.contacts.ni_contact_distance + ConfigManager.contacts.ni_contact_undecidable_range + 9

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """NIContact _is_in_contact, extended SetDistanceCriterion method."""
        contact_value = SetDistanceCriterion._is_in_contact(self,
                                                            monomer_1_obj,
                                                            monomer_2_obj,
                                                            **kwargs)
        return contact_value

    def _calculate_distance(self, *args, **kwargs):
        return super(NIContact, self)._calculate_distance(*args, **kwargs) - \
               args[1].get_radius()


class DihedralAngleCriterion(ContactCriterion, metaclass=ABCMeta):
    """Abstract class, criteria instances."""

    def __init__(self, criterion_distance, undecidable_range=0):
        """Contact criterion constructor.

        Arguments:
        distance_threshold -- radius of the sphere of a center at a given
        point, inside of which all points are in contact with the central
        point.
        undecidable_range -- distance from the surface of the sphere,
        at which all points assume contact value 1 according to the
        three-valued logic. Innitially set to 0.

        See also config file docstring.
        """
        self.criterion_distance = criterion_distance
        self.undecidable_range = undecidable_range

    def __repr__(self):
        return '<Contact criterion based on dihedral angle between %ss>' % (
            self.monomer_hallmark,)  # pylint: disable=no-member
        # monomer_hallmark is an abstract attr

    def __str__(self):
        return 'dihedral angle between %ss' % (
            self.monomer_hallmark,)  # pylint: disable=no-member
        # monomer_hallmark is an abstract attr

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- ignored. See CombinedCriteria._is_in_contact to get more
        information.
        """
        # pylint: disable=no-member
        angle = abs(
            (getattr(monomer_1_obj, self.monomer_hallmark)).dihedral_angle_cos(
                getattr(monomer_2_obj, self.monomer_hallmark)))
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


@for_monomer_type_only(pydesc.mers.Nucleotide)
class RpaContact(DihedralAngleCriterion):
    """Angle between nucleotide ring planes criterion."""

    def __init__(self):
        """RpaContact constructor, extended DihedralAngleCriterion method."""
        self.monomer_hallmark = "ring_plane"
        DihedralAngleCriterion.__init__(
            self, ConfigManager.contacts.rpa_contact_distance,
            ConfigManager.contacts.rpa_contact_undecidable_range)  # pylint: disable=no-member

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RpaContact _is_in_contact, extended DihedralAngleCriterion
        method."""
        return DihedralAngleCriterion._is_in_contact(self, monomer_1_obj,
                                                     monomer_2_obj, **kwargs)


class HorizontalBisectorDistanceCriterion(ContactCriterion, metaclass=ABCMeta):
    """Abstract class, criteria instances.

    Methods:
    _is_in_contact - returns contact value of a given contact distance.
    calculate_distance - calculates distance between two given mers.
    """

    def __init__(self, criterion_distance, undecidable_range=0):
        """Contact criterion constructor.

        Arguments:
        distance_threshold -- radius of the sphere of a center at a given
        point, inside of which all points are in contact with the central
        point.
        undecidable_range -- distance from the surface of the sphere,
        at which all points assume contact value 1 according to the
        three-valued logic. Innitially set to 0.

        See also config file docstring.
        """
        self.criterion_distance = criterion_distance
        self.undecidable_range = undecidable_range

    def __repr__(self):
        return '<Contact criterion based on bisector distance between %ss>' % (
            self.monomer_hallmark_plane,)  # pylint: disable=no-member

    def __str__(self):
        return ' bisector distance between %ss' % (
            self.monomer_hallmark_plane,)  # pylint: disable=no-member
        # monomer_hallmark_x is an abstract attr

    def _calculate_distance(self, monomer_1_obj, monomer_2_obj):
        """Calculates the distance between two given mers.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        """
        # pylint: disable=no-member
        bisector = (getattr(monomer_1_obj,
                            self.monomer_hallmark_plane)).bisection_plane(
            getattr(monomer_2_obj, self.monomer_hallmark_plane))
        ort_projection_point1 = bisector.ort_projection(
            getattr(monomer_1_obj, self.monomer_hallmark_point))
        ort_projection_point2 = bisector.ort_projection(
            getattr(monomer_2_obj, self.monomer_hallmark_point))
        # pylint: enable=no-member
        # monomer_hallmark_x is an abstract attr
        return (
                ort_projection_point1 - ort_projection_point2).calculate_length()

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- ignored. See CombinedCriteria._is_in_contact to get more
        information.
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


@for_monomer_type_only(pydesc.mers.Nucleotide)
class RcbpDistance(HorizontalBisectorDistanceCriterion):
    """Horizontal distance between ring centers contact criterion (for
    pairing)."""

    def __init__(self):
        """RcbpDistance constructor, extended
        HorizontalBisectorDistanceCriterion method for pairing criterion."""
        self.monomer_hallmark_point = "ring_center"
        self.monomer_hallmark_plane = "ring_plane"
        HorizontalBisectorDistanceCriterion.__init__(
            self, ConfigManager.contacts.rcb_pairing_contact_distance,
            ConfigManager.contacts.rcb_pairing_contact_undecidable_range)  # pylint: disable=no-member
        self.max_rc_dist = ConfigManager.contacts.rcb_pairing_contact_distance + ConfigManager.contacts.rcb_pairing_contact_undecidable_range + 6

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RcbpDistance _is_in_contact, extended
        HorizontalBisectorDistanceCriterion method."""
        return HorizontalBisectorDistanceCriterion._is_in_contact(self,
                                                                  monomer_1_obj,
                                                                  monomer_2_obj,
                                                                  **kwargs)


@for_monomer_type_only(pydesc.mers.Nucleotide)
class RcbsDistance(HorizontalBisectorDistanceCriterion):
    """Horizontal distance between ring centers contact criterion (for
    stacking)."""

    def __init__(self):
        """RcbsDistance constructor, extended
        HorizontalBisectorDistanceCriterion method for stacking criterion."""
        self.monomer_hallmark_point = "ring_center"
        self.monomer_hallmark_plane = "ring_plane"
        HorizontalBisectorDistanceCriterion.__init__(
            self, ConfigManager.contacts.rcb_stacking_contact_distance,
            ConfigManager.contacts.rcb_stacking_contact_undecidable_range)  # pylint: disable=no-member
        self.max_rc_dist = ConfigManager.contacts.rcb_stacking_contact_distance + ConfigManager.contacts.rcb_stacking_contact_undecidable_range + 6

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RcbsDistance _is_in_contact, extended
        HorizontalBisectorDistanceCriterion method."""
        return HorizontalBisectorDistanceCriterion._is_in_contact(self,
                                                                  monomer_1_obj,
                                                                  monomer_2_obj,
                                                                  **kwargs)


class VerticalBisectorDistanceCriterion(ContactCriterion, metaclass=ABCMeta):
    """Abstract class, criteria instances.

    Methods:
    _is_in_contact - returns contact value of a given contact distance.
    calculate_distance - calculates distance between two given mers.
    """

    def __init__(self, criterion_distance, undecidable_range=0):
        """Contact criterion constructor.

        Arguments:
        distance_threshold -- radius of the sphere of a center at a given
        point, inside of which all points are in contact with the central
        point.
        undecidable_range -- distance from the surface of the sphere,
        at which all points assume contact value 1 according to the
        three-valued logic. Innitially set to 0.

        See also config file docstring.
        """
        self.criterion_distance = criterion_distance
        self.undecidable_range = undecidable_range

    def __repr__(self):
        return '<Vertical bisector of %ss distance criterion>' % (
            self.monomer_hallmark_plane,)  # pylint: disable=no-member

    def __str__(self):
        return 'bisector distance of %ss' % (
            self.monomer_hallmark_plane,)  # pylint: disable=no-member
        # monomer_hallmark_x is abstract attr

    def _calculate_distance(self, monomer_1_obj, monomer_2_obj):
        """Calculates the distance between two given mers.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        """
        # pylint: disable=no-member
        bisector = (getattr(monomer_1_obj,
                            self.monomer_hallmark_plane)).bisection_plane(
            getattr(monomer_2_obj, self.monomer_hallmark_plane))
        ort_projection_point1 = bisector.ort_projection(
            getattr(monomer_1_obj, self.monomer_hallmark_point))
        ort_projection_point2 = bisector.ort_projection(
            getattr(monomer_2_obj, self.monomer_hallmark_point))
        return (ort_projection_point1 - getattr(monomer_1_obj,
                                                self.monomer_hallmark_point)).calculate_length() + (
                       ort_projection_point2 - getattr(monomer_2_obj,
                                                       self.monomer_hallmark_point)).calculate_length()
        # pylint: enable=no-member
        # monomer_hallmark_x is abstract attr

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """Returns three-valued logic contact value.

        Arguments:
        monomer_1_obj -- first monomer instance.
        monomer_2_obj -- second mers instance.
        lazy -- ignored. See CombinedCriteria._is_in_contact to get more
        information.
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


@for_monomer_type_only(pydesc.mers.Nucleotide)
class RcpDistance(VerticalBisectorDistanceCriterion):
    """Vertical distance between ring centers contact criterion (for
    pairing)."""

    def __init__(self):
        """RcpDistance constructor, extended
        VerticalBisectorDistanceCriterion method for pairing criterion."""
        self.monomer_hallmark_point = "ring_center"
        self.monomer_hallmark_plane = "ring_plane"
        VerticalBisectorDistanceCriterion.__init__(
            self, ConfigManager.contacts.rc_pairing_contact_distance,
            ConfigManager.contacts.rc_pairing_contact_undecidable_range)  # pylint: disable=no-member
        self.max_rc_dist = ConfigManager.contacts.rc_pairing_contact_distance + ConfigManager.contacts.rc_pairing_contact_undecidable_range + 16

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RcpDistance _is_in_contact, extended
        VerticalBisectorDistanceCriterion method."""
        return VerticalBisectorDistanceCriterion._is_in_contact(self,
                                                                monomer_1_obj,
                                                                monomer_2_obj,
                                                                **kwargs)


@for_monomer_type_only(pydesc.mers.Nucleotide)
class RcsDistance(VerticalBisectorDistanceCriterion):
    """Vertical distance between ring centers contact criterion (for
    stacking)."""

    def __init__(self):
        """RcsDistance constructor, extended
        VerticalBisectorDistanceCriterion method for stacking criterion."""
        self.monomer_hallmark_point = "ring_center"
        self.monomer_hallmark_plane = "ring_plane"
        VerticalBisectorDistanceCriterion.__init__(
            self,
            ConfigManager.contacts.rc_stacking_contact_distance,
            ConfigManager.contacts.rc_stacking_contact_undecidable_range
        )
        self.max_rc_dist = ConfigManager.contacts.rc_stacking_contact_distance + \
                           ConfigManager.contacts.rc_stacking_contact_undecidable_range + \
                           16

    def _is_in_contact(self, monomer_1_obj, monomer_2_obj, **kwargs):
        """RcsDistance _is_in_contact, extended
        VerticalBisectorDistanceCriterion method."""
        return VerticalBisectorDistanceCriterion._is_in_contact(self,
                                                                monomer_1_obj,
                                                                monomer_2_obj,
                                                                **kwargs)
