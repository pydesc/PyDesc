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
"""Structural alignment tools.

Contains classes required to solve structure alignment.

created: 4.12.2014 - Tymoteusz 'vdhert' Oleniecki
"""

from pydesc.structure import AbstractDescriptor
from pydesc.cydesc.compdesc import compdesc

import numpy
import operator


def make_basic_evaluation(clicque):
    """Basic evaluation method based on length of clicque vertices list.
    """
    return len(clicque.vertices)


class Graph(object):

    """Class that stores graphs as set of vertices and edges.
    """

    def __init__(self, vrt, edg):
        """Graph constructor.

        Arguments:
        vrt - any iterable object containing objects to be vertices of created graph.
        edg - list of tuples containing pairs of objects from vrt argument, between which there is an edge.
        """
        self.vertices = [i for i in vrt]
        self._vhash = {v: i for i, v in enumerate(vrt)}
        self.edges = [i for i in edg]
        edg_dict = dict((v, []) for v in vrt)
        for v1, v2 in edg:
            edg_dict[v1].append(v2)
            edg_dict[v2].append(v1)
        self.edges_dict = dict((v, tuple(set(l))) for v, l in list(edg_dict.items()))

    def __getitem__(self, vertex):
        return self.edges_dict[vertex]

    @property
    def adjacency_matrix(self):
        if not hasattr(self, "_adj_mtx"):
            self._mk_adj_mtx()
        return self._adj_mtx

    def _mk_adj_mtx(self):
        tmp_mtx = numpy.zeros([len(self.vertices)] * 2)
        for v1 in self.vertices:
            for v2 in self[v1]:
                i1 = self._vhash[v1]
                i2 = self._vhash[v2]
                tmp_mtx[i1, i2] = tmp_mtx[i2, i1] = 1
        self._adj_mtx = numpy.matrix(tmp_mtx)
        # ~ self._adj_mtx = numpy.matrix([[1 if v2 in self[v1] else 0 for v1 in self.vertices] for v2 in self.vertices])


class AlignedDescriptorsGraph(Graph):

    """Class representing graphs of aligned pairs of descriptors needed to find optimal alignment of two structures.
    """

    def __init__(self, vrt, edg):
        """Extended constructor of superclass Graph.

        Arguments:
        vrt -- list of Alignments representing aligned descriptors.
        edg -- list of tuples containing vertices connected with edge. Connected alignments are noncontradictory.
        """
        Graph.__init__(self, vrt, edg)
        self.desc_dict = {}
        for alg in self.vertices:
            for s in alg.structures:
                try:
                    self.desc_dict[s].append(alg)
                except KeyError:
                    self.desc_dict[s] = [alg]

    def get_alignments(self, desc1, desc2):
        """Returns all alignments of two given descriptors."""
        try:
            return list(
                set(self.desc_dict[desc1]).intersection(set(self.desc_dict[desc2]))
            )
        except KeyError:
            raise ValueError(
                "Given descriptores are not aligned in any vertex of the graph."
            )

    @classmethod
    def build_algdesc_graph(cls, structure1, structure2):
        """Class method that prepares graph of all possible PairAlignments of descriptors from two given structures.

        Arguments:
        structure1, structures2 -- instances of any pydesc.AbstractStructure subclasses to be aligned using descriptor approach.
        """
        descs1 = AbstractDescriptor.create_descriptors(structure1)
        descs2 = AbstractDescriptor.create_descriptors(structure2)
        return cls.build_consistency_graph(
            list(filter(bool, descs1)), list(filter(bool, descs2))
        )

    @staticmethod
    def build_consistency_graph(descs1, descs2):
        """Static method that builds AlignedDescriptorsGraph from two lists of descriptors.

        Arguments:
        descs1, descs2 -- lists of pydesc.Descriptors from 1st and 2nd structure to be aligned respectively.
        """

        def add_rmsd(res):
            res[1].rmsd = res[0]
            return res[1]

        def cpd_if_good(de1, de2):
            res = compdesc(de1, de2)
            ret = []
            for n, i in enumerate(res):
                ss1, ss2 = [
                    sel.create_structure(stc)
                    for sel, stc in zip(i[1].get_selections(), i[1].structures)
                ]
                if (
                    not min(
                        list(
                            map(
                                operator.methodcaller("adjusted_number"),
                                (de1, de2, ss1, ss2),
                            )
                        )
                    )
                    < 2
                ):
                    i[1].mnf = n
                    ret.append(i)
            return ret

        vrts = [
            add_rmsd(cpdres)
            for d1 in descs1
            for d2 in descs2
            for cpdres in cpd_if_good(d1, d2)
        ]
        edgs = [
            (v1, v2)
            for i, v1 in enumerate(vrts)
            for v2 in vrts[i + 1 :]
            if v1.is_consistent(v2)
        ]
        return AlignedDescriptorsGraph(vrts, edgs)

    @staticmethod
    def build_simple_consistency_graph(descs1, descs2):
        """Static method that builds AlignedDescriptorsGraph from two lists of descriptors.

        Arguments:
        descs1, descs2 -- lists of pydesc.Descriptors from 1st and 2nd structure to be aligned respectively.
        """

        def add_rmsd(res):
            try:
                res[0][1].rmsd = res[0][0]
                return res[0][1]
            except IndexError:
                return

        def cpd_if_good(de1, de2):
            res = compdesc(de1, de2)
            ret = []
            for n, i in enumerate(res):
                ss1, ss2 = [
                    sel.create_structure(stc)
                    for sel, stc in zip(i[1].get_selections(), i[1].structures)
                ]
                if (
                    not min(
                        list(
                            map(
                                operator.methodcaller("adjusted_number"),
                                (de1, de2, ss1, ss2),
                            )
                        )
                    )
                    < 2
                ):
                    i[1].mnf = n
                    ret.append(i)
            return ret

        vrts = list(
            filter(
                bool, [add_rmsd(cpd_if_good(d1, d2)) for d1 in descs1 for d2 in descs2]
            )
        )
        edgs = [
            (v1, v2)
            for i, v1 in enumerate(vrts)
            for v2 in vrts[i + 1 :]
            if v1.is_consistent(v2)
        ]
        return AlignedDescriptorsGraph(vrts, edgs)


class AlignedMersGraph(Graph):
    def __init__(self, vrt, edg, mrs_dict):
        """Extended Graph constructor.

        Arguments:
        vrt -- list of tuples containing two aligned mers from different descriptors and alignment as 3rd element.
        edg -- list of tuples containing two connected vertices from vrt list.
        mrs_dict -- dictionary with mers as keys and lists of vertices connected with all vertices containing key mer as values.
        """
        Graph.__init__(self, vrt, edg)
        self.mer_dict = mrs_dict

    @staticmethod
    def build_from_desc_graph(desc_graph):
        """Static method that builds AlignedMersGraph out of AlignedDescriptorsGraph object.

        Argument:
        desc_graph -- instance of AlignedDescriptorsGraph.
        """

        def add_edge(v1, v2):
            """"""
            edgs.append((v1, v2))
            for i, j in ((v1, v2), (v2, v1)):
                try:
                    edg_dict[i[0]].append(j)
                except KeyError:
                    edg_dict[i[0]] = [j]

        vrts = []
        edgs = []
        edg_dict = {}
        for alg in desc_graph.vertices:
            for i, aa_pair in enumerate(alg):
                vrt = aa_pair + (alg,)
                vrts.append(vrt)
                for aa_pair2 in alg:
                    vrt2 = aa_pair2 + (alg,)
                    add_edge(vrt, vrt2)

        for v1, v2 in desc_graph.edges:
            # v1 and v2 are PairAlignments of descriptors
            # nv1, nv2 - stands for "new vertex"
            for aa_pair1 in v1:
                nv1 = aa_pair1 + (v1,)
                for aa_pair2 in v2:
                    nv2 = aa_pair2 + (v2,)
                    add_edge(nv1, nv2)

        return AlignedMersGraph(vrts, edgs, edg_dict)

    @staticmethod
    def get_AM_from_ADG_vertices(stc1, stc2, adg):
        """Returns list of vertices for adjacency matrix created from AlignedDescriptorsGraph in order used during creation of adjacency matrix.

        Arguments:
        stc1 -- 
        stc2 -- 
        adg -- aligned descriptors graph.
        """
        return AlignedMersGraph._getAMfADGv(stc1, stc2, adg)[0]

    @staticmethod
    def _getAMfADGv(stc1, stc2, adg):
        """Returns list of vertices for adjacency matrix created from AlignedDescriptorsGraph and lenght of blocks of submatrices for vertices related to same pairs of mers.

        Arguments:
        stc1 -- 
        stc2 -- 
        adg -- aligned descriptors graph.
        """
        lens = {}

        def get_blocks_len(m1, m2, i):
            try:
                lens[(m1, m2)] += 1
            except KeyError:
                lens[(m1, m2)] = 1
            return (m1, m2, i)

        return (
            [
                get_blocks_len(m1, m2, i)
                for m1 in stc1
                for m2 in stc2
                for i, d in enumerate(adg.vertices)
                if (m1, m2) in d
            ],
            lens,
        )

    @staticmethod
    def make_adjacency_matrix_from_ADG(stc1, stc2, adg):
        """Returns adjacency matrix for aligned mers based on AlignedDescriptorsGraph and structures.

        Arguments:
        stc1 -- 
        stc2 -- 
        adg -- aligned descriptors graph.
        """
        mprs, lens = AlignedMersGraph._getAMfADGv(stc1, stc2, adg)
        dtm_mtx = numpy.zeros([len(mprs), len(adg.vertices)])
        for i, (m1, m2, j) in enumerate(mprs):
            dtm_mtx[i, j] = 1
        dtm_mtx = numpy.matrix(dtm_mtx)

        ad_mtx = numpy.matrix(adg.adjacency_matrix)
        numpy.fill_diagonal(ad_mtx, 1)
        res = dtm_mtx * ad_mtx * dtm_mtx.T

        start_i = 0
        ord = [i[:2] for i in mprs]
        for pair in sorted(set(ord), key=lambda x: ord.index(x)):
            end_i = start_i + lens[pair]
            res[start_i:end_i, start_i:end_i] = 0
            start_i = end_i

        return res


class ClicqueSolver(object):

    """Static class that provides methods to search for maximum clicques in graphs.
    """

    @staticmethod
    def find_maximum_clq(graph, eval_fx):
        """Method that finds maximum clicques and evaluates them with given function.

        Argument:
        graph -- Graph instance to be searched for maximum clicques.
        eval_fx -- function that evaluates clicque (Graph instance). It must take Graph instance as argument and return int or float.
        """
        # ? zamiast tego wstawic solver
        v = graph.vertices
        e = graph.edges
        g = Graph(v, e)
        clicques = [g]
        # ?
        points = list(map(eval_fx, clicques))
        return sorted(zip(points, clicques), key=lambda x: x[0])
