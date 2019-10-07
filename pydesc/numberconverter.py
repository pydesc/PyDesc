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
Class for conversion of pdb mers numbers to PyDesc inds.

created: 13.03.2014, Tymoteusz 'hert' Oleniecki
"""

import re

import numpy as np

from pydesc.mers import ConfigManager


def perform_smith_waterman(models_ids):  # pylint: disable=invalid-name
    """Aligns given PDB_ids using Smith-Waterman algorithm.

    Argument:
    models_ids -- list of lists containing PDB_ids to be aligned.

    Used by NumberConverter to align PDB_ids provided by different models
    (e.g. NMR) of the same structure.
    """
    aligned_ids = models_ids[0]
    for ids in models_ids[1:]:
        sw_matrix = build_smith_waterman_matrix(ids, aligned_ids)
        aligned_ids = go_backwards(sw_matrix)
    return [aligned_ids]


def build_smith_waterman_matrix(ids, aligned_ids):
    """Builds and returns Smith-Waterman matrix.

    Arguments:
    ids -- sequence of horizontal PDB_ids to be aligned.
    aligned_ids -- sequence of vertical PDB_ids to be aligned.

    Returns SW_matrix - list of lists containing scoring tuples. Scoring
    tuple is a tuple containing value of alignment represented by current
    matrix item, traceback (coordinates of next matrix item to choose) and
    item preferred to choose if all adjacent items have the same score.
    """

    def compare(column, row_, id1_, id2_):
        previous_value = sw_matrix[row_ - 1, column - 1, 0]
        return previous_value + (match if id1_ == id2_ else 0)

    def key(insertion_tuple):
        """
        Produces tuple of scoring tuple features that are to be submitted to
        evaluation features to evaluate are, respectively:
        -- Smith-Waterman value,
        -- chain character of residue chain
        -- residue pdb number
        -- insertion code
        """
        return (insertion_tuple[0],) + insertion_tuple[3:]

    def pack_tuple(tup):
        val, (x, y), pdb_id = tup
        chain = pdb_id[0]
        ind = pdb_id[1]
        icode = pdb_id[2]
        return val, x, y, ord(chain), ind, ord(icode or '\x00')

    match = 1
    sw_matrix = np.full((len(aligned_ids) + 1, len(ids) + 1, 6), None)
    sw_matrix[0, :, 0] = 1
    sw_matrix[:, 0, 0] = 1
    # dimes: seq1 + 1, seq2 + 1,
    #  (score, back_trace_x_coord, back_trace_y_coord, chain, id, icode)
    for row, id1 in enumerate(aligned_ids, 1):
        for col, id2 in enumerate(ids, 1):
            # note that row and cols are shifted by 1
            v1 = sw_matrix[row - 1, col, 0]
            v2 = sw_matrix[row, col - 1, 0]
            v3 = compare(col, row, id1, id2)
            sw_values = (v1, v2, v3)
            # Smith-Waterman values for, respectively: vertical
            # move, horizontal move, match (diagonal move)
            traceback = ((row - 1, col), (row, col - 1), (row - 1, col - 1))
            # traceback is a tuple of matrix coordinates of
            # appropriate matrix entries
            id_to_choose = (id2, id1, id1)
            # ids connected choose of, respectively, vertical,
            # horizontal and diagonal move
            zipper = zip(sw_values, traceback, id_to_choose)
            scoring_tuples = [pack_tuple(insertion_tuple) for insertion_tuple
                              in zipper]
            sw_matrix[row, col] = max(scoring_tuples, key=key)
    return sw_matrix


def go_backwards(sw_matrix):
    """Aligns two sequences contained by Smith-Waterman matrix.

    Arguments:
    sw_matrix -- matrix created by build_smith_waterman_matrix function.

    Returns list of aligned items.
    """

    def unpack_tuple(tup):
        return chr(tup[0]), tup[1], None if not tup[2] else chr(tup[2])

    row, col, dummy = sw_matrix.shape
    row, col = row - 1, col - 1
    new_aligned_ids = []
    while True:
        new_row, new_col, *tup = sw_matrix[row, col, 1:]
        new_aligned_ids.append(unpack_tuple(tup))
        if (new_row, new_col) == (0, 0):
            break
        row, col = new_row, new_col
    return [i for i in reversed(new_aligned_ids)]


class UnknownPDBid(Exception):
    pass


class PDBid(tuple):
    """Tuple type subclass, stores monomer PDB id.

    PDB_id is a tuple containing subsequent values:
    -- chain id (string)
    -- pdb integer (integer)
    -- insertion code (string; ' ' if there is no insertion code).
    """

    def __str__(self):
        chain = '?' if self.chain is None else self.chain
        icode = '' if self.icode is None else self.icode
        return (chain + str(self.ind) + icode).strip()

    @property
    def chain(self):
        """Chain ID"""
        return self[0]

    @property
    def ind(self):
        """Index"""
        return self[1]

    @property
    def icode(self):
        """Insertion code"""
        return self[2]

    @staticmethod
    def create_from_string(pdb_id):
        """Returns PDB id tuple.

        Argument:
        pdb_id -- string in format <chain><pdb_number><pdb_insertion_code>,
        e.g. C12A.
        """
        match = re.match('^(.)([0-9]*)([^0-9])?$', pdb_id)

        if match is None:
            raise ValueError("Unexpected id string %s\n" % pdb_id)

        tuple_ = match.groups()
        icode = None if tuple_[2] in (' ', '') else tuple_[2]
        return pdb_id([tuple_[0], int(tuple_[1]), icode])

    @staticmethod
    def create_from_pdb_residue(pdb_residue):
        """Returns PDB id tuple.

        Argument:
        pdb_residue -- BioPython residue.
        """
        id_tuple = pdb_residue.get_full_id()
        icode = None if id_tuple[3][2] == ' ' else id_tuple[3][2]
        return PDBid((id_tuple[2], int(id_tuple[3][1]), icode))


class NumberConverter(object):
    """Class of objects responsible for converting PDB names of mers to PyDesc
    integers (inds).
    """

    def __init__(self, pdb_models):
        """NumberConverter constructor.

        Arguments:
        pdb_models -- list of pdb models from Bio.PDB parser.

        Sets attributes ind2pdb and pdb2ind which provides translation between
        PDB_id and PyDesc integers (inds).

        For structures containing few models - performs Smith-Waterman
        algorithm to establish correct order of mers common for all models.
        """

        def is_not_water(pdb_residue):
            return not pdb_residue.get_resname() in ConfigManager.mers.solvent

        factory_mth = PDBid.create_from_pdb_residue
        models_ids = []
        for pdb_model in pdb_models:
            no_solvent_ids = filter(is_not_water, pdb_model.get_residues())
            models_ids.append([factory_mth(id_) for id_ in no_solvent_ids])

        if len(models_ids) > 1:
            models_ids = perform_smith_waterman(models_ids)

        self.ind2pdb = models_ids[0]
        self.pdb2ind = {tuple(v): i for i, v in enumerate(self.ind2pdb)}

    def get_pdb_id(self, ind):
        """Returns list of PDB-id tuples.

        PDB-id tuple consist of:
        - chain character (index 0) - could be any ASCII alphabetic character.
        NOTE: case sensitive!
        - PDB number (index 1) - int!
        - insertion code (index 2) - additional alphabetic character that
        distinguishes mers described with the same PDB number.

        Argument:
        ind -- PyDesc integer.
        """
        return self.ind2pdb[ind]

    def get_list_of_inds(self, pdb_id_tuples):
        """Returns list of PyDesc integers (inds) of mers corresponding to
        given list of PDB_id or tuples containing chain character (string),
        monomer integer (ind) and insertion code (string; ' ' if none).

        None is returned for any key absent in structure.
        Arguments:
        pdb_id_tuples -- list of PDB id instances or proper tuples (see
        PDB_id object docstring for more information).
        """
        inds = []
        for key in pdb_id_tuples:
            try:
                ind = self.get_ind(key)
            except UnknownPDBid:
                ind = None
            inds.append(ind)
        return inds

    def get_ind(self, pdb_id):
        """Takes PDB_id or appropriate tuple and returns proper monomer
        PyDesc integer (ind).

        Argument:
        pdb_id -- PDB id object or tuple congaing proper pdb id.
        """
        try:
            return self.pdb2ind[pdb_id]
        except KeyError:
            raise UnknownPDBid('Given PDB was not present in structure.')
