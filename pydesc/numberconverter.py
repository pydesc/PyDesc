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
Class for convertion of pdb mers numberd to PyDesc inds.

created: 13.03.2014, Tymoteusz 'hert' Oleniecki
"""

import re

from pydesc.warnexcept import warn

from pydesc.mers import ConfigManager


def convert_to_id(pdb_residue):
    """Returns PDB id tuple.

    Argument:
    pdb_residue -- BioPython residue.
    """
    warn("Function convert_to_id is deprecated. Use PDB_id.create_from_pdb_residue instead.", DeprecationWarning, 1)

    return PDBid.create_from_pdb_residue(pdb_residue)


def perform_smith_waterman(models_ids):  # pylint: disable=invalid-name
    """Aligns given PDB_ids using Smith-Waterman algorithm.

    Argument:
    models_ids -- list of lists containing PDB_ids to be aligned.

    Used by NumberConverter to align PDB_ids provided by different models (e.g. NMR) of the same structure.
    """
    aligned_ids = models_ids[0]
    for ids in models_ids[1:]:
        # iteration over subsequent models pdb ids
        sw_matrix = build_smith_waterman_matrix(ids, aligned_ids)
        aligned_ids = go_backwards(sw_matrix, iter(reversed(tuple(ids))), iter(reversed(tuple(aligned_ids))))
    return [aligned_ids]


def build_smith_waterman_matrix(ids, aligned_ids):  # pylint: disable=invalid-name, too-many-locals
    # first letters of names should be written with upper case
    # three locals are defined in order to maintain code in clarity
    """Builds and returns Smith-Waterman matrix.

    Arguments:
    ids -- sequence of horizontal PDB_ids to be aligned.
    aligned_ids -- sequence of vertical PDB_ids to be aligned.

    Returns SW_matrix - list of lists containing scoring tuples. Scoring tuple is a tuple containing value of alignment
    represented by current matrix item, traceback (coordinates of next matrix item to choose) and item prefered
    to choose if all adjacent items have the same score.
    """
    gap_open = 0
    match = 1
    sw_matrix = [[(0, (0, 0))] + [(0, (0, i)) for i, dummy in enumerate(ids)]]
    # pylint: disable=cell-var-from-loop
    for row, id1 in enumerate(aligned_ids):
        # notice that row is delayed by 1 due to additional row
        sw_matrix.append([(0, (row, 0))])
        for col, id2 in enumerate(ids):
            # the same happened with column, while additional column
            # is constructed
            compare = lambda: sw_matrix[-2][col][0] + match if id2 == id1 else sw_matrix[-2][col][0]
            sw_values = (sw_matrix[-1][col][0] + gap_open, sw_matrix[-2][col + 1][0] + gap_open, compare())
            # Smith-Waterman values for, respectively: vertical
            # move, horizontal move, match (diagonal move)
            traceback = ((row + 1, col), (row, col + 1), (row, col))
            # traceback is a tuple of matrix coordinates of
            # appropriate matrix entries
            id_to_choose = (id2, id1, id1)
            # ids connected choose of, respectively, vertical,
            # horizontal and diagonal move
            scoring_tuples = zip(sw_values, traceback, id_to_choose)
            # scoring_tuple consists of: appropriate Smith-Waterman value, traceback and
            # id to choose if the value is the highest
            key = lambda scoring_tuple: (scoring_tuple[0],) + scoring_tuple[2][:2]
            # key produces tuple of scoring tuple features that are to be submitted to evaluation
            # features to evaluate are, respectively: --Smith-Waterman value,
            # -- chain character of residue chain
            # -- residue pdb number
            # -- insertion code
            sw_matrix[-1].append(min(scoring_tuples, key=key))
    return sw_matrix
    # pylint: enable=cell-var-from-loop
    # creating locals in loop is the point of above loop


def go_backwards(sw_matrix, vertical_iterator, horizontal_iterator):
    """Aligns two sequences contained by Smith-Waterman matrix.

    Arguments:
    sw_matrix -- matrix created by build_smith_waterman_matrix function.
    vertical_iterator, horizontal_iterator -- iterators or generators that iterates over appropriate sequences.

    Returns list of aligned items.
    """
    row = len(sw_matrix) - 1
    col = len(sw_matrix[0]) - 1
    new_aligned_ids = []
    while True:
        break_now = False
        new_row, new_col = sw_matrix[row][col][1]
        addition = None
        if row - 1 == new_row and col - 1 == new_col:
            addition = next(vertical_iterator)
            next(horizontal_iterator)
        elif row - 1 == new_row:
            addition = next(horizontal_iterator)
        elif col - 1 == new_col:
            addition = next(vertical_iterator)
        if new_row == 0 and new_col == 0:
            break_now = True
        row, col = new_row, new_col
        new_aligned_ids.append(addition)
        if break_now:
            break
    return [i for i in reversed(new_aligned_ids)]


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
        pdb_id -- string in format <chian><pdb_number><pdb_insertion_code>, e.g. C12A.
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
    """Class of objects responsible for converting PDB names of mers to PyDesc integers (inds)."""

    def __init__(self, pdb_models):
        """NumberConverter constructor.

        Arguments:
        pdb_models -- list of pdb models from Bio.PDB parser.

        Sets attributes dict_ind_to_pdb and dict_pdb_to_ind which provides translation between PDB_id and
        PyDesc integers (inds).

        For structures containing few models - performs Smith-Waterman algorithm to establish correct order
        of mers common for all models.
        """
        def is_not_water(pdb_residue):
            import pdb; pdb.set_trace()
            return False if pdb_residue.get_resname() in ConfigManager.mers.solvent else True

        models_ids = []
        for pdb_model in pdb_models:
            no_solvent_ids = filter(is_not_water, pdb_model.get_residues())
            models_ids.append([PDBid.create_from_pdb_residue(id_) for id_ in no_solvent_ids])

        if len(models_ids) > 1:
            models_ids = perform_smith_waterman(models_ids)

        self.dict_ind_to_pdb = dict([(pair[0] + 1, pair[1]) for pair in enumerate(models_ids[0])])
        self.dict_pdb_to_ind = dict(map(tuple, map(reversed, self.dict_ind_to_pdb.items())))

    def get_pdb_id(self, ind):
        """Returns list of PDB-id tuples.

        PDB-id tuple consist of:
        - chain character (index 0) - could be any ASCII alphabetic character. NOTE: case sensitive!
        - PDB number (index 1) - int!
        - insertion code (index 2) - additional alphabetic character that distinguishes mers described with the same PDB number.

        Argument:
        ind -- PyDesc integer.
        """
        return self.dict_ind_to_pdb[ind]

    def get_list_of_inds(self, pdb_id_tuples, distinguish_chains=True):
        """Returns list of PyDesc integers (inds) of mers corresponding to given list of PDB_id
        or tuples containing chain character (string), monomer integer (ind) and insertion code (string; ' ' if none).

        None is returned for any key absent in structure.
        Arguments:
        pdb_id_tuples -- list of PDB id instances or proper tuples (see PDB_id object docstring for more information).
        distinguish_chains -- initially set to True; if so, chain character is considered, otherwise integers of monomer of given pdb number are returned, regardles for chain character.
        """
        if not distinguish_chains:
            sliced_ids = [id_[1:] for id_ in pdb_id_tuples]
            for key in self.dict_pdb_to_ind:
                if key[1:] in sliced_ids and key not in pdb_id_tuples:
                    pdb_id_tuples.append(key)
        return [self.dict_pdb_to_ind[pdb_id] if pdb_id in self.dict_pdb_to_ind else None for pdb_id in pdb_id_tuples]

    def get_ind(self, pdb_id):
        """Takes PDB_id or apropriate tuple and returns proper monomer PyDesc integer (ind).

        Argument:
        pdb_id -- PDB id object or tuple containg proper pdb id.
        """
        try:
            return self.dict_pdb_to_ind[pdb_id]
        except KeyError:
            if pdb_id[2] is None:
                pdb_id = pdb_id[:-1] + (' ',)
                msg = 'Values in structure dict contain space instead of None'
            elif pdb_id[2] == ' ':
                pdb_id = pdb_id[:-1] + (None,)
                msg = "Given PDB_id contains ' ' instead of None"
            else:
                raise KeyError(pdb_id)
            warn(
                DeprecationWarning("' ' as empty insertion code is no longer supperted. Use None insetad. (%s)." % msg))

            return self.dict_pdb_to_ind[pdb_id]
