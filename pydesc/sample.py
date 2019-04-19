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
Tymoteusz
nie usuwac, bedzie rozbudowany przy tworzeniu tutoriala
"""

from pydesc.structure import StructureLoader


PROTEIN_PDBS = ['1A24', '2BLL', '2JRM', '2LJP', '3BIP', '3FV6', '3G67', '3J96', '3NPU', '3PSC', '4LTT', '4NJ6', '4ONK', '4YYN', '4ZTD', '5ERB', '5IFH', '5LF9', '5MPV', '5X55']


def get(stc_name="1no5", ind=0):
    """Returns structure of given PDB id.

    Arguments:
    stc_name -- PDB id; str. 1NO5 by default.
    ind -- 0 by default. Number of model to be returned. Pass 'all' to get all models.
    """
    sl = StructureLoader()
    res = sl.load_structures(stc_name)
    if 'all' != ind:
        res = res[ind]
    return res
