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
Package providing functions for visualization of some PyDesc features in PyMOL.

Usage:
    - make sure Python used by PyMOL has PyDesc installed (how? see below)
    - simply import pydesc in PyMOL console and create objects or use PyMOL to
    execute scripts using pydesc library
    - import pydesc.api.pymol and use its functions to visualize pydesc objects

How to install PyDesc in PyMOL's Python?
Either use pip in PyMOL console or build PyMOL (
https://github.com/schrodinger/pymol-open-source/blob/master) with PREFIX_PATH
variable set to path storing already install PyDesc.

created: 19.09.2013 - , Tymoteusz 'hert' Oleniecki
"""

from pymol import cmd
from pydesc.structure.files import PDBWriter
from pydesc.selection import Everything


def _fmt(string):
    return string.replace(" ", "_").replace("<", "").replace(">", "")


def draw_structures(structures, split_states=False):
    """Create representation of given structures in PyMOL.

    Argument:
    structures -- list of pydesc structures each of which will be loaded as separate
    state of new object in PyMOL.

    """
    for structure_index, structure_obj in enumerate(structures):
        pdb_stream = PDBWriter.create_pdb_string(structure_obj, True)
        name = structure_obj.name
        state_n = structure_index + 1
        if split_states:
            name = f"{name}_{state_n}"
            state_n = 0
        Registry.register(structure_obj, name, state_n)
        cmd.read_pdbstr(pdb_stream, name, state=state_n)


def draw_contact(structure, ind1, ind2, point="rc", contact_name=None, gap=0.5):
    name, state_n = Registry.get_structure_or_parent_name(structure)
    trt = structure.trt_matrix
    p1 = getattr(structure[ind1], point)
    p2 = getattr(structure[ind2], point)
    p1n = _fmt(f"{name}_{ind1}_{point}")
    p2n = _fmt(f"{name}_{ind2}_{point}")
    cmd.pseudoatom(p1n, pos=p1.get_coord(trt), state=state_n)
    cmd.pseudoatom(p2n, pos=p2.get_coord(trt), state=state_n)
    if contact_name is None:
        contact_name = p1n + "_" + p2n
    cmd.distance(_fmt(contact_name), p1n, p2n, gap=gap, label=0, state=state_n)
    cmd.delete(p1n)
    cmd.delete(p2n)


def draw_contact_maps(contact_maps, split_contacts=False, point="rc"):
    def mk_map(contacts, map_name=None, gap=0.5):
        contacts = set([frozenset(i) for i in contacts])
        for ind1, ind2 in contacts:
            try:
                draw_contact(structure, ind1, ind2, point, map_name, gap=gap)
            except AttributeError:
                msg = "Skipping contact between %s and %s due to lack of pseudoatom."
                inp = ind1, ind2
                print(msg % inp)
                continue

    for contact_map in contact_maps:
        structure = contact_map.structure
        name, _ = Registry.get_structure_or_parent_name(structure)
        contacts_dct = {}
        for ind_pair, value in contact_map:
            contacts_dct.setdefault(value, []).append(ind_pair)

        sure_map_name = None if split_contacts else f"{name}_cm"
        uncertain_map_name = None if split_contacts else f"{name}_pcm"
        mk_map(contacts_dct.get(1, []), uncertain_map_name, gap=0.75)
        mk_map(contacts_dct.get(2, []), sure_map_name, gap=0.25)

        cmd.color("red", sure_map_name)
        cmd.color("orange", uncertain_map_name)


def select(substructure, selection_name="sele"):
    name, state_n = Registry.get_name(substructure.derived_from)
    pdb_ids = list(Everything().specify(substructure))

    def _fmt_id(pdb_id):
        icode = '' if pdb_id.icode is None else pdb_id.icode
        chain = '' if pdb_id.chain is None else f'chain {pdb_id.chain} and'
        return f"({chain} resi {pdb_id.ind}{icode})"

    for mer_id in pdb_ids:
        selection_string = f"{name} and ({_fmt_id(mer_id)})"
        cmd.select(selection_name, selection_string, merge=1)


class Registry:
    desc2mol = {}

    @classmethod
    def register(cls, pydesc_obj, pymol_name, pymol_state):
        pymol_id = (pymol_name, pymol_state)
        if pymol_id in cls.desc2mol.values():
            raise KeyError("Name and stated already in registered.")
        cls.desc2mol[pydesc_obj] = pymol_id

    @classmethod
    def remove(cls, pydesc_obj):
        del cls.desc2mol[pydesc_obj]

    @classmethod
    def get_name(cls, pydesc_obj):
        try:
            return cls.desc2mol[pydesc_obj]
        except KeyError:
            raise KeyError(f"Structure {pydesc_obj} is not registered (was it never "
                           f"drawn or drawn and removed?).")

    @classmethod
    def get_structure_or_parent_name(cls, pydesc_obj):
        try:
            return cls.get_name(pydesc_obj)
        except KeyError:
            return cls.get_name(pydesc_obj.derived_from)


if __name__ == "pydesc.descmol":
    pass
