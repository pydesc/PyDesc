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
PyMOL plugin with classes and functions to visualize PyDesc descriptors.

created: 19.09.2013 - , Tymoteusz 'hert' Oleniecki
"""

from pymol import cmd
from pydesc.structure.files import PDBWriter


def _fmt(string):
    return string.replace(" ", "_").replace("<", "").replace(">", "")


def draw_structures(structures):
    """Create representation of given structures in PyMOL.

    Argument:
    structures -- list of pydesc structures each of which will be loaded as separate
    state of new object in PyMOL.

    """
    for structure_index, structure_obj in enumerate(structures):
        pdb_stream = PDBWriter.create_pdb_string(structure_obj, True)
        name = structure_obj.name
        state_n = structure_index + 1
        cmd.read_pdbstr(
            pdb_stream, name, state=state_n,
        )
        Registry.register(structure_obj, name, state_n)


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


class Registry:
    desc2mol = {}

    @classmethod
    def register(cls, pydesc_obj, pymol_name, pymol_state):
        cls.desc2mol[pydesc_obj] = (pymol_name, pymol_state)

    @classmethod
    def get_name(cls, pydesc_obj):
        return cls.desc2mol[pydesc_obj]

    @classmethod
    def get_structure_or_parent_name(cls, pydesc_obj):
        try:
            return cls.get_name(pydesc_obj)
        except KeyError:
            return cls.get_name(pydesc_obj.derived_from)


if __name__ == "pydesc.descmol":
    pass
