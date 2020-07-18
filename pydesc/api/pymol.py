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
    """Format name in manner safe for PyMOL."""
    return string.replace(" ", "_").replace("<", "").replace(">", "")


def draw_structures(structures, split_states=False):
    """ Draw given pydesc structure and add it to registry.

    Args:
        structures: list of structure instances.
        split_states(bool): determines if structures should be loaded as separate
        objects in PyMOL (True) or as states of the same object (default; False).

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
    """Draw single contact as unlabeled distance between mers.

    Args:
        structure: (sub)structure instance. In case of structures -- they have to be
        already drawn in PyMOL (registered). Same requirement applies to parent
        structures of substructure.
        ind1(int): 1st mer ind.
        ind2(int): 2nd mer ind.
        point(string): name of atom or pseudoatom present in both mers at ends of
        contact line to be drawn. Optional.
        contact_name(string): PyMOL name of distance object to be created. Note that
        if this object already exists, newly created distance will be merged with
        existing one.
        gap(float): length of gaps between lines.

    """
    name, state_n = Registry.get_structure_or_parent_id(structure)
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


def draw_contact_maps(contact_maps, structures=None, split_contacts=False, point="rc"):
    """Draw given contact maps.

    Draws two contact maps: one for contacts with contact value 2 only (sure contacts),
    and second one with with contact of value 1 (uncertain contacts). Former is
    colored red and has less gaps between dashes, while latter is colored orange and
    has greater gaps.

    Args:
        contact_maps: list of pydesc contact maps.
        structures: (optional) list of structures corresponding to contact maps.
        split_contacts(bool): determines if each contact is to be drawn as separate
        object. False by default.
        point(str): name of atom or pseudoatom to start and end contact lines at (
        "rc" by default).

    """

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

    for i, contact_map in enumerate(contact_maps):
        if structures is None:
            structure = contact_map.structure
        else:
            structure = structures[i]
        name, _ = Registry.get_structure_or_parent_id(structure)
        contacts_dct = {}
        for ind_pair, value in contact_map:
            contacts_dct.setdefault(value, []).append(ind_pair)

        sure_map_name = None if split_contacts else f"{name}_cm"
        uncertain_map_name = None if split_contacts else f"{name}_pcm"
        mk_map(contacts_dct.get(1, []), uncertain_map_name, gap=0.75)
        mk_map(contacts_dct.get(2, []), sure_map_name, gap=0.25)

        cmd.color("red", sure_map_name)
        cmd.color("orange", uncertain_map_name)


def draw_pseudoatoms(structure, pseudoatom_name, anchor_name=None, split_objects=False):
    """Draw pseudoatoms of choice.

    Mers not having given pseudoatom are skipped.

    Args:
        structure: (sub)structure for which pseudoatoms should be drawn.
        pseudoatom_name(str): name of pseudoatom to be drawn.
        anchor_name(str): optional; name of (pseudoatom) to be an anchor. If anchor
        name is given there will be drawn additional dash starting at anchor and
        ending at pseudoatom.
        split_objects(bool): False by default. If so, all created pseudoatoms are
        single pymol object, otherwise each gets different name.

    """
    name, state = Registry.get_structure_or_parent_id(structure)
    obj_name = f"{name}_{pseudoatom_name}"
    for mer in structure:
        pdb_id = structure.derived_from.converter.get_pdb_id(mer.ind)
        pdb_str = pdb_id.format()
        chain = pdb_id.chain
        if split_objects:
            obj_name = f"{name}_{pdb_id}_{pseudoatom_name}"
        try:
            vector = tuple(getattr(mer, pseudoatom_name).vector)
        except AttributeError:
            continue
        pseudoatom_dct = {"state": state, "resi": pdb_str, "chain": chain}
        psd_name = f"PSD{state}"
        cmd.pseudoatom(obj_name, pos=vector, name=psd_name, **pseudoatom_dct)
        if anchor_name is not None:
            try:
                vector = tuple(getattr(mer, anchor_name).vector)
            except AttributeError:
                continue
            anc_name = f"ANC{state}"
            cmd.pseudoatom(obj_name, pos=vector, name=anc_name, **pseudoatom_dct)
            pseudoatom_selection = (
                f"{obj_name} and resi {pdb_str} and name {psd_name} "
                f"and chain {chain}"
            )
            anchor_selection = (
                f"{obj_name} and resi {pdb_str} and name {anc_name} "
                f"and chain {chain}"
            )
            cmd.bond(pseudoatom_selection, anchor_selection)


def select(substructure, selection_name="sele"):
    """Create PyMOL selection based on given substructure.

    Note that selections do not have states, so selecting the same set of mers from
    different states of NMR or trajectories will create selection for all states.

    Args:
        substructure: substructure derived from drawn (registered) structure.
        selection_name(str): name of PyMOL selection to be created (optional; "sele"
        by default).

    """
    name, state_n = Registry.get_id(substructure.derived_from)
    pdb_ids = list(Everything().specify(substructure))

    def _fmt_id(pdb_id):
        icode = "" if pdb_id.icode is None else pdb_id.icode
        chain = "" if pdb_id.chain is None else f"chain {pdb_id.chain} and"
        return f"({chain} resi {pdb_id.ind}{icode})"

    for mer_id in pdb_ids:
        selection_string = f"{name} and ({_fmt_id(mer_id)})"
        cmd.select(selection_name, selection_string, merge=1)


class Registry:
    """Static class keeping track of drawn structures."""

    desc2mol = {}

    @classmethod
    def register(cls, pydesc_obj, pymol_name, pymol_state):
        """Add given pydesc object to registry with given pymol name and state
        number. That pair forms PyMOL id."""
        pymol_id = (pymol_name, pymol_state)
        if pymol_id in cls.desc2mol.values():
            raise KeyError("Name and stated already in registered.")
        cls.desc2mol[pydesc_obj] = pymol_id

    @classmethod
    def remove(cls, pydesc_obj):
        """Remove given object from registry."""
        del cls.desc2mol[pydesc_obj]

    @classmethod
    def get_id(cls, pydesc_obj):
        """Return PyMOL id of given object or raise KeyError."""
        try:
            return cls.desc2mol[pydesc_obj]
        except KeyError:
            raise KeyError(
                f"Structure {pydesc_obj} is not registered (was it never "
                f"drawn or drawn and removed?)."
            )

    @classmethod
    def get_structure_or_parent_id(cls, pydesc_obj):
        """Return id of given object or structure its derived from (or raise
        KeyError)."""
        try:
            return cls.get_id(pydesc_obj)
        except KeyError:
            return cls.get_id(pydesc_obj.derived_from)
