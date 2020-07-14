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


def draw_structures(structures):
    """Create representation of given structures in PyMOL.

    Argument:
    structures -- list of pydesc structures each of which will be loaded as separate
    state of new object in PyMOL.

    """
    for structure_index, structure_obj in enumerate(structures):
        pdb_stream = PDBWriter.create_pdb_string(structure_obj, True)
        cmd.read_pdbstr(
            pdb_stream, structure_obj.name, state=structure_index + 1,
        )


def color_segments(self, light_color="gray80", dark_color="gray40"):
    """PyMOL method that colors descriptor segments.

    Arguments:
    light_color -- string, PyMOL color name for 5-mer segments.
    dark_color -- string, PyMOL color name for longer segments.
    """
    for segment in self.segments:
        temp_selection = segment.select()
        temp_selection.select(self.derived_from, name="temporary__pydesc_selection")
        Registry.delete(temp_selection)
        if len(segment) == 5:
            cmd.color(light_color, "temporary__pydesc_selection")
        else:
            cmd.color(dark_color, "temporary__pydesc_selection")
    temp_selection = self.central_element.select()
    temp_selection.select(self.derived_from, name="temporary__pydesc_selection")
    Registry.delete(temp_selection)
    cmd.color("tv_orange", "temporary__pydesc_selection")
    cmd.delete("temporary__pydesc_selection")


def create_model(self, name=None):
    """Creates PyMOL object connected with current structure (adds entry to Registry)"""
    if name is None:
        name = (
            self.derived_from.name
            + "_"
            + self.__class__.__name__[:4].lower()
            + str(Registry.get_next_number())
        )
    Registry.add(self, name, 0, structure_flag=True)
    cmd.read_pdbstr(
        self.create_pdb_string(enumerate_atoms=True, transformed=True).read(),
        name,
        state=1,
    )  # pylint: disable=protected-access
    cmd.set_title(name, 0, "PD")


def get_pymol_matrix(self):
    """Returns transformation-rotation-transformation matrix in PyMOL format."""
    matrix = (
        [
            self.rotation_matrix[0][0],
            self.rotation_matrix[0][1],
            self.rotation_matrix[0][2],
            self.translation_vector[0],
            self.rotation_matrix[1][0],
            self.rotation_matrix[1][1],
            self.rotation_matrix[1][2],
            self.translation_vector[1],
            self.rotation_matrix[2][0],
            self.rotation_matrix[2][1],
            self.rotation_matrix[2][2],
            self.translation_vector[2],
        ]
        + list(self.prerotational_translation_vector)
        + [1.0]
    )
    return matrix


def poke_pymol(self, matrix):
    """Forces PyMOL object to move.

    Argument:
    matrix -- transformation-rotation-transformation matrix in PyMOL format.
    """
    for attached_structure in Registry.objects[self.derived_from].attachment:
        if cmd.get_title(*Registry.get(attached_structure, "name", "state")) == "PD":
            cmd.transform_object(Registry.get(attached_structure, "name"), matrix)


def select_in_pymol(self, structure_obj, name="sele", distinguish_chains=None):
    """Creates PyMOL selection related to self."""
    if distinguish_chains is None:
        distinguish_chains = (
            self._distinguish_chains
        )  # pylint: disable=protected-access
    get_pdb_ind = (
        lambda id_: str(id_.ind) + str(id_.icode)
        if id_.icode is not None
        else str(id_.ind)
    )
    cmd.select(
        name,
        "%s and (%s)"
        % (
            Registry.get(structure_obj, "name"),
            " + ".join(
                [
                    "(chain " + id_.chain + " and resi " + get_pdb_ind(id_) + ")"
                    for id_ in self.specify(structure_obj, distinguish_chains).ids
                ]
            ),
        ),
    )
    cmd.enable(name)
    Registry.add(self, name, structure_obj, selection_flag=True)


def show(self, point_name, color_segments=True):  # pylint: disable=redefined-outer-name
    # color_segments is a patch-method, not a function, so it is not covered
    """Creates PyMOL visualisation of descriptor.

    Argiment:
    point_name -- name of the atom or pseudoatom to be joined with dashes in pymol.
    color_segments -- initially set to True, if so - calls descriptors method color_segments.
    """
    # CREATING PYMOL OBJECTS
    for contact in self.contacts:
        try:
            obj_name = self.derived_from.name
        except AttributeError:
            obj_name = "PyDesc_obj"
        central_name = "%s_%s_%i_cent" % (
            obj_name,
            point_name,
            self.central_element.central_monomer.ind,
        )
        centers_name = "%s_%s_%i" % (
            obj_name,
            point_name,
            self.central_element.central_monomer.ind,
        )
        distances_name = "%s_desc_%i_%s" % (
            obj_name,
            self.central_element.central_monomer.ind,
            point_name,
        )
        cmd.pseudoatom(
            central_name,
            pos=list(
                getattr(self.central_element.central_monomer, point_name).get_coord(
                    self.trt_matrix
                )
            ),
        )
        for i, element in enumerate(contact.elements):
            if element != self.central_element:
                try:
                    cmd.pseudoatom(
                        centers_name,
                        pos=getattr(element.central_monomer, point_name).get_coord(
                            self.trt_matrix
                        ),
                    )
                except AttributeError:
                    # for mers without proper attribute that are incorporated
                    pass
            cmd.pseudoatom(
                "temp_desc_dist%i" % i,
                pos=getattr(element.central_monomer, point_name).get_coord(
                    self.trt_matrix
                ),
            )
        cmd.distance(distances_name, "temp_desc_dist0", "temp_desc_dist1")
        cmd.delete("temp_desc_dist0")
        cmd.delete("temp_desc_dist1")
        cmd.hide("labels", distances_name)
        cmd.color("tv_orange", distances_name)
    # CREATING PYDESC OBJECTS
    for item in [central_name, centers_name]:
        item_structure = SupportingStructure(item, self.derived_from)
        Registry.add(item_structure, item)
        cmd.set_title(item, 0, "PD")
    temp_selection = self.select()
    temp_selection.select(self.derived_from, name="temporary__pydesc_selection")
    Registry.delete(temp_selection)
    cmd.orient("temporary__pydesc_selection")
    cmd.delete("temporary__pydesc_selection")
    if color_segments is True:
        self.color_segments()


def show_contacts(self, point="rc", split_contacts=False, add_name="_cmap"):
    """Shows distances between all contacted mers.

    Argument:
    point -- string, name of mers pseudoatom to be merged with distance lines.
    """
    # pylint: disable=invalid-name
    trt = self.substructure.trt_matrix
    points = []
    for ind1 in self.contacts:
        for ind2 in self.contacts[ind1]:
            try:
                p1 = getattr(self.substructure[ind1], point)
                p2 = getattr(self.substructure[ind2], point)
            except AttributeError:
                print(
                    (
                        "Skipping contact between %s and %s due to lack of pseudoatom."
                        % (self.substructure[ind1].pid, self.substructure[ind1].pid)
                    )
                )
            p1n = fmt(str(self.substructure) + "_" + str(ind1) + "_" + point)
            p2n = fmt(str(self.substructure) + "_" + str(ind2) + "_" + point)
            cmd.pseudoatom(p1n, pos=p1.get_coord(trt))
            cmd.pseudoatom(p2n, pos=p2.get_coord(trt))
            name = (
                p1n + "_" + p2n if split_contacts else self.substructure.name + add_name
            )
            cmd.distance(fmt(name), p1n, p2n)
            points.extend([p1n, p2n])
    if not split_contacts:
        for i in points:
            cmd.delete(i)
    cmd.hide("labels")
    # pylint: enable=invalid-name


if __name__ == "pydesc.descmol":
    cmd.hide("lines")
    cmd.show("cartoon")
    cmd.show("sticks")
    cmd.set("cartoon_side_chain_helper")
