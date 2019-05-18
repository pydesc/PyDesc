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
# from pymol.cgo import load_cgo
import tkSimpleDialog
import tkMessageBox
import os
from pydesc import structure
from pydesc import mers
from pydesc import selection
from pydesc import geometry
from pydesc.contacts import contactmap


def fmt(string):
    return string.replace(' ', '_').replace("<", "").replace(">", "")

def register_pydesc():
    import sys
    sys.path.append(os.path.dirname(__file__))


class Attachment(object):

    """Container for objects attached to any structure.

    Methods:
    attach -- adds an attached object.
    """

    def __init__(self, structure_obj):
        """Attachment constructor.

        Argument:
        structure_obj -- instance of pydesc.structure.AbstractStructure.
        """
        structure_obj.trt_matrix.attachment = self
        self.main_object = structure_obj
        self.attached_objects = list()

    def __iter__(self):
        """Returns iterator that iterates over Attachement main object and and all attached objects."""
        return iter([self.main_object] + self.attached_objects)

    def attach(self, substructure):
        """Attaches given object to main object.

        Attached objects moves, when main object is moved and forces main object and all attached object to move, when they are moved.
        """
        if substructure not in self:
            self.attached_objects.append(substructure)
        else:
            raise Warning("Object %s is already attached to %s" % (str(substructure), str(self.main_object)))


class Entry(object):

    """Entry to the Registry.

    Stores data that provide connection between PyDesc and PyMOL.
    """

    def __init__(self, pymol_object_name, pymol_state, attachment, selection_flag=False):
        self.name = pymol_object_name
        self.state = pymol_state
        self.attachment = attachment
        self.selection = selection_flag

    def __repr__(self):
        return "<Registry entry connected with PyMOL %s>" % self.name


class Registry(object):

    """Static class, that stores all objects created in PyDesc and visualized in PyMOL.

    Registry is also responsible for generating names for some substructures that are created during work.
    Attributes:
    _name_counter -- int, used to create automattically generated names.
    objects -- the dictionary connecting PyDesc structures as keys and entries as values; entries stores PyMOL object names for instace.
    attachemnts -- the list of attachments instances.
    Methods:
    add -- adds object (structure or selection) to the registry.
    delete -- deletes object from registry.
    get -- returns data stored in entries.
    get_next_number -- returns new number used in automatically generated names.
    """

    _name_counter = 0
    objects = {}
    attachments = []

    @staticmethod
    def add(structure_obj, name, state=1, selection_flag=False, structure_flag=False):
        """Adds entry to the registry.

        Arguments:
        structure_obj -- a structure or selection from PyDesc.
        name -- name of the connected PyMOL object.
        state -- PyMOL object state connected with given PyDesc objec.
        selection_falg -- flag for selections.
        structure_flag -- flag that provides distinguish between structures and substructures.
        """
        attachment = None
        if selection_flag:
            Registry.objects[structure_obj] = Entry(name, state, None, True)
        else:
            if structure_obj.derived_from in Registry.objects and structure_obj.derived_from != structure_obj:
                attachment = Registry.objects[structure_obj.derived_from].attachment
                attachment.attach(structure_obj)
            if structure_obj not in Registry.objects:
                if attachment is None:
                    Registry.objects[structure_obj] = Entry(name, state, Attachment(structure_obj), selection_flag)
                else:
                    Registry.objects[structure_obj] = Entry(name, state, attachment, selection_flag)
            else:
                raise Warning("This object is already registred")

    @staticmethod
    def delete(structure_obj):
        """Removes entry conected with given structure instance form registry.

        Argument:
        structure_obj -- an instace of PyDesc structure.
        """
        try:
            del Registry.objects[structure_obj]
        except KeyError:
            raise Warning("Given structure is not registred")

    @staticmethod
    def get(structure_obj, *args):
        """Returns attributes of entry connected with given structure instance.

        Arguments:
        args -- down to one string - name of attributes that are to be returned.
        """
        if len(args) == 1:
            return getattr(Registry.objects[structure_obj], args[0])
        return [getattr(Registry.objects[structure_obj], argument) for argument in args]

    @staticmethod
    def get_next_number():
        """Returns name_counter value and raises it by 1"""
        Registry._name_counter += 1
        return Registry._name_counter - 1


class SupportingStructure(structure.AbstractStructure):

    """Structures that represents PyMOL objects in PyDesc."""

    def __init__(self, name, derived_from):
        structure.AbstractStructure.__init__(self, derived_from)
        self.name = name
        self._monomers = []


# ======
# PATCHE
# ======


def load_structure_decorator(current_load_structure):
    """Returns wrapped load_structures if given method is not already wrapped.

    Wrapped method returns list of PyDesc structures and creates PyMOL objects.
    Argument:
    current_load_structure -- structure loader's method load_structures from PyDesc structure module.
    """
    if current_load_structure.__name__ == "wrapped_load_structure":
        return current_load_structure

    def wrapped_load_structure(self, code, *args, **kwargs):
        """Wrapped load_structures method: loads structures from databases or from local copy of PDB file.

        Creates PyMOL objects and returns a list of PyDesc structures connected with those objects.
        Arguments:
        code -- str, pdb code; check PyDesc documentation for more information.
        """
        list_of_structures = current_load_structure(self, code, *args, **kwargs)
        for structure_index, structure_obj in enumerate(list_of_structures):
            cmd.read_pdbstr(structure_obj.create_pdb_string(transformed=False).read(), structure_obj.name, state=structure_index + 1)  # pylint: disable=protected-access
            Registry.add(structure_obj, structure_obj.name, structure_index + 1, structure_flag=True)
            cmd.set_title(structure_obj.name, structure_index + 1, "PD")
        return list_of_structures
    return wrapped_load_structure


def color_segments(self, light_color="gray80", dark_color="gray40"):
    """PyMOL method that colors descriptor segments.

    Arguments:
    light_color -- string, PyMOL color name for 5-mer segments.
    dark_color -- string, PyMOL color name for longer segments.
    """
    for segment in self.segments:
        temp_selection = segment.select()
        temp_selection.select(
            self.derived_from, name="temporary__pydesc_selection")
        Registry.delete(temp_selection)
        if len(segment) == 5:
            cmd.color(light_color, "temporary__pydesc_selection")
        else:
            cmd.color(dark_color, "temporary__pydesc_selection")
    temp_selection = self.central_element.select()
    temp_selection.select(
        self.derived_from, name="temporary__pydesc_selection")
    Registry.delete(temp_selection)
    cmd.color("tv_orange", "temporary__pydesc_selection")
    cmd.delete("temporary__pydesc_selection")


def create_model(self, name=None):
    """Creates PyMOL object connected with current structure (adds entry to Registry)"""
    if name is None:
        name = self.derived_from.name + "_" + \
            self.__class__.__name__[
                :4].lower() + str(Registry.get_next_number())
    Registry.add(self, name, 0, structure_flag=True)
    cmd.read_pdbstr(self.create_pdb_string(enumerate_atoms=True, transformed=True).read(), name, state=1)    # pylint: disable=protected-access
    cmd.set_title(name, 0, "PD")


def get_pymol_matrix(self):
    """Returns transformation-rotation-transformation matrix in PyMOL format."""
    matrix = [
        self.rotation_matrix[0][0], self.rotation_matrix[0][
            1], self.rotation_matrix[0][2], self.translation_vector[0],
        self.rotation_matrix[1][0], self.rotation_matrix[1][
            1], self.rotation_matrix[1][2], self.translation_vector[1],
        self.rotation_matrix[2][0], self.rotation_matrix[2][
            1], self.rotation_matrix[2][2], self.translation_vector[2]
    ] + list(self.prerotational_translation_vector) + [1.0]
    return matrix


def iter_pseudoatoms(self):
    """Returns iterator that iterates over all pseudoatoms in current object."""
    return iter(self.pseudoatoms.values())


def poke_pymol(self, matrix):
    """Forces PyMOL object to move.

    Argument:
    matrix -- transformation-rotation-transformation matrix in PyMOL format.
    """
    for attached_structure in Registry.objects[self.derived_from].attachment:
        if cmd.get_title(*Registry.get(attached_structure, "name", "state")) == "PD":
            cmd.transform_object(Registry.get(attached_structure, "name"), matrix)


def rotate_decorator(current_rotate_method):
    """Returns wrapped rotate method, if given method is not already wrapped.

    Wrapped method affect PyDesc trt_matrix and PyMOL related objects matrices.
    Argument:
    current_rotate_method -- structure's rotate method from PyDesc structure module.
    """
    if current_rotate_method.__name__ == "wrapped_rotate":
        return current_rotate_method

    def wrapped_rotate(self, rotation_matrix):
        """Changes trt matrix and affects PyMOL related objects."""
        self.trt_matrix._update_transformation()    # pylint: disable=protected-access
        current_rotate_method(self, rotation_matrix)
        matrix = [
            rotation_matrix[0][0], rotation_matrix[
                0][1], rotation_matrix[0][2], 0.0,
            rotation_matrix[1][0], rotation_matrix[
                1][1], rotation_matrix[1][2], 0.0,
            rotation_matrix[2][0], rotation_matrix[
                2][1], rotation_matrix[2][2], 0.0,
            0.0, 0.0, 0.0, 1.0
        ]
        self.poke_pymol(matrix)
    return wrapped_rotate


def select_in_pymol(self, structure_obj, name="sele", distinguish_chains=None):
    """Creates PyMOL selection related to self."""
    if distinguish_chains is None:
        distinguish_chains = self._distinguish_chains   # pylint: disable=protected-access
    get_pdb_ind = lambda id_: str(id_.ind) + str(id_.icode) if id_.icode is not None else str(id_.ind)
    cmd.select(name, "%s and (%s)" % (Registry.get(structure_obj, 'name'), " + ".join(["(chain " + id_.chain + " and resi " + get_pdb_ind(id_) + ")" for id_ in self.specify(structure_obj, distinguish_chains).ids])))
    cmd.enable(name)
    Registry.add(self, name, structure_obj, selection_flag=True)


def show(self, point_name, color_segments=True):    # pylint: disable=redefined-outer-name
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
            obj_name = 'PyDesc_obj'
        central_name = "%s_%s_%i_cent" % (obj_name, point_name, self.central_element.central_monomer.ind)
        centers_name = "%s_%s_%i" % (obj_name, point_name, self.central_element.central_monomer.ind)
        distances_name = "%s_desc_%i_%s" % (obj_name, self.central_element.central_monomer.ind, point_name)
        cmd.pseudoatom(central_name, pos=list(getattr(self.central_element.central_monomer, point_name).get_coord(self.trt_matrix)))
        for i, element in enumerate(contact.elements):
            if element != self.central_element:
                try:
                    cmd.pseudoatom(centers_name, pos=getattr(element.central_monomer, point_name).get_coord(self.trt_matrix))
                except AttributeError:
                    # for mers without proper attribute that are incorporated
                    pass
            cmd.pseudoatom('temp_desc_dist%i' % i, pos=getattr(element.central_monomer, point_name).get_coord(self.trt_matrix))
        cmd.distance(distances_name, 'temp_desc_dist0', 'temp_desc_dist1')
        cmd.delete('temp_desc_dist0')
        cmd.delete('temp_desc_dist1')
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

# wlasne cylindry
# def show(self, color_segments = True):
    # obj = dict((criterion.monomer_hallmark, []) for contact in self.primary_contacts for criterion in contact.criterion.criteria)
    # for contact in self.primary_contacts:
        # for element in contact.elements:
        # if element != self.central_element:
        # for subcriterion in contact.criterion.criteria:
        # obj[subcriterion.monomer_hallmark] += [CYLINDER]
        # obj[subcriterion.monomer_hallmark].extend(getattr(element.central_monomer, subcriterion.monomer_hallmark).get_coord(self.trt_matrix))
        # obj[subcriterion.monomer_hallmark].extend(getattr(self.central_element.central_monomer, subcriterion.monomer_hallmark).get_coord(self.trt_matrix))
        # obj[subcriterion.monomer_hallmark] += [0.04, 1,1,1,1,1,1]
        # #clylinder radius, 2x RGB colors for both ends
        # cmd.pseudoatom("%s_%s_%i" % (self.derived_from.name, subcriterion.monomer_hallmark, self.central_element.central_monomer.ind) , pos = getattr(element.central_monomer, subcriterion.monomer_hallmark).get_coord(self.trt_matrix))
    # for hallmark in obj:
        # cmd.load_cgo(obj[hallmark], "%s_desc_%i_%s" % (self.derived_from.name, self.central_element.central_monomer.ind, hallmark))
    # temp_selection = self.select()
    # temp_selection.select(self.derived_from, name = "temporary__pydesc_selection")
    # cmd.orient("temporary__pydesc_selection")
    # cmd.delete("temporary__pydesc_selection")
    # if color_segments == True:
        # self.color_segments()

# linie zamiast cylindrow
# def show(self):
    # for index, contact in enumerate(self.primary_contacts):
        # obj = [BEGIN, LINES, COLOR, 1, 0, 0]
        # for element in contact.elements:
        # obj.extend([VERTEX] + getattr(element.central_monomer, contact.criterion.monomer_hallmark).get_coord(self.trt_matrix))
        # obj += [END]
        # cmd.load_cgo(obj, "desc%i" % index)


def show_contacts(self, point='rc', split_contacts=False, add_name='_cmap'):
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
                print('Skipping contact between %s and %s due to lack of pseudoatom.' % (self.substructure[ind1].pid, self.substructure[ind1].pid))
            p1n = fmt(str(self.substructure) + "_" + str(ind1) + '_' + point)
            p2n = fmt(str(self.substructure) + "_" + str(ind2) + '_' + point)
            cmd.pseudoatom(p1n, pos=p1.get_coord(trt))
            cmd.pseudoatom(p2n, pos=p2.get_coord(trt))
            name = p1n + '_' + p2n if split_contacts else self.substructure.name + add_name
            cmd.distance(fmt(name), p1n, p2n)
            points.extend([p1n, p2n])
    if not split_contacts:
        for i in points:
            cmd.delete(i)
    cmd.hide('labels')
    # pylint: enable=invalid-name


def transform_decorator(current_transform_method):
    """Returns wrapped transform method, if given method is not already wrapped.

    Argument:
    current_transform_method -- TRTMatrix transform method from PyDesc structure module.
    """
    if current_transform_method.__name__ == "wrapped_transform":
        return current_transform_method

    # pylint: disable=invalid-name
    def wrapped_transform(self, *args, **kwargs):
        self._update_transformation()   # pylint: disable=protected-access
        return current_transform_method(self, *args, **kwargs)
    wrapped_transform.__doc__ = current_transform_method.__doc__
    return wrapped_transform
    # pylint: enable=invalid-name


def translate_decorator(current_translate_method):
    """Returns wrapped translate if given method is not already wrapped.

    Wrapped method translates PyDesc structure and related PyMOL object.
    Argument:
    current_translate_method -- structure's method translate from PyDesc structure module.
    """
    if current_translate_method.__name__ == "wrapped_translate":
        return current_translate_method

    def wrapped_translate(self, vector):
        """Translates PyDesc structure and related PyMOL object.

        Argument:
        vector -- list of three coordinates.
        """
        self.trt_matrix._update_transformation()    # pylint: disable=protected-access
        current_translate_method(self, vector)
        matrix = [
            1.0, 0.0, 0.0, vector[0],
            0.0, 1.0, 0.0, vector[1],
            0.0, 0.0, 1.0, vector[2],
            0.0, 0.0, 0.0, 1.0
        ]
        self.poke_pymol(matrix)
    return wrapped_translate


def _update_transformation(self):
    """Checks if PyMOL object matrix has changed and updates self, if so."""
    try:
        if cmd.get_title(*Registry.get(self.attachment.main_object, "name", "state")) == "PD":
            matrix = cmd.get_object_matrix(Registry.get(self.attachment.main_object, "name"))
            self.reset()
            self.add_rotation([matrix[:3], matrix[4:7], matrix[8:11]])
            self.add_translation([matrix[3], matrix[7], matrix[11]])
            self.add_prerotational_translation(list(matrix[12:15]))
    except AttributeError:
        return


def convert_to_pml_mtx(self):
    """Returns current TRTMatrix as flat python list"""
    mtx = [i for i in self.rotation_matrix.flat]
    vec = list(self.translation_vector)
    return mtx[:3] + vec[0:1] + mtx[3:6] + vec[1:2] + mtx[6:] + vec[2:3] + list(self.prerotational_translation_vector) + [1]

def transform_object(self, trtm):
    cmd.transform_object(Registry.get(self, "name"), trtm.convert_to_flat_pymol_matrix())


def w_load(app):
    """Creates structure loader and calls load_structures method with parameters given by user in dialog window."""
    loader = structure.StructureLoader()
    pdb_code = tkSimpleDialog.askstring("PyDesc", "Please enter a 4-char pdb code:", parent=app.root)
    try:
        loader.load_structures(pdb_code)
    except:
        tkMessageBox.showerror("Invalid code", "The code you entered is incorrect:" + pdb_code)
        raise


def __init__(self):
    self.menuBar.addcascademenu("Plugin", "PyDesc")
    self.menuBar.addmenuitem("PyDesc", "command", "Load structure",
                             label="Load structure", command=lambda x=self: w_load(x))
# ===
# INIT
# ===
# pylint: disable=protected-access
if __name__ == "pydesc.descmol":
    cmd.hide("lines")
    cmd.show("cartoon")
    cmd.show("sticks")
    cmd.set("cartoon_side_chain_helper")
    # ABSTRACT STRUCTURE
    structure.AbstractStructure.create_model = create_model
    structure.AbstractStructure.poke_pymol = poke_pymol
    structure.AbstractStructure.rotate = rotate_decorator(structure.AbstractStructure.rotate)
    structure.AbstractStructure.translate = translate_decorator(structure.AbstractStructure.translate)
    structure.AbstractStructure.transform_object = transform_object
    # LOAD STRUCTURE
    structure.StructureLoader.load_structures = load_structure_decorator(structure.StructureLoader.load_structures)
    # DESCRIPTOR
    structure.AbstractDescriptor.show = show
    structure.AbstractDescriptor.color_segments = color_segments
    # SELECTION
    selection.Selection.select = select_in_pymol
    # TRTMatrix
    geometry.TRTMatrix._update_transformation = _update_transformation
    geometry.TRTMatrix.get_pymol_matrix = get_pymol_matrix
    geometry.TRTMatrix.transform = transform_decorator(geometry.TRTMatrix.transform)
    geometry.TRTMatrix.convert_to_flat_pymol_matrix = convert_to_pml_mtx
    # MONOMER
    mers.Mer.iter_pseudoatoms = iter_pseudoatoms
    # CONTACTMAP
    contactmap.ContactMap.show_contacts = show_contacts
# pylint: enable=protected-access
