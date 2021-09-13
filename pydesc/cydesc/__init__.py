# Copyright 2017 Maciek Dziubinski
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
Interface to descriptor procedures written in C.


created: 25.12.2013 - Maciek Dziubinski
"""

import ctypes
import ctypes.util
import os

import pkg_resources

import pydesc
import pydesc.chemistry as monomer
import pydesc.contacts.contactmap as contactmap
import pydesc.contacts.contacts as contacts
import pydesc.structure as structure
import pydesc.util.typesdictionary as typesdictionary


def load_library(name):
    """
    Finds and loads a dynamic library. Returns an instance of ctypes.CDLL.


    If pydesc has been imported from a package it uses pkg_resources module to unpack and retrieve
    libraries. Otherwise it searches lib directory in the directory pydesc package is located.

    To provide portability it searches for files with extensions .dll, .so and .dylib.
    """

    # pkg_resorces causes spurious PyLint errors.
    extensions = [".dll", ".so", ".dylib"]

    try:
        req = pkg_resources.get_distribution(
            "pydesc"
        ).as_requirement()  # pylint: disable=E1103
        if pydesc.__file__ != pkg_resources.resource_filename(
            req, "pydesc/maps.py"
        ):  # pylint: disable=E1101
            req = None

    except pkg_resources.DistributionNotFound:
        req = None

    for ext in extensions:
        try:
            if req is not None:
                fname = pkg_resources.resource_filename(
                    req, "/lib/lib%s%s" % (name, ext)
                )  # pylint: disable=E1101
            else:
                fname = os.path.join(
                    os.path.dirname(os.path.dirname(pydesc.__file__)),
                    "lib/lib%s%s" % (name, ext),
                )
        except KeyError:
            continue

        if os.path.isfile(fname):
            return ctypes.CDLL(fname)

    raise Exception("Could not load lib%s." % name)


# This is not a constant.
libcydesc = load_library("cydesc")  # pylint: disable=C0103


def use_library(lib):
    """
    Decorator for CInDelMeta metaclass setting a C library to be searched for functions.
    """

    def class_wrapper(mcs):
        """
        Binds CInDelMeta to %s.
        """

        class WrappingClass(mcs):
            """ This class wraps CInDelMeta to change library object. """

            _lib = lib

        return WrappingClass

    class_wrapper.__doc__ = class_wrapper.__doc__ % lib._name  # pylint: disable=W0212
    return class_wrapper


class CInDelMeta(ctypes.Structure.__class__):
    # Metaclass needs access to protected members of a class.
    # pylint: disable=W0212

    """
    A metaclass applicable to subclasses of ctypes.Structure. It is responsible for
    calling C functions when creating and deleting instances, and for resising memory buffer
    to accomodate components visible only in C.

    It calls the following C functions:
        * sizeof_<<name>>
        * init_<<name>>
        * del_<<name>>

    where <<name>> is the name of a Python class being defined.

    By default this metaclass is bound to libcydesc. To change this use use_library(lib) decorator.
    """

    _lib = libcydesc

    def __call__(cls, *args, **kwargs):
        """
        Creates and initializes class instance.

        After executing superclass __call__ calls C function init_<<name>>.
        """
        res = super(CInDelMeta, cls).__call__(*args, **kwargs)

        init_func_name = "init_" + cls.__name__
        try:
            init_func = getattr(cls.__class__._lib, init_func_name)
        except AttributeError:
            pass
        else:
            init_func(ctypes.byref(res))

        return res

    def __new__(mcs, name, bases, nmspc):
        def del_dec(del_met):
            """ Decorator for __del__ method calling delete function in C """

            def wrapper(self):
                """%s
                Also calls %s in C.
                """
                del_met(self)
                if not getattr(self, "_already_deleted", False) and self._b_needsfree_:
                    del_func(ctypes.byref(self))
                    self._already_deleted = True

            wrapper.__doc__ = wrapper.__doc__ % (del_met.__doc__, del_func_name)
            return wrapper

        del_func_name = "del_" + name
        try:
            del_func = getattr(mcs._lib, del_func_name)
        except AttributeError:
            pass
        else:
            if "__del__" in nmspc:
                old_del = nmspc["__del__"]
            else:

                def old_del(self):
                    """ Calls __del__ in superclass (if exists). """
                    if hasattr(super(res, self), "__del__"):
                        # res refers to a class created by enclosing method
                        super(res, self).__del__()  # pylint: disable=E1003

            nmspc["__del__"] = del_dec(old_del)

        size = ctypes.sizeof(
            ctypes.Structure.__class__.__new__(mcs, name, bases, nmspc)
        )

        # setting the size of the structure
        sizeof_func_name = "sizeof_" + name
        try:
            sizeof_func = getattr(mcs._lib, sizeof_func_name)
        except AttributeError:
            csize = size
        else:
            csize = sizeof_func()

        if csize > size:
            nmspc["_fields_"].append(("_cdata_buf", ctypes.c_char * (csize - size)))

        res = ctypes.Structure.__class__.__new__(mcs, name, bases, nmspc)
        assert ctypes.sizeof(res) == csize

        return res


class CPoint(ctypes.Structure, metaclass=CInDelMeta):
    """
    Class for holding 3D coordinates. Corresponds to CPoint in cydesc.h
    """

    _fields_ = [("x", ctypes.c_float), ("y", ctypes.c_float), ("z", ctypes.c_float)]

    def __init__(self, coo):
        # pylint: disable=W0231, C0103
        # __init__ supplied by ctypes.Structure should not be called, if there
        # is an __init__ supplied in a subclass.
        # Single letter attribute names are acceptable in this context.
        """ Accepts an iterable containing at least three number (floats). """
        self.x = coo[0]
        self.y = coo[1]
        self.z = coo[2]

    def __iter__(self):
        """ Returns an iterator over three coordinates. """
        return iter([self.x, self.y, self.z])

    def __repr__(self):
        return "<CPoint (%.3f, %.3f, %.3f)" % (self.x, self.y, self.z)


class CMer(ctypes.Structure, metaclass=CInDelMeta):
    """
    Class for holding momomers. Corresponds to CMer in cydesc.h
    """

    _fields_ = [
        ("ind", ctypes.c_int),
        ("type", ctypes.c_int),
        ("next_ind", ctypes.c_int),
        ("n_points", ctypes.c_int),
        ("points", ctypes.POINTER(CPoint)),
        ("point_names", ctypes.POINTER(ctypes.c_char_p)),
        ("type_name", ctypes.c_char_p),
    ]

    # python types are converted into integers according to this dictionary
    types_dict = typesdictionary.TypesDictionary()
    types_dict[monomer.Residue] = 1
    types_dict[monomer.Nucleotide] = 2
    types_dict[monomer.Ligand] = 3
    types_dict[monomer.Ion] = 4

    def __init__(self, m):  # pylint: disable=W0231
        # __init__ supplied by ctypes.Structure should not be called, if there
        # is an __init__ supplied in a subclass.
        """ Accepts instances of monomer.Monomer. """
        self.ind = m.ind
        self.type = self.types_dict[type(m)]
        self.type_name = type(m).__name__

        atoms = m.representation
        points = [CPoint(list(a)) for a in atoms]
        point_names = [a.name for a in atoms]
        self.n_points = len(points)

        if m.next_mer:
            self.next_ind = m.next_mer.ind
        else:
            self.next_ind = 0

        # final conversion
        self.points = ctypes.cast(
            (CPoint * self.n_points)(
                *points
            ),  # converts <<points>> into a ctypes array
            ctypes.POINTER(CPoint),
        )  # casts an array of CPoint onto a POINTER (to CPoint)
        self.point_names = ctypes.cast(
            (ctypes.c_char_p * self.n_points)(*point_names),
            ctypes.POINTER(ctypes.c_char_p),
        )

    def __repr__(self):
        return "<CMer: ind: %d; type: %d(%s)>" % (self.ind, self.type, self.type_name)


class CSeg(ctypes.Structure, metaclass=CInDelMeta):
    """
    Class for holding segments. Corresponds to CSeg in cydesc.h
    """

    _fields_ = [("start", ctypes.c_int), ("end", ctypes.c_int)]

    def __iter__(self):
        """ Returns iterator over a part if start and end indices. """
        return iter((self.start, self.end))

    def __repr__(self):
        return "<CSeg: start: %d, end: %d>" % (self.start, self.end)


class CStructure(ctypes.Structure, metaclass=CInDelMeta):
    """
    Class for holding structures. Corresponds to CStructure in cydesc.h
    """

    _adjusted_number_ftype = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_int)

    _fields_ = [
        ("name", ctypes.c_char_p),
        ("n_monomers", ctypes.c_int),
        ("mers", ctypes.POINTER(CMer)),
        ("n_segs", ctypes.c_int),
        ("segs", ctypes.POINTER(CSeg)),
        ("adjusted_number_p", _adjusted_number_ftype),
    ]

    def __init__(self, struct):  # pylint: disable=W0231
        # __init__ supplied by ctypes.Structure should not be called, if there
        # is an __init__ supplied in a subclass.
        """ Accepts instances of structure.AbstractStructure. """
        mers = list(struct)
        cmers = list(map(CMer, mers))

        segs = []

        start = 0

        for (i, (mer1, mer2)) in enumerate(zip(mers, mers[1:])):
            if not mer1.next_mer == mer2:
                cmers[i].next_ind = 0
                segs.append((cmers[start].ind, cmers[i].ind))
                start = i + 1

        segs.append((cmers[start].ind, cmers[-1].ind))

        cmers[-1].next_ind = 0

        csegs = [CSeg(*s) for s in segs]

        try:
            self.name = struct.name
        except AttributeError:
            self.name = ""

        def adjusted_number(start, end):
            """ Returns an adjusted number of segments between given mers."""
            res = struct[start:end].adjusted_number()
            return res

        self.n_monomers = len(mers)
        self.monomers = ctypes.cast(
            (CMer * self.n_monomers)(*cmers), ctypes.POINTER(CMer)
        )
        self.n_segs = len(segs)
        self.segs = ctypes.cast((CSeg * self.n_segs)(*csegs), ctypes.POINTER(CSeg))

        self.adjusted_number_p = self._adjusted_number_ftype(adjusted_number)
        self.structure = struct

    def __repr__(self):
        return "<CStructure: %s; number of mers: %d, number of segments: %d>" % (
            self.name,
            self.n_monomers,
            self.n_segs,
        )


class CContact(ctypes.Structure, metaclass=CInDelMeta):
    """
    A C equivalent of a contact (structure.Contact).
    """

    _fields_ = [
        ("mer1", ctypes.c_int),
        ("mer2", ctypes.c_int),
        ("val", ctypes.c_int),
        ("type", ctypes.c_int),
    ]

    def __init__(self, mer1, mer2, val, type_=0):  # pylint: disable=W0231
        # __init__ supplied by ctypes.Structure should not be called, if there
        # is an __init__ supplied in a subclass.
        """
        Accepts PyDesc indices of mers, contact value and type.

        Type does not exist in PyDesc directly but may be used to transfer information to CyDesc.
        """

        self.mer1 = mer1
        self.mer2 = mer2
        self.val = val
        self.type = type_

    def __repr__(self):
        return "<CContact: mer1: %d, mer2: %d val: %d type: %d>" % (
            self.mer1,
            self.mer2,
            self.val,
            self.type,
        )


class CContactMap(ctypes.Structure, metaclass=CInDelMeta):
    """
    A C equivalent of a contactmap.ContactMap.

    Can be used along with a CStructure to to convey a descriptor (structure.AbstractDescriptor).
    """

    _fields_ = [
        ("structure", ctypes.POINTER(CStructure)),
        ("n_contacts", ctypes.c_int),
        ("contacts", ctypes.POINTER(CContact)),
    ]

    def __init__(self, struct, contact_map):  # pylint: disable=W0231
        # __init__ supplied by ctypes.Structure should not be called, if there
        # is an __init__ supplied in a subclass.
        """ Accepts instances of contactmap.ContactMap. """

        contacts_list = []

        for mer1, mer_dict in list(contact_map.contacts.items()):
            for mer2, val in list(mer_dict.items()):
                contacts_list.append(CContact(mer1, mer2, val))

        self.structure = ctypes.pointer(struct)
        self.n_contacts = len(contacts_list)

        # final conversion
        self.contacts = ctypes.cast(
            (CContact * self.n_contacts)(
                *contacts_list
            ),  # converts <<contacts>> into a ctypes array
            ctypes.POINTER(CContact),
        )  # casts an array of CContact onto a POINTER (to CContact)

    def __repr__(self):
        return "<CContactMap for %s:: n_contacts:%d>" % (
            repr(self.structure),
            self.n_contacts,
        )


class CElement(ctypes.Structure, metaclass=CInDelMeta):
    """
    A C equivalent of a structure.Element.
    """

    _fields_ = [
        ("center", ctypes.c_int),
        ("start", ctypes.c_int),
        ("end", ctypes.c_int),
        ("type", ctypes.c_int),
        ("status", ctypes.c_int),
        ("optional", ctypes.c_int),
    ]

    # python types are converted into integers according to this dictionary
    types_dict = {structure.ElementChainable: 1, structure.ElementOther: 2}

    def __init__(self, element):  # pylint: disable=W0231
        # __init__ supplied by ctypes.Structure should not be called, if there
        # is an __init__ supplied in a subclass.
        """ Accepts instances of structure.Element.

        This method always sets optional attribute to False, since this information
        is not available in Element object itself.

        """

        self.center = element.central_monomer.ind
        self.start = element[0].ind
        self.end = element[-1].ind
        self.type = self.types_dict[type(element)]
        self.optional = 0

    def __repr__(self):
        return "<CElement center: %d type: %d>" % (self.center, self.type)


class CDescriptor(ctypes.Structure, metaclass=CInDelMeta):
    """
    A C equivalent of a structure.Descriptor.

    This structure is in fact a tuple of CStructure, CContactMap, and an array
    of elements. This is due to the obvious fact, that C doesn't support inheritance.
    """

    _fields_ = [
        ("structure", ctypes.POINTER(CStructure)),
        ("contact_map", ctypes.POINTER(CContactMap)),
        ("n_elements", ctypes.c_int),
        ("elements", ctypes.POINTER(CElement)),
        ("central_element", ctypes.c_int),
    ]

    # python types are converted into integers according to this dictionary
    types_dict = {structure.ProteinDescriptor: 1, structure.NucleotideDescriptor: 2}

    def __init__(self, descriptor):  # pylint: disable=W0231
        # __init__ supplied by ctypes.Structure should not be called, if there
        # is an __init__ supplied in a subclass.
        """ Accepts instances of structure.AbstractDescriptor. """

        contact_map = contactmap.ContactMap(
            descriptor, contacts.DescriptorCriterion(descriptor)
        )
        contact_map.calculate_contacts()

        struct = CStructure(descriptor)
        ccmap = CContactMap(struct, contact_map)

        # elements = [CElement(structure.Element.build(descriptor[k])) for k, v in sorted(contact_map.contacts.items()) if len(v) > 0]
        elements = list(map(CElement, descriptor.elements))

        for element, celement in zip(descriptor.elements, elements):
            celement.optional = (
                1
                if descriptor.elements_values[element.central_monomer.ind].optional
                else 0
            )

        self.structure = ctypes.pointer(struct)
        self.contact_map = ctypes.pointer(ccmap)
        self.n_elements = len(elements)
        self.elements = ctypes.cast(
            (CElement * self.n_elements)(*elements), ctypes.POINTER(CElement)
        )
        self.central_element = descriptor.central_element.central_monomer.ind

    def __repr__(self):
        return "<CDescriptor for %s#%d:: n_elements:%d>" % (
            repr(self.structure),
            self.central_element,
            self.n_elements,
        )
