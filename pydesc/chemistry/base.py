# Copyright 2019 Tymoteusz Oleniecki
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

"""Base classes for mers, ligands and generic group of atoms."""

from functools import wraps

import numpy
import scipy.linalg

import pydesc.geometry
import pydesc.geometry
from pydesc.chemistry import ConfigManager
from pydesc.warnexcept import IncompleteParticle
from pydesc.warnexcept import UnknownParticleName
from pydesc.warnexcept import WrongAtomDistances
from pydesc.warnexcept import warn

norm = scipy.linalg.get_blas_funcs("nrm2")
NotSet = object()


def register_pseudoatom(method):
    """Property with able to cache results in instances "pseudoatoms" attribute.

    Args:
        method: method to be turned into property.

    Returns:
        : property

    """

    @wraps(method)
    def get_cached_or_calculate(self):
        cache = self.pseudoatoms
        property_name = method.__name__
        try:
            return cache[property_name]
        except KeyError:
            value = method(self)
            return cache.setdefault(property_name, value)

    return property(get_cached_or_calculate)


def register_dynamic_feature(method):
    """Property with able to cache results in instances "dynamic_features" attribute.

    Args:
        method: method to be turned into property.

    Returns:
        : property

    """

    @wraps(method)
    def get_cached_or_calculate(self):
        cache = self.dynamic_features
        property_name = method.__name__
        try:
            return cache[property_name]
        except KeyError:
            value = method(self)
            return cache.setdefault(property_name, value)

    return property(get_cached_or_calculate)


class Atom(pydesc.geometry.Coord):
    """Point representation of atoms.

    Args:
        coords(array): atom's coordinates.
        element(str): name of atom's element.
        serial_number(int): atom's serial number.
        occupancy(float): atom's occupancy value.
        b_factor(float): B-factor value.

    """

    def __init__(
        self, coords, element, serial_number=None, occupancy=0.0, b_factor=0.0
    ):
        super().__init__(numpy_vec=coords)
        self.element = element
        self.occupancy = occupancy
        self.b_factor = b_factor
        self.serial_number = serial_number

    def __repr__(self):
        return "<Atom at %s>" % " ".join(["%.2f" % i for i in self.vector])

    def copy(self):
        """Return copy of self"""
        klass = type(self)
        data = (
            self.vector,
            self.element,
            self.serial_number,
            self.occupancy,
            self.b_factor,
        )
        return klass(*data)


class AtomProxy(pydesc.geometry.Coord):
    """Proxy implementing Atom interface, but taking coords from given Trajectory
    object."""

    def __init__(self, atom, trajectory):
        self.atom = atom
        self.trajectory = trajectory

    def copy(self):
        """Cast self to atom."""
        return self.cast_to_atom()

    def cast_to_atom(self):
        """Return new Atom instance holding state of represented trajectory atom."""
        atom_copy = self.atom.copy()
        atom_copy.vector = self.vector
        return atom_copy

    @property
    def vector(self):
        """Current coordinates."""
        return self.trajectory.get_atom_coords(self.atom)

    def __repr__(self):
        serial = self.atom.serial_number
        return f"<AtomProxy to #{serial} from {repr(self.trajectory)}>"


class Pseudoatom(pydesc.geometry.Coord):
    """Auxiliary points acting as atoms in some context.

    Takes coordinates as 3 separate floats or as single vector of floats.

    Args:
        x(float): x coordinate.
        y(float): y coordinate.
        z(float): z coordinate.
        numpy_vec(array): coordinates.
        name(str): pseudoatom name.

    """

    def __repr__(self):
        name = (" " + self.name) or ""
        coords = "%f %f %f" % tuple(self.vector)
        return "<Pseudoatom%s at %s>" % (name, coords)

    def __init__(self, x=0.0, y=0.0, z=0.0, numpy_vec=None, name=None):
        self.name = name
        super().__init__(x, y, z, numpy_vec)


class AtomSet:
    """Abstract class, representation of mers and particles present in molecular
    structures.

    Args:
        ind(int): id of set of atoms.
        name(str): name of set of atoms.
        chain(str): name of chain this set of atoms belongs to.
        atoms(dict): map of atoms names (str) to Atom instances.

    """

    @classmethod
    def reset_config_cache(cls):
        """Resets cache of configuration settings in this class and all subclasses.
        Should be called after relevant changes in ConfigManager."""
        cls._config_cache = {}
        for sub in cls.__subclasses__():
            sub.reset_config_cache()

    @classmethod
    def get_config(cls, prop_name):
        """Return given setting from configuration manager."""

        try:
            return cls._config_cache[prop_name]
        except KeyError:
            res = cls._get_config(prop_name)
            cls._config_cache[prop_name] = res
            return res

    @classmethod
    def _get_config(cls, prop_name):
        """Return given setting from configuration manager."""
        try:
            cls_name = cls.__name__.lower()
            branch = ConfigManager.chemistry  # pylint: disable=no-member
            if cls_name != "atomset":
                branch = getattr(branch, cls_name)

            res = getattr(branch, prop_name)
        except AttributeError:
            if issubclass(cls.__base__, AtomSet):
                res = cls.__base__.get_config(prop_name)
            else:
                raise

        return res

    @staticmethod
    def is_chainable():
        """Return False."""
        return False

    def __init__(self, ind, name, chain, atoms):
        self.name = name
        self.chain = chain
        self.ind = ind
        self.atoms = atoms
        self.pseudoatoms = {}
        self.dynamic_features = {}

        self._ss = "="

    def __len__(self):
        return len(self.atoms)

    def __repr__(self):
        try:
            return "<%s: %s no. %i>" % (self.__class__.__name__, self.name, self.ind)
        except (TypeError, KeyError):
            return "<%s: %s>" % (self.__class__.__name__, self.name)

    def __iter__(self):
        """Return iterator over atoms (not pseudoatoms)."""
        return iter(self.atoms.values())

    def __getattr__(self, name):
        name = name.lstrip()
        try:
            return self.atoms[name]
        except KeyError:
            repr_ = self.ind or id(self)
            class_name = type(self).__name__
            raise AttributeError(f"{class_name} {repr_} has no attribute {name}.")

    def has_bond(self, atoms_set):
        """Return False. Meant to be overwritten."""
        return False

    def reset_dynamic_cache(self):
        """Reset pseudoatoms and dynamic features."""
        self.dynamic_features = {}
        self.pseudoatoms = {}

    @register_pseudoatom
    def rc(self):
        """Geometrical center of AtomSet."""
        coordinates = [a.vector for a in self]

        vector = numpy.average(coordinates, 0)
        return Pseudoatom(numpy_vec=vector, name="rc")

    @property
    def representation(self):
        """Return PyDesc representation as list of atoms."""
        return [getattr(self, indicator) for indicator in self.get_config("indicators")]


class Mer(AtomSet):
    """Extension of AtomSet representing mers of biopolymers."""

    @staticmethod
    def is_chainable():
        """Return True."""
        return True

    def __init__(self, ind, name, chain, atoms):
        AtomSet.__init__(self, ind, name, chain, atoms)

        try:
            if self.get_config("check_distances"):
                backbone_atoms = dict(
                    (atom_name, None) for atom_name in self.get_config("backbone_atoms")
                )
                for atom_pair in self.get_config("crucial_atom_distances"):
                    self._check_distance(backbone_atoms, *atom_pair)
            self._check_backbone_atoms()
        except (AttributeError, KeyError):
            data = type(self).__name__, self.ind
            msg = "Backbone atoms lacking, unable to create %s from residue %s"
            raise IncompleteParticle(msg % data)

        self.next_mer = None
        self.prev_mer = None

    def _check_backbone_atoms(self):
        tuple(self.iter_bb_atoms())

    def _check_distance(self, atoms, name_1, name_2, min_dist, max_dist):
        """Raises WrongAtomDistances if atoms distance doesn't meet class criteria.

        Args:
            atoms: map atom names to Atom instances.
            name_1(str): 1st atom name.
            name_2(str): 2nd atom name.
            min_dist(float): min acceptable distance.
            max_dist(float): max acceptable distance.

        """
        distance = norm(atoms[name_1].vector - atoms[name_2].vector)
        if not min_dist <= distance <= max_dist:
            raise WrongAtomDistances(name_1.strip(), name_2.strip(), self)

    def has_bond(self, atoms_set):
        """Return True if given set of mers comes after current in biopolymer chain.

        Depends on setting ConfigManager.chemistry.mer_acceptable_distance.

        Args:
            atoms_set: instance of Mer class.

        Returns:
            bool: True if *atoms_set* is next to *self*.

        """
        # TODO: should depend on class-dependedn setting (if on any)
        if type(atoms_set) != type(self):
            return False
        bb_atoms = self.get_config("backbone_atoms")
        last_atom = self.atoms[bb_atoms[-1]]
        next_atom = atoms_set.atoms[bb_atoms[0]]
        try:
            distance = (last_atom - next_atom).calculate_length()
            return distance <= ConfigManager.chemistry.mer_acceptable_distance
        except UnboundLocalError:
            return False

    def iter_bb_atoms(self):
        """Get iterator over backbone atoms (see configuration)."""
        bb_atoms = self.get_config("backbone_atoms")
        return iter([self.atoms[attr_name] for attr_name in bb_atoms])

    def iter_nbb_atoms(self):
        """Get iterator over non-backbone atoms (see configuration)."""
        bb_atoms = self.get_config("backbone_atoms")
        return iter(
            [
                atom
                for atom_name, atom in list(self.atoms.items())
                if atom_name not in bb_atoms
            ]
        )

    @property
    def seq(self):
        try:
            code = self.get_config("code")
        except KeyError:
            warn(UnknownParticleName(self))
            return "?"
        try:
            return code[self.name]
        except KeyError:
            pass
        code = self.get_config("additional_code")
        return code[self.name]

    @register_pseudoatom
    def rc(self):
        """Geometrical center of side chain (non-backbone atoms)."""
        non_backbone_coordinates = [a.vector for a in self.iter_nbb_atoms()]
        vector = numpy.average(non_backbone_coordinates, 0)
        return Pseudoatom(numpy_vec=vector, name="rc")


class Ligand(AtomSet):
    """Abstract class, representation for ligands."""

    def __init__(self, ind, name, chain, atoms):
        AtomSet.__init__(self, ind, name, chain, atoms)
