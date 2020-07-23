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

"""Auxiliary classes responsible for creation of different AtomSet subclasses."""

from abc import ABCMeta
from abc import abstractmethod

import numpy

from pydesc.config import ConfigManager
from pydesc.chemistry.base import Atom
from pydesc.chemistry.base import AtomProxy
from pydesc.chemistry.base import AtomSet
from pydesc.chemistry.full_atom import MonoatomicIon
from pydesc.chemistry.full_atom import Compound
from pydesc.chemistry.full_atom import Nucleotide
from pydesc.chemistry.full_atom import Residue
from pydesc.numberconverter import PDBid
from pydesc.warnexcept import IncompleteParticle
from pydesc.warnexcept import WarnManager
from pydesc.warnexcept import WrongAtomDistances
from pydesc.warnexcept import WrongAtomSetType


class AtomSetFactory(metaclass=ABCMeta):
    """Abstract factory class for AtomSet subclass instances.

    Subclasses produces AtomSet instances from different types of input data.

    Args:
        classes: list of subclasses of AtomSet to be used.

    """

    def __init__(self, classes=None):
        if classes is None:
            classes = [Residue, Nucleotide, MonoatomicIon, Compound]
        self.classes = classes

    @abstractmethod
    def create(self, *args):
        pass

    @property
    def chainable(self):
        return [i for i in self.classes if i.is_chainable()]

    @property
    def other(self):
        return [i for i in self.classes if not i.is_chainable()]

    def _create_possible_instances(self, base_atoms, warnings_, classes):
        """Return dictionary of different instance AtomSet subclasses as values and
        their classes as keys.

        Args:
            base_atoms: AtomSet instance initialized with some atoms.
            warnings_: context manager catching warnings.
            classes: list of AtomsSet subclasses to be tried out.

        Returns:
            tuple: first element is map of classes (keys) to instances (values);
            second is warning manager.

        """
        base_data = self.unpack_base(base_atoms)
        atoms = {}
        for subclass in classes:
            try:
                with warnings_(subclass):
                    instance = self._create_subclass_instance(subclass, base_data)
                    atoms[subclass] = instance
            except (IncompleteParticle, WrongAtomDistances, WrongAtomSetType):
                pass

        return atoms, warnings_

    @staticmethod
    def _create_subclass_instance(subclass, base_data):
        return subclass(*base_data)

    @staticmethod
    def unpack_base(base):
        """Return structure, PyDesc index, name, chain and atoms from given base"""
        return base.ind, base.name, base.chain, base.atoms


class BioPythonAtomSetFactory(AtomSetFactory):
    """AtomSet factory producing sets of atoms from parsed BioPython objects."""

    def create(
        self,
        pdb_residue,
        structure_obj=None,
        warn_in_place=True,
        warnings_=None,
        base=None,
    ):
        """Create all possible subclasses of AtomSet from given BioPython residue.

        Args:
            pdb_residue: instance of BioPython's Residue.
            structure_obj: host structure with attribute "converter" set.
            warn_in_place(bool): determines if warnings should be raised
            instantaneous, or if the should be suppressed and put to WarnManager.
            warnings_: instance of WarnManager.
            base: AtomSet instance with attribute "atoms" set.

        Returns:
            tuple: first element is map of classes (keys) to instances (values);
            second is warning manager.

        """
        name = self.get_pdb_residue_name(pdb_residue)
        if name in ConfigManager.chemistry.solvent:
            return None, None

        if warnings_ is None:
            warnings_ = WarnManager(pdb_residue)

        try:
            ind = structure_obj.converter.get_ind(
                PDBid.create_from_pdb_residue(pdb_residue)
            )
        except (AttributeError, KeyError):
            ind = None

        if base is None:
            base = AtomSet(ind, *self.unpack_pdb_residue(pdb_residue, name))

        atoms_dct, warnings_ = self._create_possible_instances(
            base, warnings_, self.classes
        )
        if warn_in_place:
            for class_ in self.classes:
                warnings_.raise_all(class_)

        atoms_dct[AtomSet] = base
        return atoms_dct, warnings_

    def unpack_pdb_residue(self, pdb_residue, name=None):
        """Return important data from pdb_residue.

        Args:
            pdb_residue: instance of BioPython's Residue.
            name(str): name of AtomSet (Residue in BioPython). None by default.

        Returns:
            : name(str), chain(str) and atoms(dict: name(str) -> atom instances)

        """
        if name is None:
            name = self.get_pdb_residue_name(pdb_residue)
        chain = pdb_residue.get_full_id()[2]
        crt = self.create_atom_from_bio_atom
        atoms = {
            pdb_atom.get_fullname().strip(): crt(pdb_atom) for pdb_atom in pdb_residue
        }
        return name, chain, atoms

    @staticmethod
    def create_atom_from_bio_atom(pdb_atom):
        """Return Atom instance from given Bio.Atom instance."""
        coords = numpy.array(pdb_atom.get_coord())
        serial_number = pdb_atom.get_serial_number()
        element = pdb_atom.element
        return Atom(coords, element, serial_number)

    @staticmethod
    def get_pdb_residue_name(pdb_residue):
        """Get residue name from given *pdb_residue* (Bio.PDBResidue instance)."""
        return pdb_residue.get_resname().strip()


class MDTrajAtomSetFactory(AtomSetFactory):
    """AtomSet factory producing sets of atoms consisting of AtomProxy objects,
    which are proxies of loaded trajectory."""

    def create(self, atoms_set, trajectory):
        """Create set of atoms linked with given pydesc Trajectory instance."""
        *base_data, atoms = self.unpack_base(atoms_set)
        atoms_proxies = {}
        for atom_name, atom in atoms.items():
            new_atom = AtomProxy(atom, trajectory)
            atoms_proxies[atom_name] = new_atom
        new_base = *base_data, atoms_proxies
        proxy_atoms = self._create_subclass_instance(type(atoms_set), new_base)
        return proxy_atoms


class CopyingFactor(AtomSetFactory):
    """AtomSet factory producing sets of atoms from already existing ones."""

    def create(self, atoms_set):
        """Return copy of given set of atoms.
        
        Args:
            atoms_set: instance of AtomsSet subclass.

        Returns:
            : copy of given argument.

        """
        *base_data, atoms = self.unpack_base(atoms_set)
        copied_atoms = {k: v.copy() for k, v in atoms.items()}
        base_data += (copied_atoms,)
        atoms_set = self._create_subclass_instance(type(atoms_set), base_data)
        return atoms_set
