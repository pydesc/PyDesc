from copy import deepcopy

import numpy

from pydesc.config import ConfigManager
from pydesc.numberconverter import PDBid
from pydesc.warnexcept import IncompleteParticle
from pydesc.warnexcept import WarnManager
from pydesc.warnexcept import WrongAtomDistances
from pydesc.warnexcept import WrongMerType
from pydesc.mers.base import Atom, AtomProxy
from pydesc.mers.base import Mer
from pydesc.mers.full_atom import Ion
from pydesc.mers.full_atom import Ligand
from pydesc.mers.full_atom import Nucleotide
from pydesc.mers.full_atom import Residue
from abc import ABCMeta, abstractmethod


class MerFactory(metaclass=ABCMeta):
    """Abstract factory class for Mer subclass instances."""

    def __init__(self, classes=None):
        """Mer factory initializer.

        Argument:
        classes -- list of classes to be used.

        This class tries to build mers classes from given.
        """
        if classes is None:
            classes = [Residue, Nucleotide, Ion, Ligand]
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

    @classmethod
    def copy_mer(cls, mer):
        """Return copy of given mer.

        Argument:
        mer -- mer subclass instance.
        """
        base_data = cls.unpack_base(mer)
        base_data = base_data[:-1] + (deepcopy(base_data[-1]),)
        mer = cls._create_mer_of_type(type(mer), base_data)
        return mer

    def _create_possible_monomers(self, base_monomer, warnings_, classes):
        """Return dictionary of different monomer types as values and
        subclasses of MonomerChainable and MonomerOther
        as keys.

        Arguments:
        mers.
        base_monomer -- an instance of Monomer class containing atoms from
        pdb_residue.
        warnings_ -- context manager for catching warnings.
        classes -- list of classes to try to initialize.
        """
        base_data = self.unpack_base(base_monomer)
        mers = {}
        for monomer_type in classes:
            try:
                with warnings_(monomer_type):
                    mers[monomer_type] = self._create_mer_of_type(
                        monomer_type, base_data
                    )
            except (IncompleteParticle, WrongAtomDistances, WrongMerType):
                pass

        return mers, warnings_

    @staticmethod
    def _create_mer_of_type(monomer_type, base_data):
        """Return monomer of given type based on given base data dictionary.

        Arguments:
            monomer_type -- monomer subclass.
            base_data -- dict of base data.
        """
        return monomer_type(*base_data)

    @staticmethod
    def unpack_base(base):
        """Return structure, PyDesc index, name, chain and atoms from given
        base (monomer.Monomer instance)."""
        return base.structure, base.ind, base.name, base.chain, base.atoms


class BioPythonMerFactory(MerFactory):
    def create(
        self,
        pdb_residue,
        structure_obj=None,
        warn_in_place=True,
        warnings_=None,
        base=None,
    ):
        """Create all possible mer subclass from given BioPython residue.

        Returns dict of created classes and WarnManager ready to raise warnings.
        """
        name = self.get_pdb_residue_name(pdb_residue)
        if name in ConfigManager.mers.solvent:
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
            base = Mer(structure_obj, ind, *self.unpack_pdb_residue(pdb_residue, name))

        mers, warnings_ = self._create_possible_monomers(base, warnings_, self.classes)
        if warn_in_place:
            for class_ in self.classes:
                warnings_.raise_all(class_)

        mers[Mer] = base
        return mers, warnings_

    def unpack_pdb_residue(self, pdb_residue, name=None):
        """Return important data from pdb_residue.

        Argument:
        pdb_residue -- instance of Bio.PDB.PdbResidue.
        name -- str; residue name; None by default. If so name is taken from
        pdb_residue with get_pdb_residue_name method.

        Returns tuple of name, chain name and dict of atoms.
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
        """Get residue name from given *pdb_residue* (Bio.PDBResidue
        instance)."""
        return pdb_residue.get_resname().strip()


class MDTrajMerFactory(MerFactory):

    """Mer factory responsible for creating mers consisting of AtomProxy objects."""

    def create(self, mer, trajectory):
        """Create mer linked with given pydesc Trajectory instance."""
        *base_data, atoms = self.unpack_base(mer)
        atoms_proxies = {}
        for atom_name, atom in atoms.items():
            new_atom = AtomProxy(atom, trajectory)
            atoms_proxies[atom_name] = new_atom
        new_base = *base_data, atoms_proxies
        proxy_mer = self._create_mer_of_type(type(mer), new_base)
        return proxy_mer
