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
Classes that represents mers present in representations of (sub)structures.

created: 11.07.2013 - 31.07.2013, Tymoteusz 'hert' Oleniecki
"""

import pydesc.geometry
import pydesc.selection
import numpy
import operator
import itertools

import scipy.linalg
norm = scipy.linalg.get_blas_funcs('nrm2')

from pydesc.config import ConfigManager as ConfigManager
from pydesc.numberconverter import PDB_id
from abc import ABCMeta
from pydesc.warnexcept import warn as warn
from pydesc.warnexcept import WarnManager as WarnManager
from pydesc.warnexcept import IncompleteParticle as IncompleteParticle
from pydesc.warnexcept import IncompleteChainableParticle as IncompleteChainableParticle
from pydesc.warnexcept import WrongAtomDistances as WrongAtomDistances
from pydesc.warnexcept import UnknownParticleName as UnknownParticleName
from pydesc.warnexcept import NoConfiguration as NoConfiguration
from pydesc.warnexcept import Info
try:
    import prody
except ImportError:
    warn(Info("No module: prody"))

# pylint: disable=no-member
ConfigManager.new_branch("monomer")
ConfigManager.monomer.set_default("monomer_acceptable_distance", 2.0)
ConfigManager.monomer.set_default("solvent", ['HOH'])
ConfigManager.monomer.new_branch("nucleotide")
ConfigManager.monomer.new_branch("residue")
ConfigManager.monomer.new_branch("monomerchainable")
ConfigManager.monomer.new_branch("ion")
ConfigManager.monomer.new_branch("ligand")
ConfigManager.monomer.set_default("backbone_atoms", ())
ConfigManager.monomer.monomerchainable.set_default("check_distances", False)
ConfigManager.monomer.residue.set_default("residue_code", {
    'ILE': 'I', 'GLN': 'Q',
    'GLX': 'Z', 'GLY': 'G',
    'GLU': 'E', 'CYS': 'C',
    'HIS': 'H', 'SER': 'S',
    'LYS': 'K', 'PRO': 'P',
    'ASX': 'B', 'ASN': 'N',
    'VAL': 'V', 'THR': 'T',
    'ASP': 'D', 'TRP': 'W',
    'PHE': 'F', 'ALA': 'A',
    'MET': 'M', 'LEU': 'L',
    'ARG': 'R', 'TYR': 'Y'})
ConfigManager.monomer.residue.set_default("residue_additional_code", {
    'DNP': 'A', 'ABI': 'A', 'ALM': 'A', 'MAA': 'A', 'TIH': 'A', 'FLA': 'A', 'DAL': 'A', 'CSD': 'A',
    'BNN': 'A', 'HAC': 'A', 'PRR': 'A', 'AYA': 'A', 'CHG': 'A', 'DHA': 'A', 'TPQ': 'A', 'SEG': 'A',
    'DIV': 'V', 'MVA': 'V', 'DVA': 'V',
    'BUG': 'L', 'DLE': 'L', 'CLE': 'L', 'NLN': 'L', 'NLE': 'L', 'NLP': 'L', 'MLE': 'L', 'LEF': 'L',
    'DIL': 'I', 'IIL': 'I',
    'DPR': 'P', 'HYP': 'P',
    'MSE': 'M', 'OMT': 'M', 'CXM': 'M', 'FME': 'M', 'MME': 'M',
    'DAH': 'F', 'PHI': 'F', 'DPN': 'F', 'HPQ': 'F', 'PHL': 'F',
    'LTR': 'W', 'TPL': 'W', 'DTR': 'W', 'TRO': 'W', 'HTR': 'W',
    'MSA': 'G', 'SAR': 'G', 'MPQ': 'G', 'GLZ': 'G', 'GSC': 'G', 'GL3': 'G', 'NMC': 'G',
    'DSN': 'S', 'SEL': 'S', 'SEP': 'S', 'SET': 'S', 'SAC': 'S', 'SVA': 'S', 'MIS': 'S', 'OAS': 'S',
    'TPO': 'T', 'ALO': 'T', 'DTH': 'T', 'BMT': 'T',
    'BCS': 'C', 'SOC': 'C', 'C5C': 'C', 'C6C': 'C', 'SCS': 'C', 'PEC': 'C', 'DCY': 'C', 'EFC': 'C',
    'SCY': 'C', 'SMC': 'C', 'CSX': 'C', 'BUC': 'C', 'CSO': 'C', 'PR3': 'C', 'CCS': 'C', 'CEA': 'C', 'CME': 'C',
    'CSP': 'C', 'CSS': 'C', 'CSW': 'C', 'CY1': 'C', 'CY3': 'C', 'CYG': 'C', 'CYM': 'C', 'CYQ': 'C', 'SCH': 'C',
    'SHC': 'C', 'OCS': 'C', 'CAS': 'C',
    'TYQ': 'Y', 'TYS': 'Y', 'TYB': 'Y', 'STY': 'Y', 'DTY': 'Y', 'IYR': 'Y', 'PAQ': 'Y', 'TYY': 'Y',
    'PTR': 'Y', 'TYI': 'Y',
    'MEN': 'N',
    'DGN': 'Q', 'MGN': 'Q',
    '2AS': 'D', 'ASB': 'D', 'DAS': 'D', 'ASK': 'D', 'ASL': 'D', 'ASQ': 'D', 'BHD': 'D', 'ASA': 'D',
    'DSP': 'D',
    '5HP': 'E', 'CGU': 'E', 'DGL': 'E', 'GMA': 'E', 'GGL': 'E', 'PCA': 'E',
    'DLY': 'K', 'LYM': 'K', 'LLY': 'K', 'LYZ': 'K', 'KCX': 'K', 'LLP': 'K', 'TRG': 'K', 'SHR': 'K',
    'ALY': 'K',
    'ARM': 'R', 'ACL': 'R', 'HAR': 'R', 'HMR': 'R', 'AGM': 'R', 'DAR': 'R',
    'HIC': 'H', '3AH': 'H', 'NEM': 'H', 'NEP': 'H', 'DHI': 'H', 'MHS': 'H', 'HIP': 'H', })
ConfigManager.monomer.residue.set_default("backbone_atoms", ('N', 'CA', 'C'))
ConfigManager.monomer.residue.set_default("check_distances", False)
ConfigManager.monomer.residue.set_default("crucial_atom_distances", (('C', 'CA', 1.35, 1.71), ('CA', 'N', 1.35, 1.75)))
ConfigManager.monomer.residue.set_default("indicators", ('CA', 'cbx'))
ConfigManager.monomer.residue.set_default("old_cbx_calculation", False)
ConfigManager.monomer.residue.set_default("adjusted_segment_length", 18.0)
ConfigManager.monomer.nucleotide.set_default("nucleotide_code", {'  G': 'G', '  C': 'C', '  U': 'U', '  A': 'A', ' DG': 'G', ' DA': 'A', ' DT': 'T', ' DC': 'C'})
ConfigManager.monomer.nucleotide.set_default("backbone_atoms", ("P", "O5'", "C5'", "C4'", "C3'", "O3'"))
ConfigManager.monomer.nucleotide.set_default("ring_atoms", ("N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"))
ConfigManager.monomer.nucleotide.set_default("check_distances", False)
ConfigManager.monomer.nucleotide.set_default("crucial_atom_distances", (('P', "O5'", 1.54, 1.66), ("O5'", "C5'", 1.34, 1.54), ("C5'", "C4'", 1.44, 1.56), ("C4'", "C3'", 1.46, 1.58), ("C3'", "O3'", 1.37, 1.49)))
ConfigManager.monomer.nucleotide.set_default("indicators", ("C3'", 'P',  'ring_center'))
ConfigManager.monomer.set_default("moving_average", 3)
ConfigManager.monomer.ion.set_default("indicators", ("rc",))
ConfigManager.monomer.ion.set_default("radii", {'BE': 0.59,
                                         'BA': 1.49,
                                         'BI': 1.17,
                                         'BK': 1.1,
                                         'BR': 1.82,
                                         'RU': 0.82,
                                         'RE': 0.77,
                                         'TM': 1.17,
                                         'RA': 1.62,
                                         'RB': 1.66,
                                         'RH': 0.805,
                                         'P': 0.58,
                                         'GE': 0.87,
                                         'GD': 1.078,
                                         'GA': 0.76,
                                         'OS': 0.77,
                                         'C': 0.3,
                                         'HO': 1.041,
                                         'HF': 0.85,
                                         'HG': 1.33,
                                         'PR': 1.13,
                                         'PT': 0.94,
                                         'PU': 1.14,
                                         'PB': 1.33,
                                         'PA': 1.16,
                                         'PD': 1.0,
                                         'PO': 1.08,
                                         'PM': 1.11,
                                         'ZN': 0.88,
                                         'K': 1.52,
                                         'O': 1.26,
                                         'S': 1.7,
                                         'W': 0.8,
                                         'EU': 1.31,
                                         'ZR': 0.86,
                                         'ER': 1.03,
                                         'MG': 0.86,
                                         'MO': 0.83,
                                         'MN': 0.97,
                                         'AU': 1.51,
                                         'FR': 1.94,
                                         'FE': 0.92,
                                         'NI': 0.83,
                                         'NA': 1.16,
                                         'NB': 0.86,
                                         'ND': 1.43,
                                         'ES': 0.928,
                                         'NP': 1.24,
                                         'B': 0.41,
                                         'CO': 0.885,
                                         'CM': 1.11,
                                         'CL': 1.67,
                                         'CA': 1.14,
                                         'CF': 1.09,
                                         'CE': 1.15,
                                         'N': 1.32,
                                         'V': 0.93,
                                         'CS': 1.81,
                                         'CR': 0.94,
                                         'CU': 0.91,
                                         'SR': 1.32,
                                         'SI': 0.54,
                                         'SN': 0.83,
                                         'SM': 1.36,
                                         'SC': 0.885,
                                         'SB': 0.9,
                                         'SE': 1.84,
                                         'YB': 1.16,
                                         'DY': 1.21,
                                         'LA': 1.172,
                                         'F': 1.19,
                                         'LI': 0.9,
                                         'TL': 1.64,
                                         'LU': 1.001,
                                         'TH': 1.08,
                                         'TI': 1.0,
                                         'TE': 2.07,
                                         'TB': 1.063,
                                         'TC': 0.785,
                                         'TA': 0.86,
                                         'AC': 1.26,
                                         'AG': 1.29,
                                         'I': 2.06,
                                         'IR': 0.82,
                                         'AM': 1.4,
                                         'AL': 0.675,
                                         'AS': 0.72,
                                         'U': 1.165,
                                         'AT': 0.76,
                                         'IN': 0.94,
                                         'Y': 1.04,
                                         'CD': 1.09,
                                         'XE': 0.62})

ConfigManager.monomer.ligand.set_default("indicators", ("rc",))
ConfigManager.new_branch("structure_mon")
ConfigManager.structure_mon.set_default("simple_secondary_structure_code", {'H': 'H', 'B': 'E', 'E': 'E', 'G': 'H', 'I': 'H', 'T': 'C', 'S': 'C', '-': 'C', '=': '='})
# pylint: enable=no-member


class MonomerFactory(object):

    def __init__(self, classes=None):
        """Monomer factory initializator.

        Argument:
        classes -- list of classes to be used.

        This class tries to build monomers classes from given.
        """
        if classes is None:
            classes = [Residue, Nucleotide, Ion, Ligand]
        self.chainable = [i for i in classes if issubclass(i, MonomerChainable)]
        self.other = [i for i in classes if issubclass(i, MonomerOther)]

    def create_from_BioPDB(self, pdb_residue, structure_obj=None, warn_in_place=True, warnings_=None, base=None, dbg=False):
        """Class method, returns Monomer instances.

        Returns dictionary of different monomer types as values, calls _create_monomers to create actual objects.
        This method facilitates checks and routines common to all monomer creations.

        The returned dictionary contains two special entries:
            'warnings' - an instance of WarnManager storing eventual warnings.
            Monomer - an instance of Monomer class containing atoms from pdb_residue.

        Arguments:
        pdb_residue -- instance of BioPython Bio.PDB.Residue based on which monomer is created.
        structure_obj -- Structure instance to which the monomer belongs. Could be None for unbonded monomers.
        Initially set to None.
        warn_in_place -- True or False. Determines if warnings are to be raised immediately or returned as a result.
        The former forces constructores to raise warning immediately.
        The latter stores raised warnings in context manager delivered as value of 'warnings' key in returned dictionary.
        warnings_ -- context manager for catching warnings. Should be supplied when restarting monomer creation.
        base -- an instance of Monomer class containing atoms from pdb_residue. Should be supplied when restarting monomer creation.
        """

        name = self.get_pdb_residue_name(pdb_residue)
        if name in ConfigManager.monomer.solvent:
            return None     # We ignore solvent

        if warnings_ is None:
            warnings_ = WarnManager(pdb_residue)

        try:
            ind = self.structure.converter.get_ind(pydesc.numberconverter.PDB_id.from_pdb_residue(pdb_residue))
        except (AttributeError, KeyError):
            ind = None

        if base is None:
            base = Monomer(structure_obj, ind, *self.unpack_pdb_residue(pdb_residue, name))

        mers = self._create_monomers(pdb_residue, structure_obj, base, warnings_, self.chainable + self.other, dbg)
        if warn_in_place:
            for class_ in self.chainable + self.other:  #TODO do we really need to keep them separately?
                warnings_.raise_all(class_)
        else:
            mers['warnings'] = warnings_

        mers[Monomer] = base
        return mers

    def _create_monomers(self, pdb_residue, structure_obj, base_monomer, warnings_, classes, dbg=False):
        """Return dictionary of different monomer types as values and subclasses of MonomerChainable and MonomerOther as keys.

        Arguments:
        pdb_residue -- instance of BioPython Bio.PDB.Residue based on which monomer is created.
        structure_obj -- Structure instance to which the monomer belongs. Could be None for unbonded monomers.
        base_monomer -- an instance of Monomer class containing atoms from pdb_residue.
        warnings_ -- context manager for catching warnings.
        classes -- list of classes to try to initialize.
        """
        base_data = self.unpack_base(base_monomer)
        mers = {}
        for monomer_type in classes:
            try:
                with warnings_(monomer_type):
                    mers[monomer_type] = monomer_type(*base_data)
            except (IncompleteParticle, WrongAtomDistances, AttributeError, ValueError, KeyError) as e:
                # AttributeError is raised by nucleotides during Residue.__init__
                # KeyError is raised by different __init__s when needed Atoms/Pseudoatoms are absent
                pass

        return mers

    def unpack_pdb_residue(self, pdb_residue, name=None):
        """Return important data from pdb_residue.

        Argument:
        pdb_residue -- instance of Bio.PDB.PdbResidue.
        name -- str; residue name; None by default. If so name is taken from pdb_residue with get_pdb_residue_name method.

        Returns tuple of name, chain name and dict of atoms.
        """
        if name is not None:
            name = self.get_pdb_residue_name(pdb_residue)
        chain = pdb_residue.get_full_id()[2]
        atoms = {pdb_atom.get_fullname().strip(): Atom(pdb_atom) for pdb_atom in pdb_residue}
        return name, chain, atoms

    def unpack_base(self, base):
        return base.structure, base.ind, base.name, base.chain, base.atoms

    def get_pdb_residue_name(self, pdb_residue):
        return pdb_residue.get_resname().strip()


class Atom(pydesc.geometry.Coord):

    """Representation of atoms described in pdb files.

    Subclass of pydesc.geometry.Coord class.

    Attributes:
    name -- string containing atom name.
    element -- string, element name.
    pdb_atom -- instance of BioPython Atom class.
    """

    def __init__(self, pdb_atom):     # pylint:disable=super-init-not-called
        # there is no need to call dict.__init__
        """Atom constructor.

        Arguments:
        pdb_atom -- instance of BioPython Atom class.
        """
        self._vector = None
        self.vector = numpy.array(pdb_atom.get_coord())
        self.element = pdb_atom.element
        self.pdb_atom = pdb_atom
        self.prody_atom = None

    def __repr__(self):
        return "<Atom at %s>" % " ".join(["%.2f" % i for i in self.vector])

    @property
    def vector(self):
        """Returns Atoms vector."""
        return self._vector if self.prody_atom is None else self.prody_atom.getCoords()

    @vector.setter
    def vector(self, value):
        """Sets Atoms vector."""
        self._vector = value


class Pseudoatom(pydesc.geometry.Coord):

    """Representation of any point related to monomer other than atom.

    Subclass of pydesc.geometry.Coord class.

    Attributes:
    name -- string containing point name.
    """

    def __repr__(self):
        return "<Pseudoatom %s: %f %f %f>" % ((self.name,) + tuple(self.vector))

    def __init__(self, x=.0, y=.0, z=.0, numpy_vec=None, name=''):     # pylint:disable=super-init-not-called
        # there is no need to call dict.__init__
        """Pseudoatom constructor.

        Arguments:
        x, y, z -- pseudoatom coordinates.
        numpy_vec -- NumPy array containing coordinates (if provided XYZ are ignored). None by default.
        name -- string, pseudoatom name. None by default.
        owner -- instance od pydesc.monomer.Monomer subclass that contains pseudoatom. None by default.
        calc_method -- owner method to calculate pseudoatom coordinates in dynamic mode. None by default.

        Sets attribute 'dynamic' to False. If set to True - coordinates are calculated each time when readed.
        """
        self.name = name
        pydesc.geometry.Coord.__init__(self, x, y, z, numpy_vec)


class DynamicPropertiesDict(dict):

    """Class of dicts to store values that need to be recalculated for every frame of molecular dynamics trajectory."""

    def __init__(self, owner):
        """DynamicProprtyDict costructor.

        Argument:
        owner -- instance of owning mer.

        Sets attribut 'dynamic' to False.
        """
        self.owner = owner
        self.dynamic = False
        dict.__init__(self)

    def __getitem__(self, key):
        """Returns value of given key.

        Argument:
        key -- name of value to be returned or recalculated.

        Returns value of given key if 'dynamic' attrbute is set to False. Otherwise tries to return value of the key if it is not None. If value is None - forces owner to recalculate the value byt calling "calculate_<key>" method.
        """
        value = dict.__getitem__(self, key)
        if self.dynamic:
            if value is None:
                object.__getattribute__(self.owner, "calculate_%s" % key)()
                value = dict.__getitem__(self, key)
        return value

    def add(self, item):
        """Adds item to self with key equal to item name."""
        self[item.name] = item

    def reset_all_values(self):
        """Sets all values to None."""
        for key in self:
            self[key] = None

    def recalculate_all_values(self):
        """Forces recalculation of all values."""
        for key in self:
            object.__getattribute__(self.owner, "calculate_%s" % key)()


class Monomer(object):

    """Abstract class, representation of mers and particles present in molecular structures.

    Subclasses:
    MonomerChainable -- residues and nucleotides.
    MonomerOther -- ligands or their type: ions.
    """

    @classmethod
    def reset_config_cache(cls):
        """Resets cache of configuration settings in this class and all subclasses. Should be called after relevant changes in ConfigManager.
        """
        cls._config_cache={}
        for sub in cls.__subclasses__():
            sub.reset_config_cache()

    @classmethod
    def get_config(cls, prop_name):
        """Returns class configuration from confiuration manager.

        All data is cached in _config_cache class attribute.

        Argument:
        prop_name -- name of configuration to be returned.
        """

        try:
            return cls._config_cache[prop_name]
        except KeyError:
            res = cls._get_config(prop_name)
            cls._config_cache[prop_name] = res
            return res

    @classmethod
    def _get_config(cls, prop_name):
        """Returns class configuration from confiuration manager.

        Argument:
        prop_name -- name of configuration to be returned.
        """

        try:
            cls_name = cls.__name__.lower()

            branch = ConfigManager.monomer  # pylint: disable=no-member
            if cls_name != 'monomer':
                branch = getattr(branch, cls_name)

            res = getattr(branch, prop_name)
        except AttributeError:
            if issubclass(cls.__base__, Monomer):   # pylint: disable=no-member
                res = cls.__base__.get_config(prop_name)  # pylint:disable=no-member, protected-access
                # __base__ is not absent
                # protected access to superclass method
            else:
                raise

        return res


    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Monomer contructor.

        Arguments:
        structure_obj -- Structure in which current monomer is included.
        name -- str; mers name.
        chain -- str; chain name.
        atoms -- dict; dict of str names of atoms as keys and Atom instances as values.

        Sets attributes:
        name -- mer or ligand name, up to three letters, according to PDB file.
        structure -- the Structure instance to which the monomer belongs.
        my_chain -- character of the chain that the monomers belong to, according to PDB file.
        atoms - dict of atoms building current monomer repersented by Atom instances.
        ind -- PyDesc integer.
        pseudoatoms -- dict of Pseudoatoms.
        dynamic_properties -- dict of other geometrical properties like planes for cyclic chemical compounds.
        _ss -- secondary structure sign.
        """

        self.structure = structure_obj
        self.name = name
        self.chain = chain
        self.ind = ind
        self.atoms = atoms

        self.pseudoatoms = DynamicPropertiesDict(self)
        self.dynamic_properties = DynamicPropertiesDict(self)
        self._ss = '='

    def __len__(self):
        """Return sum of lenghts of monomer's atoms and pseudoatoms."""
        return len(list(iter(self)))

    def __repr__(self):
        try:
            return '<%s: %s no. %i, PDB: %s>' % (self.__class__.__name__, self.name, self.ind, str(self.get_pdb_id()))
        except (TypeError, KeyError):
            return '<%s: %s>' % (self.__class__.__name__, self.name)

    def __iter__(self):
        """Return monomer iterator.

        Monomer iterator iterates over its atoms and pseudoatoms dictionaries.
        """
        return iter([self.atoms[atom] for atom in sorted(self.atoms)] + [self.pseudoatoms[point] for point in sorted(self.pseudoatoms)])

    def __getattr__(self, name):
        """Returns proper attribute value.

        Argument:
        name -- string, attribute name.
        """
        name = name.lstrip()
        try:
            return self.atoms[name]
        except KeyError:
            try:
                return self.pseudoatoms[name]
            except (AttributeError, KeyError):
                repr_ = self.ind if self.ind is not None else str(self)
                raise AttributeError("Monomer %s has no attribute %s" % (repr_, name))

    def __getitem__(self, name):
        """Deprecated method. Returns proper attribute value.

        Argument:
        name -- string, attribute name.
        getitem -- True by default, False if called by __getattr method.
        """
        warn(DeprecationWarning(
        """Atom eventually won't inherit from dict type, so avoid getting to attributes via getitem.
        Use getattr instead, e.g.
        instead of
        >>> print my_atom['rc']
        use
        >>> print my_atom.rc
        or access atoms or pseudoatoms dicts directly:
        >>> print my_atom.pseudoatoms['rc']
        """), 1)
        name = name.lstrip()
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            try:
                return object.__getattribute__(self, 'atoms')[name]
            except KeyError:
                try:
                    return self.pseudoatoms[name]
                except (AttributeError, KeyError):
                    repr_ = self.ind if self.ind is not None else str(self)
                    raise AttributeError("Monomer %s has no attribute %s" % (repr_, name))

    def _finalize(self):
        """Method called by structures to calculate and set attributes that need structural information to be calculated."""
        self.calculate_rc()

    def calculate_rc(self):
        """Sets Monomer's attribute rc (geometrical center).

        Adds pydesc.geometry.Coord instance representing the geometrical center of a mer to mers pseudoatoms dict.
        If possible, only sidechain atoms are taken into account.
        """
        nbbcrds = [a.vector for a in self.iter_atomsnbb()]

        if nbbcrds != []:
            vector = numpy.average(nbbcrds, 0)
        else:
            try:
                vector = self.ca.vector
            except AttributeError:
                vector = self.atoms['P  '].vector
        self.pseudoatoms['rc'] = Pseudoatom(numpy_vec=vector, name='rc')

    def iter_atoms(self):
        """Returns iterator that iterates over monomer's atoms."""
        return iter(self.atoms.values())

    def iter_atomsbb(self):
        """Returns iterator that iterates over monomer's backbone atoms."""
        return iter([])

    def iter_atomsnbb(self):
        """Returns iterator that iterates over monomer's all atoms except backbone."""
        return self.iter_atoms()

    @property
    def indicators(self):
        """Returns names of monomer attributes that store coordinates used in RMSD calculation."""
        return self.get_config('indicators')
        # temp = getattr(ConfigManager.monomer, str(type(self).__name__).lower())
        # return  getattr(temp, "indicators")

    @classmethod
    def seq_3to1(cls, seq):
        """Returns a one letter code for a given 3-letter code."""
        try:
            cls_name = cls.__name__.lower()
            code_dictionary = getattr(getattr(ConfigManager.monomer, cls_name), cls_name + "_code")     # pylint:disable=no-member
            try:
                additional_dictionary = getattr(getattr(ConfigManager.monomer, cls_name), cls_name + "_additional_code")     # pylint:disable=no-member
            except AttributeError:
                additional_dictionary = {}
            return code_dictionary[seq] if seq in code_dictionary else additional_dictionary[seq]
        except AttributeError:
            if issubclass(cls.__base__, Monomer):     # pylint:disable=no-member
                # ??? Monomer has no __base__
                return cls.__base__.seq_3to1(seq)     # pylint:disable=no-member
                # ??? same here
            raise AttributeError("No dictionary defined for class %s", str(cls))

    @classmethod
    def seq_1to3(cls, let):
        """Returns a three letter code for a given 1-letter code. In ambiguous cases the first matching code is returned."""
        try:
            cls_name = cls.__name__.lower()
            code_dictionary = getattr(getattr(ConfigManager.monomer, cls_name), cls_name + "_code")     # pylint:disable=no-member
            for seq3, seq1 in code_dictionary.items():
                if seq1 == let:
                    return seq3
            raise KeyError('Cannot translate %s to 3 letter code' % (cls_name + " symbol " + let,))
        except AttributeError:
            if issubclass(cls.__base__, Monomer):     # pylint:disable=no-member
                return cls.__base__.seq_1to3(let)     # pylint:disable=no-member
                # Monomer has __base__ attr
            raise AttributeError("No dictionary defined for class %s", str(cls))

    @property
    def seq(self):
        """Returns one letter code for mer if possible ("?" if name is unknown)."""
        try:
            return self.seq_3to1(self.name)
        except KeyError:
            warn(UnknownParticleName(self))
            return "?"
        except AttributeError:
            warn(NoConfiguration("class %s has no dictionary in configuration manager, thus '=' inserted into sequence. to turn this exception into harmless warning - set NoConfiguration in ConfigManager.warnings_and_exceptions.class_filters to 'ignore' or ;always'" % self.__class__.__name__))
            return "="

    @property
    def seq3(self):         # ??? trzyliterowy kod to imie w wypadku aa, a reszta?
        """Returns mer three letter pdb name."""
        return self.name

    @property
    def representation(self):
        """Returns indicators of current monomer set in configuration manager."""
        return [getattr(self, indicator) for indicator in self.indicators]

    def get_representation(self):
        """Deprecated method. Returns self.representation."""
        warn(DeprecationWarning("Method get_representation is deprecated. Please, use property representation instead."), 1)
        return self.representation

    def get_pdb_id(self):
        """Returns pdb id if possible, otherwise returns None."""
        try:
            return self.structure.converter.get_pdb_id(self.ind)
        except AttributeError:
            return None

    @property
    def pid(self):
        return str(self.get_pdb_id())

    @property
    def secondary_structure(self):
        """Secondary structure obtained with DSSP for maternal structure.

        If DSSP was not found or secondary structure was not calculated - returns '=' sign.

        To calculate secondary structure - use maternal structure method set_secondary_structure.

        See Bio.PDB.DSSP documentary for information about code explanation.
        """
        return self._ss

    @property
    def simple_secondary_structure(self):
        """Secondary structure in simple 3-letter code for secondary structures.

        H -- helix
        E -- extended strand
        C -- coil
        """
        temp = ConfigManager.structure_mon.simple_secondary_structure_code  # pylint:disable=no-member
        # configuration manager is dynamic with member that cannot be recognized by pylint
        return temp[self._ss]

    def select(self):
        """Returns selection - set containing only self."""
        return pydesc.selection.Set([self.structure.converter.get_pdb_id(self.ind)])
        # ??? monomery bez struktury?

    def is_next(self, monomer):
        """Deprecated method. Returns True if given mer follows current mer."""
        warn(DeprecationWarning("is_next method is deprecated, use attribute next_monomer instead."))
        return self.next_monomer == monomer

    def is_prev(self, monomer):
        """Deprecated method. Returns True if given mer preceeds current mer."""
        warn(DeprecationWarning("is_prev method is deprecated, use attribute prevoius_monomer instead."))
        return self.previous_monomer == monomer


class MonomerChainable(Monomer):

    """Abstract class, representation of residue or nucleotide.

    Subclasses:
    Residue
    Nucleotide
    """

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Chainable monomer constructor.

        Extends superclass method.
        """
        Monomer.__init__(self, structure_obj, ind, name, chain, atoms)
        try:
            backbone_atoms = dict((atom_name, None) for atom_name in self.get_config('backbone_atoms'))
        except AttributeError:
            raise IncompleteParticle(self.pid)
        
        if self.get_config('check_distances'):
            for atom_pair in self.get_config("crucial_atom_distances"):
                self._check_distance(backbone_atoms, *atom_pair)
        self._asa = None

    def _has_bond(self, monomer):
        """Returns True if the Monomer is followed by a given Monomer.

        Argument:
        monomer -- MonomerChainable instance.

        Calculates distance between backbone atoms of Monomers. Returns True or False according to the configurable monomer_acceptable_distance.
        """
        if type(monomer) != type(self):
            return False
        last_atom = self.backbone[-1]
        next_atom = monomer.backbone[0]
        try:
            distance = (last_atom - next_atom).calculate_length()
            return distance <= ConfigManager.monomer.monomer_acceptable_distance     # pylint:disable=no-member
        except UnboundLocalError:
            return False

    def _check_name(self, afix):
        """Method that raises warning if unknown particle name was found in pdb file."""
        try:
            self.seq_3to1(self.name)
        except KeyError:
            warn(UnknownParticleName(self))

    def _check_distance(self, atoms, name_1, name_2, min_dist, max_dist):
        """Raises WrongAtomDistances if atoms distance doesn't meet class criteria.

        Arguments:
        atoms -- dictionary containig atom names as keys and Atom instance as values.
        name_1, name_2 -- 1st and 2nd atom names.

        NOTE: method requires attributes in congifuration manager. They must be integers and their names should match pattern:
        min_<first lower and stripped name>_<second lower and stripped name>_dist
        and
        max_<first lower and stripped name>_<second lower and stripped name>_dist
        names should be given in default order for strings.
        """
        distance = norm(atoms[name_1].vector - atoms[name_2].vector)
        if not min_dist <= distance <= max_dist:
            raise WrongAtomDistances(name_1.strip(), name_2.strip(), self)

    @property
    def next_monomer(self):
        """Property that returns monomer following current mer in its structure."""
        try:
            return self._next_monomer
        except AttributeError:
            return None

    @next_monomer.setter
    def next_monomer(self, value):
        """Property that returns monomer following current mer in its structure."""
        self._next_monomer = value     # pylint:disable=attribute-defined-outside-init
        # there is no other way to complet property setter

    @next_monomer.deleter
    def next_monomer(self):
        """Property that returns monomer following current mer in its structure."""
        del self._next_monomer

    @property
    def previous_monomer(self):
        """Property that returns monomer preceding current mer in its structure."""
        try:
            return self._previous_monomer
        except AttributeError:
            return None

    @previous_monomer.setter
    def previous_monomer(self, value):
        """Property that returns monomer preceding current mer in its structure."""
        self._previous_monomer = value     # pylint:disable=attribute-defined-outside-init
        # same as in next_monomer.setter

    @previous_monomer.deleter
    def previous_monomer(self):
        """Property that returns monomer preceding current mer in its structure."""
        del self._previous_monomer

    @property
    def backbone(self):
        """Returns backbone atoms"""
        bb_atoms = self.get_config('backbone_atoms')
        return [self.atoms[attr_name] for attr_name in bb_atoms]

    def iter_atomsbb(self):
        """Returns iterator that iterates over monomer's backbone atoms."""
        return iter(self.backbone)

    def iter_atomsnbb(self):
        """Returns iterator that iterates over monomer's all atoms except backbone."""
        bb_atoms = self.get_config('backbone_atoms')
        return iter([i for i in self.atoms.values() if i not in self.backbone])

    def adjusted_length(self):
        """Returns distance between backbone_average pseudoatoms of this and the next monomer or None if distance cannot be computed."""

        try:
            return abs(self.backbone_average - self.next_monomer.backbone_average)
        except AttributeError:
            return None


class Residue(MonomerChainable):

    """Representation of a residue."""

    @staticmethod
    def calculate_angles_static(structure_obj):
        """Calculates all torsion angles of residues in given structure.

        Argument:
        structure_obj -- instance of AbstractStructure subclass.

        Fills 'angles' property in all residues in given (sub)structure. Calculates them using numpy, much faster than non-static Residue method.
        """
        residues = [mer for mer in structure_obj if isinstance(mer, Residue)]
        nres = len(residues)
        if nres == 0: return
        n, ca, c = numpy.transpose(numpy.array([[a.vector for a in r.backbone] for r in residues]), (1,0,2))[[0,1,2]]
        pc = numpy.empty((nres, 3), dtype=numpy.float32)
        nn = numpy.empty((nres, 3), dtype=numpy.float32)

        pc[1:] = c[:-1]
        nn[:-1] = n[1:]

        no_prev = numpy.fromiter((r.previous_monomer is None for r in residues), dtype=bool)
        no_next = numpy.fromiter((r.next_monomer is None for r in residues), dtype=bool)

        pc[no_prev] = n[no_prev]
        nn[no_next] = c[no_next]

        cca = c - ca
        cnn = c - nn
        nca = n - ca
        npc = n - pc

        pl1 = numpy.cross(cca, cnn)    # vectors perpendicular to plane 1
        pl2 = numpy.cross(nca, cca)    # vectors perpendicular to plane 2
        pl3 = numpy.cross(npc, nca)    # vectors perpendicular to plane 3

        with numpy.errstate(divide='ignore',invalid='ignore'):
            pl1, pl2, pl3 = (pl / numpy.sqrt(numpy.einsum('ij,ij->i', pl, pl)).reshape(-1,1) for pl in (pl1, pl2, pl3))

        angs = []
        for planes, direction in (((pl1, pl2), -cca), ((pl2, pl3), nca)):
            cos = numpy.einsum('ij,ij->i', *planes)
            cpr = numpy.cross(*planes)
            sin = numpy.sqrt(numpy.einsum('ij,ij->i', cpr, cpr))
            sign = numpy.sign(numpy.einsum('ij,ij->i', direction, cpr))

            t2 = numpy.arctan2(sin, cos) * sign
            t1 = numpy.nan_to_num(t2)

            angs.append(t1)

        for res, (psi, phi) in zip(residues, zip(*angs)):
            res.dynamic_properties['angles'] = (psi, phi)

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Residue constructor.

        Arguments:
        pdb_residue -- BioPython Bio.PDB.Residue instance based on which the Residue is being created.
        structure_obj -- the Structure instance which the Residue belongs to.

        Raises Warning if a given pdb_residue does not contain proper atoms or if its atoms occur in wrong distances.
        Extended MonomerChainable method.
        See also config file docstring.

        Config parameteres in branch ConfigManager.monomer.residue:
        min_c_ca_dist
        max_c_ca_dist
        min_ca_n_dist
        max_ca_n_dist
        min_c_o_dist
        max_c_o_dist
        old_cbx_calculation -- True or False
        """
        MonomerChainable.__init__(self, structure_obj, ind, name, chain, atoms)
        self.calculate_cbx()
        self._check_name('residue')

    def _finalize(self):
        """Method called by structures to calculate and set attributes that need structural information to be calculated."""
        super(Residue, self)._finalize()
        self.calculate_backbone_average()

    def calculate_backbone_average(self):
        """Calculates coordinates of average ca pseudoatom and adds it to current residue pseudotaoms.

        Average ca is calculated as moving average for configurable number of residues around current residue.
        """
        steps = self.get_config('moving_average')     # pylint:disable=no-member
        if not steps % 2 == 1:
            raise ValueError("Wrong Number of steps for moving average. Configure correnct length with ConfigManager")
        average_ca = numpy.array(self.ca.vector)
        next_mer = last_mer = self
        cnt = 1
        try:
            for step in range(steps / 2):     # pylint:disable=unused-variable
                # step is a dummy
                next_mer = next_mer.next_monomer
                last_mer = last_mer.previous_monomer
                average_ca += next_mer.ca.vector + last_mer.ca.vector
                cnt += 2
        except AttributeError:
            # AttributeError is raised by mers at the beginning and at the end of chain
            # they have no next/previous monomers
            pass

        self.pseudoatoms['backbone_average'] = Pseudoatom(numpy_vec=(average_ca / cnt))

    @property
    def angles(self):
        """Property that returns torsion angles (in order: psi and phi) of residue."""
        try:
            return self.dynamic_properties['angles']
        except KeyError:
            self.calculate_angles()
            return self.dynamic_properties['angles']

    def calculate_angles(self):
        """Calculates torsion angles of residue and fills 'angles' property."""
        ang_psi, ang_phi = 0., 0.

        try:
            pd_resid = self.structure.prody_structure['', self.my_chain, self.get_pdb_id()[1]]
            try:
                ang_psi = prody.calcPsi(pd_resid, radian=True)
            except ValueError:
                pass

            try:
                ang_phi = prody.calcPhi(pd_resid, radian=True)
            except ValueError:
                pass
        except AttributeError:
            prm = self.previous_monomer
            nxm = self.next_monomer

            atoms = [self.atoms['N'], self.atoms['CA'], self.atoms['C']]

            pl2 = pydesc.geometry.Plane.build(*atoms)

            if prm is not None:
                pl3 = pydesc.geometry.Plane.build(*([prm.atoms['C']] + atoms[:2]))
                ang_phi = pl2.dihedral_angle(pl3)

            if nxm is not None:
                pl1 = pydesc.geometry.Plane.build(*(atoms[1:]+[nxm.atoms['N']]))
                ang_psi = pl1.dihedral_angle(pl2)



        self.dynamic_properties['angles'] = (ang_psi, ang_phi)

    @property
    def ca(self):     # pylint:disable=invalid-name
        # 'ca' is valid for us
        """Property that returns current residue alpha carbon Atom object."""
        #~ return self.atoms['CA ']
        return self.atoms['CA']

    def calculate_cbx(self):
        """Adds Pseudoatom containing coordinates of the point that lies 1A farther from carbon alfa, than does carbon beta; or carbon alfa coordinates for GLY."""
        if self.get_config('old_cbx_calculation'):
            self.calculate_cbx_legacy()
            return
        if self.name == "GLY":
            #~ self.pseudoatoms['cbx'] = Pseudoatom(numpy_vec=self.atoms['CA '].vector, name='cbx')
            self.pseudoatoms['cbx'] = Pseudoatom(numpy_vec=self.atoms['CA'].vector, name='cbx')
        else:
            #~ ca = self.atoms['CA '].vector
            ca = self.atoms['CA'].vector
            #~ cb = self.atoms['CB '].vector
            cb = self.atoms['CB'].vector
            vec = cb - ca
            nrm = norm(vec)
            vec = vec * ((nrm + 1) / nrm)

            cbx = ca + vec

            self.pseudoatoms['cbx'] = Pseudoatom(numpy_vec=cbx, name='cbx')

    def calculate_cbx_legacy(self):
        """Creates pydesc.geometry.Coord instance containing coordinates of cbx calculated in legacy mode and assigns it to residue cbx property.

        Needs numpy to proceed. Uses Kabsch algorithm to superpose patternal set of C, CA, N and CBX, eeven for GLY.
        """
        pattern = [
            [1.26462, -0.673997, -3.024425],
            [0, 0, -2.5],
            [0, 0, 0],
            [-1.23670, -0.656232, -3.010602]]
        # positions of atoms/points C, C alfa, C beta extended by 1 A and N, respectively
        try:
            coords = [self.atoms[i] for i in ("C", "CA", "CB", "N")]
        except KeyError:
            coords = (self.atoms['C'], self.atoms['CA'], (self.atoms['CA'] - self.atoms['C']) + (self.atoms['CA'] - self.atoms['N']), self.atoms['N'])
        bb_coords = [coord_obj.get_coord() for coord_obj in coords]
        #

        # ===============
        # foregin code starts here
        # ===============
        # all subsequent comments untile notice were made by author
        # assertions replaced by ValueErrors
        # pylint:disable=invalid-name, no-member

        # check for consistency
        if len(bb_coords) != len(pattern):
            raise ValueError('Wrong lenght of backbone: mer %s' % str(self))
        L = len(bb_coords)

        # must alway center the two proteins to avoid
        # affine transformations.  Center the two proteins
        # to their selections.
        COM1 = numpy.sum(bb_coords, axis=0) / float(L)
        COM2 = numpy.sum(pattern, axis=0) / float(L)
        bb_coords = bb_coords - COM1
        pattern = pattern - COM2

        # This beautiful step provides the answer. V and Wt are the orthonormal
        # bases that when multiplied by each other give us the rotation matrix, U.
        # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
        V, S, Wt = numpy.linalg.svd(numpy.dot(numpy.transpose(pattern), bb_coords))

        # we alredy have our solution, in the aaults from SVD.
        # we just need to check for reflections and then produce
        # the rotation.  V and Wt are orthonormal, so their det's
        # are +/-1.
        reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        # U is simply V*Wt
        U = numpy.dot(V, Wt)

        # rotate and translate the molecule
        pattern = numpy.dot((pattern), U) + COM1
        pattern = pattern.tolist()

        # pylin: enable=invalid-name, no-member
        # =============
        # end of foregin code
        # =============
        self.pseudoatoms['cbx'] = Pseudoatom(*pattern[2], name='cbx')


class Nucleotide(MonomerChainable):  # TODO: Improve ConfigManager access

    """Representation of a nucleotide."""

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Nucleotide constructor.

        Arguments:
        pdb_residue -- BioPython Bio.PDB.Residue instance based on which Nucleotide is being created.
        structure -- the Structure instance to which Nucleotide belongs.

        Raises Warning if given pdb_residue does not contain proper atoms or if its atoms occur in wrong distances.
        Extended MonomerChainable method.
        See also config file docstring.

        Config parameteres in branch ConfigManager.monomer.nucleotide:
        min_o5'_p_dist
        max_o5'_p_dist
        min_c5'_o5'_dist
        max_c5'_o5'_dist
        min_c4'_c5'_dist
        max_c4'_c5'_dist
        min_c3'_c4'_dist
        max_c3'_c4'_dist
        min_c3'_o3'_dist
        max_c3'_o3'_dist
        """
        MonomerChainable.__init__(self, structure_obj, ind, name, chain, atoms)

        rats = self.get_config('ring_atoms')

        def flag(name, atom):
            if name in rats:
                atom.ring_flag = True
                return True
            atom.ring_flag = False

        self.ring_atoms = {name: atom for name, atom in self.atoms.items() if flag(name, atom)}

        self.calculate_ring_center()
        self.calculate_proximate_ring_center()
        self.calculate_ring_plane()
        self.calculate_nx()
        self.ion_neighbours = []
        self._check_name('nucleotide')

    def calculate_ring_center(self):
        """Adds pseudoatom representing base ring center."""
        #~ vec = (self.ring_atoms['N1 '].vector + self.ring_atoms['C4 '].vector) * 0.5
        vec = (self.ring_atoms['N1'].vector + self.ring_atoms['C4'].vector) * 0.5
        self.pseudoatoms['ring_center'] = Pseudoatom(numpy_vec=vec, name='ring_center')

    def calculate_ring_plane(self):
        """Adds pydesc.geometry.Plane object representing base to current nucleotide pseudoatom dictionary."""
        #~ at1, at2, at3 = self.ring_atoms['C2 '], self.ring_atoms['C4 '], self.ring_atoms['C6 ']
        at1, at2, at3 = self.ring_atoms['C2'], self.ring_atoms['C4'], self.ring_atoms['C6']
        self.ring_plane = pydesc.geometry.Plane.build(at1, at2, at3)     # pylint:disable=attribute-defined-outside-init
        # current method is called by init

    def calculate_proximate_ring_center(self):
        """Adds pseudoatom representing center of the base ring beeing closer to glycosidic bond."""
        try:
            vec = numpy.array([0., 0., 0.])
            for at in ('C4', 'C5', 'N7', 'C8', 'N9'):
                vec += self.atoms[at].vector
            vec /= 5.
            self.pseudoatoms['prc'] = Pseudoatom(numpy_vec=vec, name='prc')
        except KeyError:
            pass

    @property
    def prc(self):
        try:
            return self.pseudoatoms['prc']
        except KeyError:
            return self.pseudoatoms['ring_center']

    def calculate_nx(self):
        """Adds pseudoatom representing extended by 1.4A vector along glycosidic bond."""
        at1 = self.atoms["C1'"]
        try:
            at2 = self.N9
        except AttributeError:
            at2 = self.N1
        vec = (at2 - at1).vector
        nrm = norm(vec)
        nvec = vec * ((nrm + 1.4) / nrm)

        nx = at1.vector + nvec
        self.pseudoatoms['nx'] = Pseudoatom(numpy_vec=nx, name='nx')


class MonomerOther(Monomer):

    """Abstract class, representation for ligands.

    Subclasses:
    Ion
    Ligand
    """

    # __metaclass__ = ABCMetamonomer

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Monomer Other constructor.

        pdb_residue -- instance of BioPython residue.
        structure_obj -- instance of pydesc structure object mer should be attached to.

        Extends superclass method.
        """
        Monomer.__init__(self, structure_obj, ind, name, chain, atoms)

# pylint:disable=no-self-use
# following methods are needed and should not be a function
    def _has_bond(self, monomer):
        """Returns False, as no mer is next for non-chainable monomers"""
        return False

    @property
    def previous_monomer(self):
        """Returns None.

        For monomers other then chainable this property cannot be set to any value other then None.
        """
        return None

    @property
    def next_monomer(self):
        """Returns None.

        For monomers other then chainable this property cannot be set to any value other then None.
        """
        return None
# pylint:enable=no-self-use


class Ion(MonomerOther):

    """Representation of an ion ligand."""

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Ion constructor.

        Sets basic attributes.

        Arguments:
        pdb_residue -- Bio.PDB.Residue instance representing ion.
        structure_obj -- instance of parental PyDesc structure.
        """
        super(Ion, self).__init__(structure_obj, ind, name, chain, atoms)
        if len(self.atoms) != 1:
            raise ValueError("Failed to create Ion, given BioPython residue consists of to many atoms.")

    def get_radius(self):
        name = max(self.atoms)
        try:
            return self.get_config('radii')[name]
        except KeyError:
            warn(NoConfiguration("No radius for %s ions." % name))
            return 2.5


class Ligand(MonomerOther):

    """Representation of any ligand except ions."""

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Ligand constructor.

        Sets basic attributes.

        Arguments:
        pdb_residue -- Bio.PDB.Residue instance representing ligands other than ions.
        structure_obj -- instance of parental PyDesc structure.
        """
        super(Ligand, self).__init__(structure_obj, ind, name, chain, atoms)



Monomer.reset_config_cache()
