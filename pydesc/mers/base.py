from functools import wraps

import numpy
import scipy.linalg

import pydesc.geometry
import pydesc.geometry
from pydesc.mers import ConfigManager
from pydesc.warnexcept import IncompleteParticle
from pydesc.warnexcept import UnknownParticleName
from pydesc.warnexcept import warn
from pydesc.warnexcept import WrongAtomDistances

norm = scipy.linalg.get_blas_funcs("nrm2")
NotSet = object()


def register_pseudoatom(method):
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
    """Representation of atoms described in pdb files.

    Subclass of pydesc.geometry.Coord class.

    Attributes:
    name -- string containing atom name.
    element -- string, element name.
    pdb_atom -- instance of BioPython Atom class.
    """

    def __init__(
        self, coords, element, serial_number=None, occupancy=0.0, b_factor=0.0
    ):
        """Initialize Atom.

        Arguments:
        coords -- np.array of x, y, z coordinates.
        element -- one letter string indicating element.
        serial_number -- serial number in pdb file, None by default.
        occupancy -- value from pdb file, 0.0 by default
        b_factor -- value from pdb file, 0.0 by default
        """
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
        """Initialize with topology atom and trajectory to read coords from."""
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
        """Shadow vector attribute with readings from trajectory object."""
        return self.trajectory.get_atom_coords(self.atom)

    def __repr__(self):
        return "<AtomProxy to #%i from %s>" % (
            self.atom.serial_number,
            repr(self.trajectory),
        )


class Pseudoatom(pydesc.geometry.Coord):
    """Representation of any point related to monomer other than atom.

    Subclass of pydesc.geometry.Coord class.

    Attributes:
    name -- string containing point name.
    """

    def __repr__(self):
        name = (" " + self.name) or ""
        coords = "%f %f %f" % tuple(self.vector)
        return "<Pseudoatom%s at %s>" % (name, coords)

    def __init__(self, x=0.0, y=0.0, z=0.0, numpy_vec=None, name=None):
        """Pseudoatom constructor.

        Arguments:
        x, y, z -- pseudoatom coordinates.
        numpy_vec -- NumPy array containing coordinates (if provided XYZ are
        ignored). None by default.
        name -- string, pseudoatom name. None by default.
        """
        self.name = name
        super().__init__(x, y, z, numpy_vec)


class Mer:
    """Abstract class, representation of mers and particles present in
    molecular structures.

    Subclasses:
    MonomerChainable -- residues and nucleotides.
    MonomerOther -- ligands or their type: ions.
    """

    @classmethod
    def reset_config_cache(cls):
        """Resets cache of configuration settings in this class and all
        subclasses. Should be called after relevant changes in
        ConfigManager."""
        cls._config_cache = {}
        for sub in cls.__subclasses__():
            sub.reset_config_cache()

    @classmethod
    def get_config(cls, prop_name):
        """Returns class configuration from configuration manager.

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
        """Returns class configuration from configuration manager.

        Argument:
        prop_name -- name of configuration to be returned.
        """

        try:
            cls_name = cls.__name__.lower()

            branch = ConfigManager.mers  # pylint: disable=no-member
            if cls_name != "mer":
                branch = getattr(branch, cls_name)

            res = getattr(branch, prop_name)
        except AttributeError:
            if issubclass(cls.__base__, Mer):
                res = cls.__base__.get_config(prop_name)
            else:
                raise

        return res

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Initialize Mer.

        Arguments:
        structure_obj -- Structure in which current mer is included.
        name -- str; mers name.
        chain -- str; chain name.
        atoms -- dict; dict of str names of atoms as keys and Atom instances
        as values.

        Sets attributes:
        name -- mer or ligand name, up to three letters, according to PDB file.
        structure -- the Structure instance to which the monomer belongs.
        chain -- character of the chain that the mers belong to, according
        to PDB file.
        atoms - dict of atoms building current monomer represented by Atom
        instances.
        ind -- PyDesc integer.
        pseudoatoms -- dict of Pseudoatoms.
        dynamic_properties -- dict of other geometrical properties like
        planes for cyclic chemical compounds.
        _ss -- secondary structure sign.
        """

        self.structure = structure_obj
        self.name = name
        self.chain = chain
        self.ind = ind
        self.atoms = atoms
        self.pseudoatoms = {}
        self.dynamic_features = {}

        self._ss = "="

    def __len__(self):
        """Return sum of lengths of monomer's atoms."""
        return len(self.atoms)

    def __repr__(self):
        try:
            return "<%s: %s no. %i, PDB: %s>" % (
                self.__class__.__name__,
                self.name,
                self.ind,
                str(self.get_pdb_id()),
            )
        except (TypeError, KeyError):
            return "<%s: %s>" % (self.__class__.__name__, self.name)

    def __iter__(self):
        """Return monomer iterator.

        Monomer iterator iterates over its atoms and pseudoatoms dictionaries.
        """
        return iter(self.atoms.values())

    def __getattr__(self, name):
        """Returns proper attribute value.

        Argument:
        name -- string, attribute name.
        """
        name = name.lstrip()
        try:
            return self.atoms[name]
        except KeyError:
            repr_ = self.ind if self.ind is not None else str(self)
            raise AttributeError("Monomer %s has no attribute %s" % (repr_, name))

    def has_bond(self, mer):
        """Tell if mer has bond with another mer. Returns False by default. Meant to
        be overwritten in different chainable mers."""
        return False

    def get_pdb_id(self):
        """Returns pdb id if possible, otherwise returns None."""
        try:
            return self.structure.converter.get_pdb_id(self.ind)
        except AttributeError:
            return None

    def reset_dynamic_cache(self):
        """Reset pseudoatoms and dynamic features."""
        self.dynamic_features = {}
        self.pseudoatoms = {}

    @register_pseudoatom
    def rc(self):
        """Sets Monomer's attribute rc (geometrical center).

        Adds pydesc.geometry.Coord instance representing the geometrical
        center of a mer to mers pseudoatoms dict.
        If possible, only side chain atoms are taken into account.
        """
        non_backbone_coordinates = [a.vector for a in self.iter_nbb_atoms()]

        if non_backbone_coordinates:
            vector = numpy.average(non_backbone_coordinates, 0)
        else:
            try:
                vector = self.ca.vector
            except AttributeError:
                vector = self.atoms["P  "].vector
        return Pseudoatom(numpy_vec=vector, name="rc")

    @property
    def seq(self):
        """Returns one letter code for mer if possible ("?" if name is
        unknown)."""
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

    @property
    def representation(self):
        """Returns indicators of current monomer set in configuration
        manager."""
        return [getattr(self, indicator) for indicator in self.get_config("indicators")]

    @property
    def pid(self):
        """Return PDB id as string."""
        return str(self.get_pdb_id())

    @property
    def secondary_structure(self):
        """Secondary structure obtained with DSSP for maternal structure.

        If DSSP was not found or secondary structure was not calculated -
        returns '=' sign.

        To calculate secondary structure - use maternal structure method
        set_secondary_structure.

        See Bio.PDB.DSSP documentary for information about code explanation.
        """
        return self._ss

    @property
    def simple_secondary_structure(self):
        """Secondary structure in simple 3-letter code.

        H -- helix
        E -- extended strand
        C -- coil
        """
        temp = ConfigManager.structure_mon.simple_secondary_structure_code
        return temp[self._ss]


class MerChainable(Mer):
    """Abstract class, representation of residue or nucleotide.

    Subclasses:
    Residue
    Nucleotide
    """

    @staticmethod
    def is_chainable():
        """Return True if mer is chainable."""
        return True

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Chainable monomer constructor.

        Extends superclass method.
        """
        Mer.__init__(self, structure_obj, ind, name, chain, atoms)

        try:
            if self.get_config("check_distances"):
                backbone_atoms = dict(
                    (atom_name, None) for atom_name in self.get_config("backbone_atoms")
                )
                for atom_pair in self.get_config("crucial_atom_distances"):
                    self._check_distance(backbone_atoms, *atom_pair)
            self._check_bbatoms()
        except (AttributeError, KeyError):
            data = type(self).__name__, self.get_pdb_id()
            msg = (
                "Backbone atoms lacking, unable to create %s from residue " "%s" % data
            )
            raise IncompleteParticle(msg)

        self._asa = None
        self._next_mer = None
        self._previous_mer = None

    def _check_bbatoms(self):
        tuple(self.iter_bb_atoms())

    def _check_distance(self, atoms, name_1, name_2, min_dist, max_dist):
        """Raises WrongAtomDistances if atoms distance doesn't meet class
        criteria.

        Arguments:
        atoms -- dictionary containig atom names as keys and Atom instance
        as values.
        name_1, name_2 -- 1st and 2nd atom names.

        NOTE: method requires attributes in configuration manager.
        They must be integers and their names should match pattern:
        min_<first lower and stripped name>_<second lower and stripped
        name>_dist
        and
        max_<first lower and stripped name>_<second lower and stripped
        name>_dist
        names should be given in default order for strings.
        """
        distance = norm(atoms[name_1].vector - atoms[name_2].vector)
        if not min_dist <= distance <= max_dist:
            raise WrongAtomDistances(name_1.strip(), name_2.strip(), self)

    def has_bond(self, monomer):
        """Returns True if the Monomer is followed by a given Monomer.

        Argument:
        monomer -- MonomerChainable instance.

        Calculates distance between backbone atoms of Monomers. Returns True or
        False according to the configurable monomer_acceptable_distance.
        """
        if type(monomer) != type(self):
            return False
        bb_atoms = self.get_config("backbone_atoms")
        last_atom = self.atoms[bb_atoms[-1]]
        next_atom = monomer.atoms[bb_atoms[0]]
        try:
            distance = (last_atom - next_atom).calculate_length()
            return distance <= ConfigManager.mers.monomer_acceptable_distance
        except UnboundLocalError:
            return False

    def iter_bb_atoms(self):
        """Returns iterator that iterates over monomer's backbone atoms."""
        bb_atoms = self.get_config("backbone_atoms")
        return iter([self.atoms[attr_name] for attr_name in bb_atoms])

    def iter_nbb_atoms(self):
        """Returns iterator that iterates over monomer's all atoms except
        backbone."""
        bb_atoms = self.get_config("backbone_atoms")
        return iter(
            [
                atom
                for atom_name, atom in list(self.atoms.items())
                if atom_name not in bb_atoms
            ]
        )

    def adjusted_length(self):
        """Returns distance between backbone_average pseudoatoms of this and
        the next mer or None if distance cannot be computed."""
        try:
            return abs(self.backbone_average - self.next_mer.backbone_average)
        except AttributeError:
            return None

    @register_pseudoatom
    def rc(self):
        """Return Pseudoatom storing geometrical center of side chain."""
        non_backbone_coordinates = [a.vector for a in self.iter_nbb_atoms()]
        vector = numpy.average(non_backbone_coordinates, 0)
        return Pseudoatom(numpy_vec=vector, name="rc")


class MerOther(Mer):
    """Abstract class, representation for ligands."""

    @staticmethod
    def is_chainable():
        """Return True if mer is chainable."""
        return False

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Monomer Other constructor.

        pdb_residue -- instance of BioPython residue.
        structure_obj -- instance of pydesc structure object mer should be
        attached to.

        Extends superclass method.
        """
        Mer.__init__(self, structure_obj, ind, name, chain, atoms)

    @register_pseudoatom
    def rc(self):
        """Return Pseudoatom storing geometrical center of mer."""
        coordinates = [a.vector for a in self]

        vector = numpy.average(coordinates, 0)
        return Pseudoatom(numpy_vec=vector, name="rc")
