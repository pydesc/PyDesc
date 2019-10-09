import numpy
import scipy.linalg

import pydesc.geometry
import pydesc.geometry
from pydesc.warnexcept import IncompleteParticle
from pydesc.warnexcept import NoConfiguration
from pydesc.warnexcept import UnknownParticleName
from pydesc.warnexcept import warn
from pydesc.warnexcept import WrongAtomDistances
from . import ConfigManager

norm = scipy.linalg.get_blas_funcs("nrm2")


class Atom(pydesc.geometry.Coord):
    """Representation of atoms described in pdb files.

    Subclass of pydesc.geometry.Coord class.

    Attributes:
    name -- string containing atom name.
    element -- string, element name.
    pdb_atom -- instance of BioPython Atom class.
    """

    def __init__(
            self, coords, element, occupancy=0.0, b_factor=0.0
    ):  # pylint:disable=super-init-not-called
        # there is no need to call dict.__init__
        """Atom constructor.

        Arguments:
        coords -- np.array of x, y, z coordinates.
        element -- one letter string indicating element.
        occupancy -- value from pdb file, 0.0 by defautl
        b_factor -- value from pdb file, 0.0 by defautl
        """
        self.vector = coords
        self.element = element
        self.occupancy = occupancy
        self.b_factor = b_factor

    def __repr__(self):
        return "<Atom at %s>" % " ".join(["%.2f" % i for i in self.vector])


class Pseudoatom(pydesc.geometry.Coord):
    """Representation of any point related to monomer other than atom.

    Subclass of pydesc.geometry.Coord class.

    Attributes:
    name -- string containing point name.
    """

    def __repr__(self):
        return "<Pseudoatom %s: %f %f %f>" % ((self.name,) + tuple(self.vector))

    def __init__(
            self, x=0.0, y=0.0, z=0.0, numpy_vec=None, name=""
    ):  # pylint:disable=super-init-not-called
        # there is no need to call dict.__init__
        """Pseudoatom constructor.

        Arguments:
        x, y, z -- pseudoatom coordinates.
        numpy_vec -- NumPy array containing coordinates (if provided XYZ are
        ignored). None by default.
        name -- string, pseudoatom name. None by default.
        owner -- instance od pydesc.monomer.Monomer subclass that contains
        pseudoatom. None by default.
        calc_method -- owner method to calculate pseudoatom coordinates in
        dynamic mode. None by default.

        Sets attribute 'dynamic' to False. If set to True - coordinates are
        calculated each time when readed.
        """
        self.name = name
        pydesc.geometry.Coord.__init__(self, x, y, z, numpy_vec)


class DynamicPropertiesDict(dict):
    """Class of dicts to store values that need to be recalculated for every
    frame of molecular dynamics trajectory."""

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

        Returns value of given key if 'dynamic' attrbute is set to False.
        Otherwise tries to return value of the key if it is not None. If
        value is None - forces owner to recalculate the value byt calling
        "calculate_<key>" method.
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


class Mer:
    """Abstract class, representation of mers and particles present in
    molecular structures.

    Subclasses:
    MonomerChainable -- residues and nucleotides.
    MonomerOther -- ligands or their type: ions.
    """

    @staticmethod
    def is_chainable():
        """Return True if mer is chainable."""
        return False

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
            if issubclass(cls.__base__, Mer):  # pylint: disable=no-member
                res = cls.__base__.get_config(
                    prop_name
                )  # pylint:disable=no-member, protected-access
                # __base__ is not absent
                # protected access to superclass method
            else:
                raise

        return res

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Monomer constructor.

        Arguments:
        structure_obj -- Structure in which current monomer is included.
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

        self.pseudoatoms = DynamicPropertiesDict(self)
        self.dynamic_properties = DynamicPropertiesDict(self)
        self._ss = "="

    def __len__(self):
        """Return sum of lengths of monomer's atoms and pseudoatoms."""
        return len(list(iter(self)))

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
        return iter(
            [self.atoms[atom] for atom in sorted(self.atoms)]
            + [self.pseudoatoms[point] for point in sorted(self.pseudoatoms)]
        )

    def __getattr__(self, name):
        """Returns proper attribute value.

        Argument:
        name -- string, attribute name.
        """
        name = name.lstrip()
        try:
            return object.__getattribute__(self, "atoms")[name]
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
        warn(
            DeprecationWarning(
                """Atom eventually won't inherit from dict type, so avoid 
            getting to attributes via getitem.
       Use getattr instead, e.g.
       instead of
       >>> print my_atom['rc']
       use
       >>> print my_atom.rc
       or access atoms or pseudoatoms dicts directly:
       >>> print my_atom.pseudoatoms['rc']
       """
            ),
            1,
        )
        name = name.lstrip()
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            try:
                return object.__getattribute__(self, "atoms")[name]
            except KeyError:
                try:
                    return self.pseudoatoms[name]
                except (AttributeError, KeyError):
                    repr_ = self.ind if self.ind is not None else str(self)
                    raise AttributeError(
                        "Monomer %s has no attribute %s" % (repr_, name)
                    )

    def finalize(self):
        """Method called by structures to calculate and set attributes that
        need structural information to be calculated.
        """
        self.calculate_rc()

    def has_bond(self, mer):
        return False

    def calculate_rc(self):
        """Sets Monomer's attribute rc (geometrical center).

        Adds pydesc.geometry.Coord instance representing the geometrical
        center of a mer to mers pseudoatoms dict.
        If possible, only sidechain atoms are taken into account.
        """
        non_backbone_coordinates = [a.vector for a in self.iter_nbb_atoms()]

        if non_backbone_coordinates:
            vector = numpy.average(non_backbone_coordinates, 0)
        else:
            try:
                vector = self.ca.vector
            except AttributeError:
                vector = self.atoms["P  "].vector
        self.pseudoatoms["rc"] = Pseudoatom(numpy_vec=vector, name="rc")

    def iter_atoms(self):
        """Returns iterator that iterates over monomer's atoms."""
        return iter(list(self.atoms.values()))

    def iter_bb_atoms(self):
        """Returns iterator that iterates over monomer's backbone atoms."""
        return iter([])

    def iter_nbb_atoms(self):
        """Returns iterator that iterates over monomer's all atoms except
        backbone."""
        return self.iter_atoms()

    @classmethod
    def seq_3to1(cls, seq):
        """Returns a one letter code for a given 3-letter code."""
        try:
            cls_name = cls.__name__.lower()
            code_dictionary = getattr(
                getattr(ConfigManager.mers, cls_name), cls_name + "_code"
            )  # pylint:disable=no-member
            try:
                additional_dictionary = getattr(
                    getattr(ConfigManager.mers, cls_name), cls_name + "_additional_code"
                )  # pylint:disable=no-member
            except AttributeError:
                additional_dictionary = {}
            return (
                code_dictionary[seq]
                if seq in code_dictionary
                else additional_dictionary[seq]
            )
        except AttributeError:
            if issubclass(cls.__base__, Mer):  # pylint:disable=no-member
                # ??? Monomer has no __base__
                return cls.__base__.seq_3to1(seq)  # pylint:disable=no-member
                # ??? same here
            raise AttributeError("No dictionary defined for class %s", str(cls))

    @classmethod
    def seq_1to3(cls, let):
        """Returns a three letter code for a given 1-letter code. In
        ambiguous cases the first
        matching code is returned.
        """
        try:
            cls_name = cls.__name__.lower()
            code_dictionary = getattr(
                getattr(ConfigManager.mers, cls_name), cls_name + "_code"
            )  # pylint:disable=no-member
            for seq3, seq1 in list(code_dictionary.items()):
                if seq1 == let:
                    return seq3
            raise KeyError(
                "Cannot translate %s to 3 letter code" % (cls_name + " symbol " + let,)
            )
        except AttributeError:
            if issubclass(cls.__base__, Mer):  # pylint:disable=no-member
                return cls.__base__.seq_1to3(let)  # pylint:disable=no-member
                # Monomer has __base__ attr
            raise AttributeError("No dictionary defined for class %s", str(cls))

    @property
    def seq(self):
        """Returns one letter code for mer if possible ("?" if name is
        unknown)."""
        try:
            return self.seq_3to1(self.name)
        except KeyError:
            warn(UnknownParticleName(self))
            return "?"
        except AttributeError:
            warn(
                NoConfiguration(
                    "class %s has no dictionary in configuration manager, "
                    "thus '=' inserted into sequence. to turn this exception "
                    "into harmless warning - set NoConfiguration in "
                    "ConfigManager.warnings_and_exceptions.class_filters to "
                    "'ignore' or ;always'" % self.__class__.__name__
                )
            )
            return "="

    @property
    def seq3(self):  # ??? trzyliterowy kod to imie w wypadku aa, a reszta?
        """Returns mer three letter pdb name."""
        return self.name

    @property
    def representation(self):
        """Returns indicators of current monomer set in configuration
        manager."""
        return [getattr(self, indicator) for indicator in self.get_config("indicators")]

    def get_pdb_id(self):
        """Returns pdb id if possible, otherwise returns None."""
        try:
            return self.structure.converter.get_pdb_id(self.ind)
        except AttributeError:
            return None

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
        self._next_monomer = None
        self._previous_monomer = None
        self._check_name()

    def _check_name(self):
        """Method that raises warning if unknown particle name was found in
        pdb file."""
        try:
            self.seq_3to1(self.name)
        except KeyError:
            data = (
                type(self).__name__.capitalize(),
                str(self.get_pdb_id()),
                self.ind or 0,
                str(self.structure),
                self.name,
            )
            warn(
                UnknownParticleName(
                    "%s %s (no. %i) from %s has incorrect name: %s." % data
                )
            )

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

    @property
    def next_mer(self):
        """Property that returns monomer following current mer in its
        structure."""
        try:
            return self._next_monomer
        except AttributeError:
            return None

    @next_mer.setter
    def next_mer(self, value):
        """Property that returns monomer following current mer in its
        structure."""
        self._next_monomer = value

    @property
    def previous_mer(self):
        """Property that returns monomer preceding current mer in its
        structure."""
        try:
            return self._previous_monomer
        except AttributeError:
            return None

    @previous_mer.setter
    def previous_mer(self, value):
        """Property that returns monomer preceding current mer in its
        structure."""
        self._previous_monomer = value

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
        the next monomer or None if distance cannot be computed."""
        try:
            return abs(self.backbone_average - self.next_mer.backbone_average)
        except AttributeError:
            return None


class MerOther(Mer):
    """Abstract class, representation for ligands.

    Subclasses:
    Ion
    Ligand
    """

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Monomer Other constructor.

        pdb_residue -- instance of BioPython residue.
        structure_obj -- instance of pydesc structure object mer should be
        attached to.

        Extends superclass method.
        """
        Mer.__init__(self, structure_obj, ind, name, chain, atoms)

    # pylint:disable=no-self-use
    # following methods are needed and should not be a function
    def _has_bond(self, monomer):
        """Returns False, as no mer is next for non-chainable mers"""
        return False

    @property
    def previous_monomer(self):
        """Returns None.

        For mers other then chainable this property cannot be set to any
        value other then None.
        """
        return None

    @property
    def next_monomer(self):
        """Returns None.

        For mers other then chainable this property cannot be set to any
        value other then None.
        """
        return None
