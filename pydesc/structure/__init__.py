# pylint: disable=too-many-lines
"""
Classes that represents molecular structures and their substructures.

created: 10.07.2013 - , Tymoteusz 'hert' Oleniecki
"""

import math
import operator
import os.path
from abc import ABCMeta
from abc import abstractmethod
from io import StringIO

from Bio.PDB import DSSP

import pydesc.dbhandler
import pydesc.geometry
from pydesc.config import ConfigManager
from pydesc.mers.base import Mer
from pydesc.mers.factories import MerFactory
from pydesc.numberconverter import NumberConverter
from pydesc.numberconverter import PDBid
from pydesc.warnexcept import DiscontinuityError
from pydesc.warnexcept import set_filters
from pydesc.warnexcept import warn
from pydesc.warnexcept import WrongElement

# pylint: disable=no-member
ConfigManager.new_branch("element")
ConfigManager.element.set_default("element_chainable_length", 5)
ConfigManager.new_branch("structure")
ConfigManager.structure.set_default("dssp_path", "dssp")


# pylint: enable=no-member


class StructureLoader(object):
    """Loads structures from the databases using a given designation."""

    def __init__(
        self,
        handler=pydesc.dbhandler.MetaHandler(),
        parser=pydesc.dbhandler.MetaParser(QUIET=True),
        mer_factory=MerFactory(),
    ):
        """Structure loader constructor.

        Argument:
        handler -- an instance of handler.
        """  # TODO fix docstring
        self.handler = handler
        self.parser = parser
        self.mer_factory = mer_factory

    def _get_files_and_path(self, code, path):
        """Return path and list of open handlers to pdb files to be read.

        Arguments:
            code -- code passed by user (str or None if local file was passed).
            path -- path to local file.
        """
        if path is None:
            open_files = self.handler.get_file(code)
        else:
            open_files = [open(path)]
        path = open_files[0].name
        return path, open_files

    @staticmethod
    def _get_code(code, path):
        """Return PDB structure code.

        Arguments:
            code -- None is it should be read from path, code otherwise.
            path -- path to a file or None.
        """
        if code is None:
            code, dummy_ext = os.path.splitext(os.path.basename(path))
        if code.find("://") != -1:
            dummy_db, code = code.split("://")
        return code

    def _get_models(self, files, code):
        """Return list of models read from PDB files in given list.

        Argument:
            files -- list of open PDB file handlers.
            code -- code to be passed to parser.
        """
        models = []
        for handler in files:
            with handler:
                pdb_structure = self.parser.get_structure(code, handler)
            models.extend([pdb_model for pdb_model in pdb_structure])
        return models

    def _create_chain(self, pdb_chain, structure):
        """Return Chain instance.

        Arguments:
            pdb_chain -- Bio.PDB chain.
            structure -- Structure object, for which Chain is to be created.
        """

        def pick_mer(dct, most_frequent, others_):
            """Try to pick *most_frequent* key from given dict *dct*,
            otherwise pick first from list of *others_*."""
            try:
                return dct[most_frequent]
            except KeyError:
                for key in others_:
                    try:
                        return dct[key]
                    except KeyError:
                        continue
            raise ValueError("Got empty dict.")

        mers = []
        hits = dict((klass, 0) for klass in self.mer_factory.chainable)

        for pdb_residue in pdb_chain:
            mer_dct, warns = self.mer_factory.create_from_biopdb(
                pdb_residue=pdb_residue, structure_obj=structure, warn_in_place=False
            )
            if mer_dct is None:
                continue

            mers.append((mer_dct, warns))
            for mer_class in hits:
                hits[mer_class] += int(mer_class in mer_dct)

        winner = max(hits, key=lambda hit: hits[hit])
        others = self.mer_factory.other
        chain_mers = []
        for mer_dct, warns in mers:
            accepted_mer = pick_mer(mer_dct, winner, others)
            chain_mers.append(accepted_mer)
            warns.raise_all(type(accepted_mer))

        return Chain(structure, pdb_chain.get_id(), chain_mers)

    def _create_structure(self, model, path, converter):
        """Return structure for given model

        Argument:
            model -- Bio.PDB model.
            path -- path to pdb file storing the structure.
            converter -- id converter for structure to be created.
        """
        name = model.get_full_id()[0]
        structure = Structure(name, path, converter)

        chains = [self._create_chain(pdb_chain, structure) for pdb_chain in model]

        structure.finalize(chains)
        return structure

    def load_structures(self, code=None, path=None):
        """Returns a list of Structure instances and the NumberConverter
        instance.

        Arguments:
        code -- string, database designation of the structure.
        path -- string; path to file to be opened.

        To choose specific database in case of MetaHandler code should be
        given in following format:
        "<database_name>://<structure_code>", e.g.
        "PDB://1no5".
        To choose specific BioUnit type "<biounit_code>/<number>", e.g.
        "Unit://1no5/1"
        in case of MetaHandler or
        "1no5/1"
        in case of UnitHandler.
        """

        set_filters()
        path, open_files = self._get_files_and_path(code, path)
        code = self._get_code(code, path)
        models = self._get_models(open_files, code)
        converter = pydesc.numberconverter.NumberConverter(models)
        structures = [
            self._create_structure(pdb_model, path, converter) for pdb_model in models
        ]

        return structures


class BackbonedMixIn:
    def _fill_mers_attrs(self):
        """Sets mers attributes normally set by init.

        Sets next/previous_mer attributes.
        """
        for pair in zip(self._mers[:-1], self._mers[1:]):
            if pair[0].has_bond(pair[1]):
                pair[0].next_mer = pair[1]
                pair[1].previous_mer = pair[0]


class AbstractStructure(metaclass=ABCMeta):
    """Abstract class, representation of all the structures and their
    derivatives.

    NOTE:
    PICKING SLICES OF STRUCTURES RETURNS LIST OF MERS INCLUDING LAST
    INDICATED MER

    Subclasses:
    Structure -- molecular structure of a protein or a nucleic acid.
    Segment -- any continuous structure of a structure.
    Element -- central mer with two following and two preceding mers.
    Contact -- two Elements in contact.
    Descriptor -- all Contacts containing central Element.
    """

    def __init__(self, derived_from):
        """(Sub)structure constructor.

        Argument:
        derived_form -- structure, which self is derived from. Structures
        loaded from files and user structures are derived from themselves.
        """
        self.derived_from = derived_from
        self._mers = ()
        if self == derived_from:
            self.trt_matrix = pydesc.geometry.TRTMatrix()
        else:
            self.trt_matrix = self.derived_from.trt_matrix
        self._hash_monomers = None

    def __add__(self, structure_obj):
        """Returns UserStructure or Segment containing all mers present in
        current and given structure.

        Argument:
        structure_obj -- instance of AbstractStructure subclass.

        If given mers contained in two added structures are subsequent mers
        - Segment is returned.
        """
        mers = sorted(
            set(list(self) + list(structure_obj)), key=operator.attrgetter("ind")
        )
        try:
            return Segment(mers=mers)
        except (ValueError, DiscontinuityError):
            # ValueError is raised by Segment.__init__
            return PartialStructure(mers, self.derived_from.converter)

    def __contains__(self, monomer_obj):
        """Checks if given mer is present in current structure (mer needs
        to have ind attr).
        """
        return monomer_obj in self._mers

    def __getitem__(self, key):
        """Returns mer or list of mers.

        Argument:
        key -- slice or int. If slice, returns list of mer INCLUDING SECOND
        INDICATED MER with specific PyDesc integers (inds), otherwise
        returns specified mer.
        """

        # pylint: disable=no-member
        def get_hash_if_possible(param):
            """Returns _mers index for given slice parameter.

            Slice parameter can be monomer ind (PyDesc integer), negative
            index,
            string convertable to PDB_id or PDB_id itself.
            For monomer inds - monomer's index on _mers list is returned.
            Negative values and 0 are not changed.
            Other values raise IndexError.
            Strings are converted to PDB_id.
            For PDB_id - corresponding ind is taken from number converter and
            again
            """
            # TODO: brush and prune this tree
            if isinstance(param, str):
                # strings are converted to PDB_id
                param = PDBid.create_from_string(param)
            if isinstance(param, PDBid) or isinstance(param, tuple):
                # if given parameter already is a PDB_id or a corresponding
                # tuple instance
                param = self.derived_from.converter.get_ind(param)
            if isinstance(param, Mer):
                return self._mers.index(param)
            try:
                # if parameter is an integer - it is probably monomer ind
                return self._hash(param)
            except KeyError:
                # if not - only negative indexes and 0 are supported
                if not param < 1:
                    raise IndexError("%s out of mers list range" % str(param))
                return param

        if isinstance(key, slice):
            if key.step is not None:
                raise TypeError("AbstractStructure does not support steps")
            if key.start is None:
                key = slice(self._mers[0].ind, key.stop, None)
            if key.stop is None:
                key = slice(key.start, self._mers[-1].ind, None)
            start, end = list(map(get_hash_if_possible, [key.start, key.stop]))
            try:
                segment = Segment(self._mers[start], self._mers[end])
                if not all(mer in self for mer in segment):
                    raise DiscontinuityError()
                return segment
            except (DiscontinuityError, ValueError):
                end = end + 1 if end != -1 else None
                return PartialStructure(
                    self._mers[start:end], self.derived_from.converter
                )
        else:
            return self._mers[get_hash_if_possible(key)]
        # pylint: enable=no-member
        # _mers is surely set outside init

    def __iter__(self):
        """Returns iterator that iterates over structure mers."""
        return iter(self._mers)

    def __len__(self):
        """Returns number of mers present in current structure."""
        return len(self._mers)

    def _hash(self, ind):
        """Returns index on _mers list corresponding to given PyDesc integer
        (ind).

        Argument:
        ind -- PyDesc integer.
        """
        if self._hash_monomers is None:
            self._set_hash()
        return self._hash_monomers[ind]

    def _set_hash(self):
        """Sets _has_monomers attribute.

        _has_monomers is dictionary containing all mers inds as keys and
        their indexes on _mers list as values.
        It is used by __getitem__ as hash list.
        """
        self._hash_monomers = dict(
            (monomer_obj.ind, index) for index, monomer_obj in enumerate(self._mers)
        )

    def create_pdb_string(self, transformed=True):  # TODO: move to separate class
        """Returns an StringIO pdb-like object.

        Argument:
        transformed -- initially set to True, if so - creates PyMOL object
        with respect for all previous movements; otherwise uses coordinates
        from pdb file.
        """
        line_n = 0
        components = []
        for monomer_obj in self:
            for atom in monomer_obj.iter_atoms():
                pdb_id = monomer_obj.get_pdb_id()
                if transformed:
                    coord = atom.get_coord(self.trt_matrix)
                else:
                    coord = atom.get_coord()
                insertion_code = pdb_id.icode or " "
                values = (
                    line_n,
                    atom.name,
                    monomer_obj.name,
                    monomer_obj.my_chain,
                    pdb_id[1],
                    insertion_code,
                    coord[0],
                    coord[1],
                    coord[2],
                    atom.occupancy,
                    atom.b_factor,
                    monomer_obj.pdb_residue.get_segid(),
                    atom.element,
                )
                components.append(values)
                line_n += 1
        pdb_line = (
            "ATOM  %5i %4s %3s%2s%4i%1s%11.3f" "%8.3f% 8.3f%6.2f %5.2f      %3s%2s"
        )
        sorted_components = sorted(components, key=lambda vals: vals[0])
        components = [pdb_line % v for v in sorted_components]
        components.append("END")
        return StringIO("\n".join(components))

    def next_mer(self, monomer_obj):
        """Returns next monomer available in current structure for given
        monomer.

        Argument:
        monomer_obj -- instance of pydesc.monomer.Monomer.
        """
        try:
            if monomer_obj.next_mer in self:
                return monomer_obj.next_mer
            return None
        except AttributeError:
            return None

    def rotate(self, rotation_matrix):
        """Rotates all points related to structure.

        Affects structure trt_matrix. Transformed coordinates of points are
        calculated when get_transformed_coord method is called on coord
        instance.
        Argument:
        rotation_matrix -- list of three lists of three floats.
        """
        self.trt_matrix.add_rotation(rotation_matrix)

    def translate(self, vector):
        """Translates all points related to structure.

        Affects structure trt_matrix. Transformed coordinates of points are
        calculated when get_transformed_coord method is called on coord
        instance.
        Argument:
        vector -- list of three floats.
        """
        self.trt_matrix.add_translation(vector)

    def adjusted_number(self):
        """
        Returns a putative number of 'straight' segments.

        In case of protein structures segments can contains hairpins and
        other motifs with sharp bends.
        In some cases it is useful to know the number of 'straight' segments in
        such a structure, assuming
        that it fits a tight space (e.g. a sphere). This trick is used in
        CompDesc to compare protein descriptors.

        This implementation first creates a UserStructure instance.
        """

        return self[:].adjusted_number()

    def _map_mers_with_attr(self, attr, skip_other=True):
        """Returns a string of the given mer attribute.

        Arguments:
        attr -- string; name of the attribute that stores string in mers.
        skip_other -- True or False; by default set on True. If so - only
        chainable mers are considered.
        """
        if skip_other:
            objs = [i for i in self if i.is_chainable()]
        else:
            objs = list(self)
        sequence = list(map(operator.attrgetter(attr), objs))
        return "".join(sequence)

    def get_sequence(self):
        """Returns (sub)structure sequence (one letter code)."""
        return self._map_mers_with_attr("seq")

    def get_chain(self, name):  # TODO: move to structure!
        """Returns chain of given name if it is available, otherwise raises
        AttributeError.
        """
        for chn in self.chains:
            if chn.chain_name == name:
                return chn
        raise AttributeError("No chain %s in %s." % (name, str(self)))

    def save_pdb(self, path):  # TODO: move to separate class
        """Writes (sub)structure into pdb file.

        Arguments:
        path -- string; path to new file.
        """
        with open(path, "w") as file_:
            file_.write(self.create_pdb_string().read())


class Structure(AbstractStructure):
    """Representation of molecular structure of the protein or the
    nucleotide acids.
    """

    def __init__(self, name, path, converter_obj):
        """Structure constructor.

        Arguments:

        Sets Structure's list of mers and list of chains.
        Sets mers' next_mer attribute and creates their elements.
        Extended AbstractStructure method.
        """  # TODO fix docstring

        AbstractStructure.__init__(self, self)
        self.path = path
        self.name = name
        self.converter = converter_obj
        self.chains = None

    def finalize(self, chains):
        self.chains = chains
        self._mers = tuple([mer for chain in chains for mer in chain])
        self._set_hash()

    def __repr__(self):
        return "<Structure %s>" % self.name

    def __str__(self):
        return self.name

    def set_secondary_structure(self, file_path=None, dssp=None):
        """Calculates secondary structure using DSSP.

        Arguments:
        file_path -- handler or path to pdb file. By default value of
        structures 'path' attribute.
        dssp -- optional; command to call DSSP ('dssp' by default).

        Method uses Bio.PDB.DSSP. See docstring for more information.
        """
        if dssp is None:
            dssp = ConfigManager.structure.dssp_path
        if file_path is None:
            file_path = self.path
        elif not isinstance(file_path, str):
            file_path = file_path.name
        sec_stc = DSSP(self.pdb_model, file_path, dssp)
        chainable = [mer for mer in self if mer.is_chainable()]
        for mer in chainable:
            pdbid = mer.get_pdb_id()
            try:
                restup = sec_stc[mer.chain, (" ", pdbid[1], pdbid[2] or " ")]
            except KeyError:
                continue
            mer._ss = restup[2]
            mer._asa = restup[3]

    def get_secondary_structure(self):
        """Returns (sub)structure sequence of secondary structure (dssp
        code).
        """
        return self._map_mers_with_attr("secondary_structure")

    def get_simple_secondary_structure(self):
        """Returns (sub)structure sequence of secondary structure (3-letter
        code).
        """
        return self._map_mers_with_attr("simple_secondary_structure")


class PartialStructure(BackbonedMixIn, AbstractStructure):
    """Representation of substructures generated by users."""

    def __init__(self, list_of_monomers, number_converter=None, name=None):
        """User's structure constructor."""
        if not name:
            name = "PyDescObj"
        self.name = name
        self.converter = None
        self.segments = None
        derived_from = set(monomer_obj.structure for monomer_obj in list_of_monomers)
        if number_converter:
            self.converter = number_converter
        n_derived = len(derived_from)
        if n_derived == 1:
            derived_from = max(derived_from)
            AbstractStructure.__init__(self, derived_from)
            if self.converter is None:
                self.converter = derived_from.converter
            self.set_mers(sorted(list_of_monomers, key=lambda mer: mer.ind))
        elif n_derived == 0:
            AbstractStructure.__init__(self, self)
        else:
            NotImplementedError(
                "PartialStructure cannot be prepared from mers coming from "
                "different mers (yet)"
            )

    def __repr__(self):
        return "<PartialStructure: %s>" % self.name

    def _set_segments(self):
        """Sets segments attribute."""
        self.segments = []
        start = self._mers[0]
        for pair in zip(self._mers, self._mers[1:]):
            set_start = True
            try:
                if self.next_mer(pair[0]) == pair[1]:
                    set_start = False
                    continue
                end = pair[0]
                self.segments.append(Segment(start, end))
            except (AttributeError, DiscontinuityError):
                pass
            finally:
                if set_start:
                    start = pair[1]
        try:
            self.segments.append(Segment(start, self._mers[-1]))
        except DiscontinuityError:
            pass

    def set_mers(self, sequence_of_mers):
        """Set _mers attribute to tuple of mers in given sequence and
        finalize structure."""
        self._mers = tuple(sequence_of_mers)
        self.finalize()

    def finalize(self):
        """Finalize after setting mers."""
        self._fill_mers_attrs()
        self._set_segments()

    def adjusted_number(self):
        """
        Returns a putative number of 'straight' segments.

        In case of protein structures segments can contains hairpins and
        other motifs with sharp bends.
        In some cases it is useful to know the number of 'straight' segments in
        such a structure, assuming that it fits a tight space (e.g. a
        sphere). This trick is used in CompDesc to compare protein descriptors.

        This implementation sums over segments comprising the structure.
        """
        return sum(seg.adjusted_number() for seg in self.segments)


class Segment(AbstractStructure):
    """Representation of continuous substructures of DNA, RNA or protein
    structure.
    """

    def __init__(self, start=None, end=None, mers=None):
        """Segment constructor.

        Arguments:
        start -- starting MonomerChainable instance.
        end -- closing MonomerChainable instance.
        mers --

        Sets the Segment's list of Monomers and checks for continuity.
        """
        if mers is None:
            mers = [start]
            current_mer = start
            while current_mer != end:
                try:
                    current_mer = current_mer.next_mer
                except AttributeError:
                    names = start.get_pdb_id(), end.get_pdb_id()
                    msg = "It impossible to reach %s starting from %s." % names
                    raise DiscontinuityError(msg)
                mers.append(current_mer)
        else:
            mers = sorted(mers, key=operator.attrgetter("ind"))
            start = mers[0]
        AbstractStructure.__init__(self, start.structure)
        self._mers = tuple(mers)
        self._check_continuity()
        if len(self._mers) == 0:
            raise ValueError("Failed to create segment, wrong mers given.")

    def _check_continuity(self):
        """Raises DiscontinuityError if segment is not continuous."""
        for (monomer1, monomer2) in zip(self._mers, self._mers[1:]):
            if monomer1.next_mer == monomer2:
                continue
            raise DiscontinuityError(monomer1, monomer2)

    def __repr__(self):
        return "<Segment %s-%s>" % (
            self.derived_from.converter.get_pdb_id(self.start.ind),
            self.derived_from.converter.get_pdb_id(self.end.ind),
        )

    @property
    def start(self):
        """Returns first segment monomer."""
        return self._mers[0]

    @property
    def end(self):
        """Returns last segment monomer."""
        return self._mers[-1]

    def adjusted_number(self):
        """
        Returns a putative number of 'straight' segments.

        In case of protein structures segments can contains hairpins and
        other motifs with sharp bends.
        In some cases it is useful to know the number of 'straight' segments in
        such a structure, assuming that it fits a tight space (e.g. a
        sphere). This trick is used in CompDesc to compare protein descriptors.
        """

        try:
            length = sum(m.adjusted_length() for m in self._mers[2:-2])
        except (AttributeError, TypeError):
            # TO AttributeError seems to be never raised while adjusted_length
            # returns None is something is wrong
            # instead sum() raises TypeError since cannot add number to None
            return 1

        if length == 0:
            return 1

        try:
            return int(
                math.ceil(length / self._mers[0].get_config("adjusted_segment_length"))
            )
        except AttributeError:
            return 1


class Chain(BackbonedMixIn, AbstractStructure):
    """Representation of a polymer chain.

    As in PDB file, chains contain both: chainable mers and ligands.
    """

    def __init__(self, structure_obj, chain_name, mers):
        """Chain constructor.

        Arguments:
        structure_obj -- pydesc.structure from which chain is derived.
        chain_name -- name of the chain.
        mers -- sequence of mers chain consists of.

        Sets Chain's list of Monomers.
        Extended Segment method.
        """  # TODO fix docstring
        AbstractStructure.__init__(self, structure_obj)
        self.chain_name = chain_name
        self._mers = tuple(mers)
        self._fill_mers_attrs()

    def __repr__(self, mode=0):
        return "<Chain %s>" % self.name

    @property
    def name(self):
        return self.derived_from.name + self.chain_name


class AbstractElement(AbstractStructure, metaclass=ABCMeta):
    """Abstract class, representation of substructures from the Descriptor.

    Subclasses:
    ElementChainable -- continuous five-mer structure.
    ElementOther -- Element settled by Ion or Ligand.
    """

    @abstractmethod
    def __init__(self, mer):
        """Element constructor.

        Argument:
        mer -- instance of appropriate pydesc.monomer.Monomer subclass.
        """
        AbstractStructure.__init__(self, mer.structure)
        self.central_monomer = mer
        self._mers = (mer,)

    def __repr__(self, mode=0):
        return "<%s of %s>" % (
            str(self.__class__.__name__),
            str(self.derived_from.converter.get_pdb_id(self.central_monomer.ind)),
        )


class ElementChainable(AbstractElement, Segment):
    """Representation of a five-mer Segment.

    It consists of five Residues or five Nucleotides: a central mer,
    two preceding and two following mers.
    """

    def __init__(self, mer):
        """ElementChainable constructor.

        Argument:
        mer -- MonomerChainable instance that became the element settler.

        Sets ElementChainable's list of Monomers.
        """
        AbstractElement.__init__(self, mer)
        length = ConfigManager.element.element_chainable_length
        if not length % 2 == 1:
            raise ValueError("Length of chainable element should be odd.")
        msg = "Cannot create chainable element for mer %i." % mer.ind
        for dummy_step in range(length // 2):
            start = self._mers[0]
            end = self._mers[-1]
            try:
                self._mers = (start.previous_mer,) + self._mers + (end.next_mer,)
            except AttributeError:
                raise ValueError(msg)
        if self._mers.count(None) != 0:
            raise ValueError(msg)


class ElementOther(AbstractElement):
    """Class corresponding to the ElementChainable, but consisting of a
    single Ion or Ligand instance.
    """

    def __init__(self, mer):
        """Element constructor.

        Argument:
        mer -- instance of any pydesc.monomer.MonomerOther subclass.
        """
        super(ElementOther, self).__init__(mer)


class Contact(AbstractStructure):
    """Representation of two close-Element instances."""

    def __init__(self, element1, element2):
        """Contact constructor.

        Arguments:
        element1, element2 -- pydesc.structure.Element instances.

        Sets Contacts's list of Monomers.
        Extended AbstractStructure method.
        """
        self.elements = {element1, element2}
        if element1.derived_from is not element2.derived_from:
            raise ValueError(
                "Impossible to create contact instance with elements derived "
                "from different structures"
            )
        if element1.central_monomer.ind == element2.central_monomer.ind:
            raise ValueError("Impossible to create contact using one element")
        AbstractStructure.__init__(self, element1.derived_from)
        self._mers = tuple(element1) + tuple(element2)

    def __sub__(self, val):
        """Deprecated."""
        warn(
            """Subtracting contacts is no longer supported. Please, 
            use get_other_element instead.""",
            DeprecationWarning,
        )
        return self.get_other_element(val)

    def __repr__(self):
        items = (i.central_monomer.ind for i in self.elements)
        return "<Contact of %s and %s elements>" % items

    def get_other_element(self, element_obj):
        """Returns other than given element.

        Argument:
        element_obj -- instance of Element class.
        """
        elements = set(self.elements)
        try:
            elements.remove(element_obj)
        except KeyError:
            raise WrongElement("Given element is not part of this Contact.")
        return elements.pop()

    def value(self, cmap):
        """
        Contact value in contact_map associated with a structure contact is
        derived from.

        This property is required by contacts.DescriptorCriterion.
        """
        return cmap.get_contact_value(*[i.central_monomer.ind for i in self.elements])
