# pylint: disable=too-many-lines
"""
Classes that represents molecular structures and their substructures.

created: 10.07.2013 - , Tymoteusz 'hert' Oleniecki
"""

import operator
import math
import os.path

from abc import ABCMeta
from abc import abstractmethod
from Bio.PDB import DSSP
from io import StringIO
from functools import reduce

import pydesc.dbhandler
import pydesc.geometry
import pydesc.mers
import pydesc.numberconverter
import pydesc.selection
import pydesc.contacts
from pydesc.config import ConfigManager
from pydesc.warnexcept import warn
from pydesc.warnexcept import set_filters
from pydesc.warnexcept import DiscontinuityError
from pydesc.warnexcept import WrongMerType
from pydesc.warnexcept import Info

try:
    import prody
except ImportError:
    warn(Info("No module: prody"))

# pylint: disable=no-member
ConfigManager.new_branch("element")
ConfigManager.element.set_default("element_chainable_length", 5)
ConfigManager.new_branch("structure")
ConfigManager.structure.set_default("dssp_path", "dssp")


# pylint: enable=no-member


class StructureLoader(object):
    """Loads structures from the databases using a given designation."""

    def __init__(self,
                 handler=pydesc.dbhandler.MetaHandler(),
                 parser=pydesc.dbhandler.MetaParser(QUIET=True),
                 mer_factory=pydesc.mers.MerFactory(),
                 ):
        """Structure loader constructor.

        Argument:
        handler -- an instance of handler.
        """  # TODO fix docstring
        self.handler = handler
        self.parser = parser
        self.mer_factory = mer_factory

    def load_bundle(self, code, paths=None, mapping=None):
        """Returns list of structures assembled from a few files.

        Arguments:
        code -- string, database designation of the structure. If paths are given, code becomes loaded structure name.
        path -- list of strings, paths to files to be loaded.
        mapping -- dict as described in PDBBundleHander.get_mapping.

        Path and mapping are optional. If they are not given handler is used to access files from cache.

        Method designed to deal with PDB bundles, but can be used to load any structure given in a few PDB files.
        """
        paths = self.handler.get_file(code) if paths is None else [open(path) for path in paths]
        mapping = self.handler.get_mapping(code) if mapping is None else mapping

        def map_chain_char(mdl, dct, pth):
            """Returns model with re-named chains."""
            fln = pth.split('/')[-1]
            for ch in mdl:
                ch.id = dct[fln][ch.id]
            return mdl

        class TempMdl(object):

            """Temporary model used to store data while loading few PDB files."""

            def __init__(self, *args):
                """TempMdl init. Stores all args in args attribute."""
                self.args = args

            def __iter__(self):
                """TempMdl iterator."""
                return iter(reduce(operator.add, map(list, map(operator.methodcaller('get_chains'), self.args))))

            def get_full_id(self):
                """Returns PDB mer full id."""
                return self.args[0].get_full_id()

        stcll = [[] for i in paths]
        chsll = []
        tmp = []
        for i, path in enumerate(paths):
            with path as fhdlr:
                pdb_structure = self.parser.get_structure(code, fhdlr)
            tmp.append(pdb_structure)
            for j, pdb_model in enumerate(pdb_structure):
                try:
                    chsll[j].append(map_chain_char(pdb_model, mapping, path.name))
                except IndexError:
                    chsll.append([map_chain_char(pdb_model, mapping, path.name)])
                stcll[i].append(list(pdb_model.get_residues()))
        mrlst = [reduce(operator.add, mdls_tup) for mdls_tup in zip(*stcll)]
        # list of lists of mers. each sublist contains mers from all pdb files from bundle
        # (structure divided into separate pdb files) for one model.
        numcon = pydesc.numberconverter.NumberConverter(pdb_mers_list=mrlst)
        return [Structure(TempMdl(*chs_tpl), numcon, [i.name for i in paths]) for chs_tpl in chsll]

    def _get_files_and_path(self, code, path):
        """Return path and list of open handlers to pdb files to be read.

        Arguments:
            code -- code passed by user (str or None if local file was passed).
            path -- path to local file.
        """
        open_files = self.handler.get_file(code) if path is None else [open(path)]
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

        def pick_mer(dct, most_frequent, others):
            """Try to pick *most_frequent* key from given dict *dct*, otherwise pick first from list of *others*."""
            try:
                return dct[most_frequent]
            except KeyError:
                for key in others:
                    try:
                        return dct[key]
                    except KeyError:
                        continue
            raise ValueError('Got empty dict.')

        mers = []
        hits = dict((klass, 0) for klass in self.mer_factory.chainable)

        for pdb_residue in pdb_chain:
            mer_dct, warns = self.mer_factory.create_from_biopdb(
                pdb_residue=pdb_residue,
                structure_obj=structure,
                warn_in_place=False
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
        """Returns a list of Structure instances and the NumberConverter instance.

        Arguments:
        code -- string, database designation of the structure.
        path -- string; path to file to be opened.

        To choose specific database in case of MetaHandler code should be given in following format:
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
        structures = [self._create_structure(pdb_model, path, converter) for pdb_model in models]

        return structures


class AbstractStructure(object):
    """Abstract class, representation of all the structures and their derivatives.

    NOTE:
    PICKING SLICES OF STRUCTURES RETURNS LIST OF MERS INCLUDING LAST INDICATED MER

    Subclasses:
    Structure -- molecular strucutre of a protein or a nucleic acid.
    Segment -- any continuous structure of a structure.
    Element -- central mer with two following and two preceding mers.
    Contact -- two Elements in contact.
    Descriptor -- all Contacts containing central Element.
    """

    __metaclass__ = ABCMeta

    def __init__(self, derived_from):
        """(Sub)structure constructor.

        Argument:
        derived_form -- structure, which self is derived from. Structures loaded from files and user structures are derived from themselvs.
        """
        self.derived_from = derived_from
        self._mers = ()
        if self == derived_from:
            self.trt_matrix = pydesc.geometry.TRTMatrix()
        else:
            self.trt_matrix = self.derived_from.trt_matrix
        self._hash_monomers = None

    def __add__(self, structure_obj):
        """Returns UserStructure or Segment containing all mers present in current and given structure.

        Argument:
        structure_obj -- instance of AbstractStructure subclass.

        If given mers contained in two added structures are subsequent mers - Segment is returned.
        """
        mers = sorted(set(list(self) + list(structure_obj)), key=operator.attrgetter('ind'))
        try:
            return Segment(mers=mers)
        except (ValueError, DiscontinuityError):
            # ValueError is raised by Segment.__init__
            return PartialStructure(mers, self.derived_from.converter)

    def __contains__(self, monomer_obj):
        """Checks if given mer is present in current structure (mer needs to have ind attr)."""
        return monomer_obj in self._mers

    def __getitem__(self, key):
        """Returns mer or list of mers.

        Argument:
        key -- slice or int. If slice, returns list of mer INCLUDING SECOND INDICATED MER with specific PyDesc integers (inds), otherwise returns specified mer.
        """

        # pylint: disable=no-member
        def get_hash_if_possible(param):
            """Returns _mers index for given slice parameter.

            Slice parameter can be monomer ind (PyDesc integer), negative index,
            string convertable to PDB_id or PDB_id itself.
            For monomer inds - monomer's index on _mers list is returned.
            Negative values and 0 are not changed.
            Other values raise IndexError.
            Strings are converted to PDB_id.
            For PDB_id - coresponding ind is taken from numberconverter and again
            """
            if isinstance(param, str):
                # strings are converted to PDB_id
                param = pydesc.numberconverter.PDBid.create_from_string(param)
            if isinstance(param, pydesc.numberconverter.PDBid) or isinstance(param, tuple):
                # if given parameter already is a PDB_id or a coresponding tuple instance
                param = self.derived_from.converter.get_ind(param)
            if isinstance(param, pydesc.mers.Mer):
                return self._mers.index(param)
            try:
                # if parameter is an integer - it is probably monomer ind
                return self._hash(param)
            except KeyError:
                # if not - only negative indexes and 0 are supported
                if not param < 1:
                    raise IndexError(
                        "%s out of mers list range" % str(param))
                return param

        if isinstance(key, slice):
            if key.step is not None:
                raise TypeError("AbstractStructure does not support steps")
            if key.start is None:
                key = slice(self._mers[0].ind, key.stop, None)
            if key.stop is None:
                key = slice(key.start, self._mers[-1].ind, None)
            start, end = map(get_hash_if_possible, [key.start, key.stop])
            try:
                segment = Segment(self._mers[start], self._mers[end])
                if not all(mer in self for mer in segment):
                    raise DiscontinuityError()
                return segment
            except (DiscontinuityError, ValueError):
                end = end + 1 if end != -1 else None
                return PartialStructure(self._mers[start: end], self.derived_from.converter)
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
        """Returns index on _mers list corespoding to given PyDesc integer (ind).

        Argument:
        ind -- PyDesc integer.
        """
        if self._hash_monomers is None:
            self._set_hash()
        return self._hash_monomers[ind]

    def _set_hash(self):
        """Sets _has_monomers attribute.

        _has_monomers is dictionary containing all mers inds as keys and their indexes on _mers list as values.
        It is used by __getitem__ as hashlist.
        """
        self._hash_monomers = dict((monomer_obj.ind, index) for index, monomer_obj in enumerate(self._mers))

    def create_pdb_string(self, transformed=True):
        """Returns an StringIO pdb-like object.

        Argument:
        transformed -- initially set to True, if so - creates PyMOL object with respect for all previous movements; otherwise uses coordinates from pdb file.
        """
        line_n = 0
        components = []
        for monomer_obj in self:
            for atom in monomer_obj.iter_atoms():
                pdb_id = monomer_obj.get_pdb_id()
                coord = atom.get_coord(self.trt_matrix) if transformed else atom.get_coord()
                insertion_code = pdb_id.icode if pdb_id.icode is not None else ' '
                values = (line_n,
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
                          atom.element)
                components.append(values)
                line_n += 1
        pdb_line = "ATOM  %5i %4s %3s%2s%4i%1s%11.3f%8.3f% 8.3f%6.2f %5.2f      %3s%2s"
        components = [pdb_line % v for v in sorted(components, key=lambda vals: vals[0])]
        components.append("END")
        return StringIO("\n".join(components))

    def next_mer(self, monomer_obj):
        """Returns next monomer avalible in current structure for given monomer.

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

        Affects structure trt_matrix. Transformed coordinates of points are calculated when get_transformed_coord method is called on coord instance.
        Argument:
        rotation_matrix -- list of three lists of three floats.
        """
        self.trt_matrix.add_rotation(rotation_matrix)

    def translate(self, vector):
        """Translates all points related to structure.

        Affects structure trt_matrix. Transformed coordinates of points are calculated when get_transformed_coord method is called on coord instance.
        Argument:
        vector -- list of three floats.
        """
        self.trt_matrix.add_translation(vector)

    def select(self):
        """Returns set selection of all related mers.

        Overridden in frequently used types of structures for creation of most intuituv selection type.
        """
        return pydesc.selection.Set(map(self.derived_from.converter.get_pdb_id, map(operator.attrgetter('ind'), self)))

    def adjusted_number(self):
        """
        Returns a putative number of 'straight' segments.

        In case of protein structures segments can contains hairpins and other motifs with sharp bends.
        In some cases it is useful to know the number of 'straight' segments in such a structure, assuming
        that it fits a tight space (e.g. a sphere). This trick is used in CompDesc to compare protein descriptors.

        This implementation first creates a UserStructure instance.
        """

        return self[:].adjusted_number()

    def _map_mers_with_attr(self, attr, skip_other=True):
        """Returns a string of the given mer attribute.

        Arguments:
        attr -- string; name of the attribute that stores string in mers.
        skip_other -- True or False; by default set on True. If so - only chainable mers are considered.
        """
        if skip_other:
            objs = [i for i in self if isinstance(i, pydesc.mers.MerChainable)]
        else:
            objs = list(self)
        sequence = map(operator.attrgetter(attr), objs)
        return "".join(sequence)

    def get_sequence(self):
        """Returns (sub)structure sequence (one letter code)."""
        return self._map_mers_with_attr('seq')

    def get_chain(self, char):
        """Returns chain of given name or None if no such chain is available."""
        for chn in self.chains:
            if chn.chain_char == char:
                return chn
        raise AttributeError('No chain %s in structure %s.' % (char, str(self)))

    def save_pdb(self, path):
        """Writes (sub)structure into pdb file.

        Arguments:
        path -- string; path to new file.
        """
        with open(path, "w") as file_:
            file_.write(self.create_pdb_string().read())


class Structure(AbstractStructure):
    """Representation of molecular structure of the protein or the nucleotide acids."""

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
        return '<Structure %s>' % self.name

    def __str__(self):
        return self.name

    def select(self):
        """Overridden select method.

        Returns union of selections: ranges for chains and set of nonechainable mers.
        """
        return pydesc.selection.Everything()

    def link_dcd_file(self, path):
        """Reads dcd file and replaces atoms coords with coords from dcd trajectory file.

        Argument:
        path -- string; path to dcd file.

        Since trajectory is linked, all atom coords are taken from current frame.
        Pseudoatoms are recalculated every time their coords are returned.
        If Contact map is attached to structure, it is recalculated befor returning contact values,
        but only if frame has been changed since previous reading.

        Sets trajectory attribute as prody.trajectory.dcdfile.DCDFile object. To get number of frames
        call len function on that object.

        To switch between frames use next_frame method (faster) or simply define structure frame attribute:
        >>> structure.frame = 5
        >>> structure.frame
        5
        """
        self.trajectory = prody.DCDFile(path)
        with open(self.path, "r") as temp_f:
            pdstr = prody.parsePDBStream(temp_f)
        self.prody_structure = pdstr
        self.trajectory.setCoords(pdstr)
        self.trajectory.link(pdstr)
        for mer in self:
            for atom in mer:
                atom.prody_atom = pdstr[mer.get_pdb_id()][atom.name.strip()]

        def set_dyn(mer):
            """Sets attribute 'dynamic' of given mer's pseudoatoms dict to True."""
            mer.pseudoatoms.dynamic = True
            mer.dynamic_properties.dynamic = True

        map(set_dyn, self)

    def disconnect_trajectory(self):
        """Removes trajectory attached to structure and turns off dynamic calculation of psedoatoms coordinates."""
        del self.trajectory
        del self.prody_structure

        def reset_mer(mer):
            """Recalculates dynamic properties of mer."""
            mer.pseudoatoms.dynamic = False
            mer.dynamic_properties.dynamic = False
            for atom in mer.iter_atoms():
                atom.prody_atom = None
            mer.pseudoatoms.recalculate_all_values()
            mer.dynamic_properties.recalculate_all_values()

        map(reset_mer, self)

    @property
    def frame(self):
        """Property that stores current frame of loaded trajectory."""
        try:
            return self.trajectory.nextIndex()
        except AttributeError:
            return None

    @frame.setter
    def frame(self, value):
        """Property that stores current frame of loaded trajectory."""
        ind = self.trajectory.nextIndex()
        if value == ind:
            return
        elif value == 0:
            trj = self.trajectory._filename
            self.disconnect_trajectory()
            self.link_dcd_file(trj)
        elif value > len(self.trajectory):  # or value == 0:
            raise IndexError("No such frame in dcd file.")
        elif value < ind:
            self.trajectory.reset()
        else:
            value = value - ind
        for dummy in range(value):
            self.trajectory.next()
        self.refresh()

    def next_frame(self):
        """Switches to next frame if structure is linked to dcd file."""
        n_fr = self.trajectory.next()
        if n_fr is None:
            raise IndexError("No such frame in dcd file.")
        self.refresh()

    def get_number_of_frames(self):
        """Returns number of frames in linked dcd trajectory file."""
        return len(self.trajectory)

    def refresh(self):
        """Forces recalculation of all dynamic properties of mers (usually after frame change in trajectory)."""

        def reset_dynamic_prop(mer):
            """Resets all dynamic properties."""
            mer.dynamic_properties.reset_all_values()
            mer.pseudoatoms.reset_all_values()

        map(reset_dynamic_prop, self)

    def set_secondary_structure(self, file_path=None, dssp=None):
        """Calculates secondary structure using DSSP.

        Arguments:
        file_path -- handler or path to pdb file. By default value of structures 'path' attribute.
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
        for mer in pydesc.selection.MonomerType(pydesc.mers.MerChainable).create_structure(self):
            pdbid = mer.get_pdb_id()
            try:
                restup = sec_stc[mer.my_chain, (' ', pdbid[1], ' ' if pdbid[2] is None else pdbid[2])]
            except KeyError:
                continue
            mer._ss = restup[2]
            mer._asa = restup[3]

    def get_secondary_structure(self):
        """Returns (sub)structure sequence of secondary structure (dssp code)."""
        return self._map_mers_with_attr('secondary_structure')

    def get_simple_secondary_structure(self):
        """Returns (sub)structure sequence of secondary structure (3-letter code)."""
        return self._map_mers_with_attr('simple_secondary_structure')


class PartialStructure(AbstractStructure):
    """Representation of substructures generated by users."""

    def __init__(self, list_of_monomers, number_converter=None, name=None):
        """User's structure constructor."""
        if not name:
            name = 'PyDescObj'
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
            NotImplementedError("PartialStructure cannot be prepared from mers coming from different mers (yet)")

    def __repr__(self):
        return "<PartialStructure: %s>" % self.name

    def _fill_mers_attrs(self):
        """Sets mers attributes normally set by init.

        Sets next/previous_mer attributes.
        """
        for pair in zip(self._mers[:-1], self._mers[1:]):
            if pair[0]._has_bond(pair[1]):  # pylint: disable=protected-access
                # method is protected from users, not from cooperating class init
                pair[0].next_mer = pair[1]
                pair[1].previous_mer = pair[0]

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
        """Set _mers attribute to tuple of mers in given sequence and finalize structure."""
        self._mers = tuple(sequence_of_mers)
        self.finalize()

    def finalize(self):
        """Finalize after setting mers."""
        self._fill_mers_attrs()
        self._set_segments()

    def adjusted_number(self):
        """
        Returns a putative number of 'straight' segments.

        In case of protein structures segments can contains hairpins and other motifs with sharp bends.
        In some cases it is useful to know the number of 'straight' segments in such a structure, assuming
        that it fits a tight space (e.g. a sphere). This trick is used in CompDesc to compare protein descriptors.

        This implementation sums over segments comprising the structure.
        """

        return sum(seg.adjusted_number() for seg in self.segments)


class Segment(AbstractStructure):
    """Representation of continuous substructures of DNA, RNA or protein structure."""

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
            try:
                while mers[-1].next_mer != end:
                    current_mer = current_mer.next_mer
                    mers.append(current_mer)
            except AttributeError:
                raise DiscontinuityError()
        else:
            mers = sorted(mers, key=operator.attrgetter('ind'))
            start = mers[0]
        AbstractStructure.__init__(self, start.structure)
        self._mers = tuple(mers)
        self._check_continuity()
        if len(self._mers) == 0:
            raise ValueError("Failed to create segment, wrong mers given.")

    def _check_continuity(self):
        """Raises DiscontinuityError if segment is not continous."""
        for (monomer1, monomer2) in zip(self._mers, self._mers[1:]):
            if monomer1.next_mer == monomer2:
                continue
            raise DiscontinuityError(monomer1, monomer2)

    def __repr__(self):
        return '<Segment %s-%s>' % (
            self.derived_from.converter.get_pdb_id(self.start.ind),
            self.derived_from.converter.get_pdb_id(self.end.ind))

    @property
    def start(self):
        """Returns first segment monomer."""
        return self._mers[0]

    @property
    def end(self):
        """Returns last segment monomer."""
        return self._mers[-1]

    def select(self):
        """Overridden select method.

        Returns range selection of all related mers.
        """
        return pydesc.selection.Range(*map(self.derived_from.converter.get_pdb_id, [self.start.ind, self.end.ind]))

    def adjusted_number(self):
        """
        Returns a putative number of 'straight' segments.

        In case of protein structures segments can contains hairpins and other motifs with sharp bends.
        In some cases it is useful to know the number of 'straight' segments in such a structure, assuming
        that it fits a tight space (e.g. a sphere). This trick is used in CompDesc to compare protein descriptors.
        """

        try:
            length = sum(m.adjusted_length() for m in self._mers[2:-2])
        except (AttributeError, TypeError):
            # TO AttributeError seems unrised while adjusted_length returns None is something is wrong
            # instead sum() raises TypeError since cannot add number to None
            return 1

        if length == 0:
            return 1

        try:
            return int(math.ceil(length / self._mers[0].get_config('adjusted_segment_length')))
        except AttributeError:
            return 1


class Chain(AbstractStructure):
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
        for mer in self._mers:
            mer.finalize()

    def __repr__(self, mode=0):
        return '<Chain %s>' % self.name

    @property
    def name(self):
        return self.derived_from.name + self.chain_name

    def select(self):
        """Returns current chain selection."""
        return pydesc.selection.ChainSelection(self.chain_name)


class Element(AbstractStructure):
    """Abstract class, representation of substructures from the Descriptor.

    Subclasses:
    ElementChainable -- continuous five-mer structure.
    ElementOther -- Element settled by Ion or Ligand.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, monomer_obj):
        """Element constructor.

        Argument:
        monomer_obj -- instance of appropriate pydesc.monomer.Monomer subclass.
        """
        AbstractStructure.__init__(self, monomer_obj.structure)
        self.central_monomer = monomer_obj
        self._mers = (self.central_monomer)

    def __repr__(self, mode=0):
        return '<%s of %s>' % (
            str(self.__class__.__name__), str(self.derived_from.converter.get_pdb_id(self.central_monomer.ind)))

    @staticmethod
    def build(monomer_obj):
        """Static builder.

        Returns appropriate Element subclass instance.

        Argument:
        monomer_obj -- instance of pydesc.monomer.Monomer.
        """
        try:
            return ElementChainable(monomer_obj)
        except WrongMerType:
            try:
                return ElementOther(monomer_obj)
            except TypeError:
                raise TypeError("Cannot create element using mer %s" % str(monomer_obj.get_pdb_id()))


class ElementChainable(Element, Segment):
    """Representation of a five-mer Segment.

    It consists of five Residues or five Nucleotides: a central mer, two preceding and two following mers.
    """

    def __init__(self, monomer_obj):  # pylint: disable=super-init-not-called
        """ElementChainable constructor.

        Argument:
        monomer_obj -- MonomerChainable instance that became the element settler.

        Sets ElementChainable's list of Monomers.
        """
        if not isinstance(monomer_obj, pydesc.mers.MerChainable):
            raise WrongMerType("Cannot create chainable element using given mer: %i." % monomer_obj.ind)
        Element.__init__(self, monomer_obj)
        length = ConfigManager.element.element_chainable_length  # pylint: disable=no-member
        if not length % 2 == 1:
            raise ValueError("Wrong element chainable length. Length should be odd.")
        for dummy_step in range(length / 2):
            # instead of Segment init
            start = self._mers[0]
            end = self._mers[-1]
            self._mers = (start.previous_mer,) + self._mers + (end.next_mer,)
            if self._mers.count(None) != 0:
                raise ValueError("Cannot create chianable element for mer %i." % monomer_obj.ind)


class ElementOther(Element):
    """Class corresponding to the ElementChainable, but consisting of a single Ion or Ligand instance."""

    def __init__(self, monomer_obj):
        """Element constructor.

        Argument:
        monomer_obj -- instance of any pydesc.monomer.MonomerOther subclass.
        """
        if not isinstance(monomer_obj, pydesc.mers.MerOther):
            raise TypeError("Wrong monomer given to create ElementOther instance: %s" % str(monomer_obj.get_pdb_id()))
        super(ElementOther, self).__init__(monomer_obj)


class Contact(AbstractStructure):
    """Representation of two close-Element instances."""

    def __init__(self, element1, element2):
        """Contact constructor.

        Arguments:
        element1, element2 -- pydesc.structure.Element instances.

        Sets Contacts's list of Monomers.
        Extended AbstractStructure method.
        """

        self.elements = sorted([element1, element2], key=lambda x: x.central_monomer.ind)
        if element1.derived_from is not element2.derived_from:
            raise ValueError("Impossible to create contact instance with elements derived from different structures")
        if element1.central_monomer.ind == element2.central_monomer.ind:
            raise ValueError("Impossible to create contact using one element")
        self._mers = list(operator.add(*self.elements))
        AbstractStructure.__init__(self, element1.derived_from)

    def __sub__(self, val):
        """Deprecated."""
        warn("""Subtracting contacts is no longer supported. Please, use get_other_element instead.""",
             DeprecationWarning)
        return self.get_other_element(val)

    def __repr__(self):
        inpt = tuple([i.central_monomer.ind for i in self.elements])
        return '<Contact of %s and %s elements>' % inpt

    def get_other_element(self, element_obj):
        """Returns other than given element.

        Argument:
        element_obj -- instance of Element class.
        """
        if element_obj.central_monomer not in [i.central_monomer for i in self.elements]:
            raise ValueError("Given element is not included in contact instance, cannot get other element")
        if element_obj.central_monomer == self.elements[0].central_monomer:
            return self.elements[1]
        else:
            return self.elements[0]

    def select(self):
        """Returns SelectionsUnion containing both contact elements selections."""
        return pydesc.selection.SelectionsUnion(map(Segment.select, self.elements))

    @property
    def value(self, cmap=None):
        """
        Contact value in contact_map associated with a structure contact is derived from.

        This property is required by contacts.DescriptorCriterion.
        """
        if not cmap:
            cmap = self.derived_from.contact_map
        return cmap.get_contact_value(*[i.central_monomer.ind for i in self.elements])


class AbstractDescriptor(AbstractStructure):
    """Representation of PyDesc spatial unit.

    It consists of the central Element and all Elements it is in contact with, according to a given ContactMap ContactCriterion.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, element_obj, list_of_contact_objs):
        """Descriptor contructor.

        Arguments:
        element_obj -- Element instance that becomes central_element.
        list_of_contact_obj -- list of Contact instances that are present in descriptor to be builded.
        """
        AbstractStructure.__init__(self, element_obj.derived_from)
        self.central_element = element_obj
        self.contacts = list_of_contact_objs
        self._mers = [monomer_obj for monomer_obj in list(element_obj)]
        for contact_obj in self.contacts:
            self._mers += [monomer_obj for element in contact_obj.elements for monomer_obj in element if
                           monomer_obj not in self._mers]
            # self._mers += [monomer_obj for element in contact_obj.elements for monomer_obj in element if monomer_obj not in self._mers]
        # self._mers = sorted(set(self._mers), key = lambda monomer_obj: monomer_obj.ind)
        # ??? albo if w liscie merow, albo set
        self._mers = sorted(self._mers, key=lambda monomer_obj: monomer_obj.ind)
        self.segments = []
        self._set_segments()
        self.elements = [self.central_element]
        self.elements_values = {}
        self._set_elements()

    @property
    def cpid(self):
        """Returns PDB id of central element central monomer as string."""
        return str(self.central_element.central_monomer.pid)

    def _set_segments(self):
        """Fills initial attr segments."""

        def add_segment(start, end):
            """Adds ElementChainible to segments list, or Segment if start and end are distant more then proper segment lenght."""
            element_length = ConfigManager.element.element_chainable_length  # pylint: disable=no-member
            segment = Segment(start, end)
            if len(segment) == element_length:
                segment = ElementChainable(segment._monomers[element_length / 2])  # pylint: disable=protected-access
            self.segments.append(segment)

        def check_segments(presegment_1, presegment_2):
            """Reduces two tuple to one, if they contain subsequent mers."""
            if presegment_1[1] == presegment_2[0]:
                return (presegment_1[0], presegment_2[1])
            add_segment(*presegment_1)
            return presegment_2

        switch = lambda mer: (mer, mer.next_mer) if mer.next_mer in self._mers else False
        last_segment = reduce(check_segments, filter(bool, map(switch, self._mers)))
        add_segment(*last_segment)

    def _set_elements(self):
        """Fills initial attrs elements and elements_values."""
        neighbours = {self.central_element.central_monomer: []}
        elements = {}
        for con in self.contacts:
            element1, element2 = con.elements
            for el1, el2 in zip((element1, element2), (element2, element1)):
                elements[el1.central_monomer] = el1
                if con.value == 2:
                    try:
                        neighbours[el1.central_monomer].append(el2.central_monomer)
                    except KeyError:
                        neighbours[el1.central_monomer] = [el2.central_monomer]
                else:
                    if el1.central_monomer in neighbours:
                        pass
                    else:
                        neighbours[el1.central_monomer] = []
        self.elements = sorted(elements.values(), key=lambda elm: elm.central_monomer.ind)

        connected = [self.central_element.central_monomer]

        def dfs(closest):
            """DFS."""
            for i in closest:
                if i in connected:
                    continue
                connected.append(i)
                dfs(neighbours[i])

        dfs(neighbours[self.central_element.central_monomer])

        for element in self.elements:
            self.elements_values[element.central_monomer.ind] = pydesc.PropertiesRecord(
                {'optional': False if element.central_monomer in connected else True})

    def __repr__(self):
        return '<%s of %s:%s>' % (self.__class__.__name__, str(self.derived_from), self.cpid)

    @staticmethod
    def build(element_obj, contact_map=None):
        """Static builder.

        Arguments:
        element_obj -- ElementChainable instance.
        contact_map -- instance of contact map. Set on None by default. If so, contact map of structure from which given element was derived from is taken.

        Returns appropriate descriptor.
        """
        try:
            if isinstance(element_obj.central_monomer, pydesc.mers.Nucleotide):
                return NucleotideDescriptor.build(element_obj, contact_map)
            return ProteinDescriptor.build(element_obj, contact_map)
        except AttributeError:
            msg = "Given element (settled by mer %s) has no descriptor. Perhaps used structure has no contact map."
            raise ValueError(msg % str(element_obj.central_monomer.get_pdb_id()))

    @staticmethod
    def create_descriptors(structure_obj, contact_map=None):
        """Static method that creates all possible Descriptor instances for a given (sub)structure.

        Arguments:
        structure_obj -- any AbstractStructure subclass instance
        contact_map -- ContactMap instance that will be submitted to Descrpitor __init__ method, initially set to None.
        """
        descs = []
        for monomer_obj in structure_obj:
            try:
                descs.append(AbstractDescriptor.build(Element.build(monomer_obj)))
            except (TypeError, ValueError, AttributeError):
                descs.append(None)
        return descs

    def create_coord_vector(self):  # pylint: disable=too-many-locals
        # temorary method
        """Creates files containing information for C++ part of library.

        Number of returned lines per monomer depends on its type and indicators values, that can be changes in configuration manager.
        Returned files:
        coo -- each line contains cooridnates of subsequent points.
        loc -- each line contains monomer type, ind (PyDesc integer) and startline in coo file.
        con -- each line contains contact type, 1st monomer ind, 2nd monomer ind and contact_value under given contact map criterions.
        seg -- each line contains inds of segments' starting and ending mers.

        NOTE: First number in a file is always a number of its lines and num of central monomer
        """
        monomer_types = sorted([mon_type for mon_class in pydesc.mers.Mer.__subclasses__() for mon_type in
                                mon_class.__subclasses__()], key=lambda x: x.__name__)  # pylint: disable=no-member
        # subclasses are avalible
        locations = [[0, monomer_obj.ind, 0] for monomer_obj in self._mers]
        i_num = 1
        for item in locations:
            item[2] = i_num
            for i, mon_type in enumerate(monomer_types):
                if isinstance(self.derived_from[item[1]], mon_type):
                    item[0] = i
                    break
            i_num += len(self.derived_from[item[1]].indicators)
        coords = [coord.get_coord() for monomer_obj in self._mers for coord in monomer_obj if
                  coord.name.lower() in monomer_obj.indicators]
        loc = open(
            "./desc" + str(self.central_element.central_monomer.ind) + ".loc", "w")
        coo = open(
            "./desc" + str(self.central_element.central_monomer.ind) + ".coo", "w")
        con = open(
            "./desc" + str(self.central_element.central_monomer.ind) + ".con", "w")
        seg = open(
            "./desc" + str(self.central_element.central_monomer.ind) + ".seg", "w")
        # ============
        # writing coords here:
        # ============
        loc.write(str(len(locations) + 1) + "\n")
        for i in locations:
            loc.write(str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + "\n")
        loc.close()
        coo.write(str(len(coords) + 1) + "\n")
        for point in coords:
            for variable in point:
                coo.write("%.3f" % variable + " ")
            coo.write("\n")
        coo.close()
        con.write(str(len(self.contacts) + 1) + "\n")
        self.contacts = sorted(self.contacts,
                               key=lambda contact: contact.elements[1].central_monomer.ind if contact.elements[
                                                                                                  0] == self.central_element else
                               contact.elements[0].central_monomer.ind)
        for contact in self.contacts:
            # when it grows up it will be contact classification
            con.write("1 ")
            con.write("%3i" % self.central_element.central_monomer.ind + " ")
            con.write("%3i" % [mer for mer in contact if mer != self.central_element.central_monomer][0].ind + " ")
            con.write(str(contact.derived_from.contact_map.get_contact_value(
                [x.central_monomer.ind for x in contact.elements])) + "\n")
        con.close()
        seg.write(str(len(self.segments) + 1) + "\n")
        self.segments = sorted(self.segments, key=lambda segment: segment.start.ind)
        for i in self.segments:
            seg.write("%i %i" % (i.start.ind, i.end.ind) + "\n")
        seg.close()

    def select(self):
        """Overridden select method.

        Returns union of range selections (for self segments) and set of all remaining mers.
        """
        components = map(Segment.select, self.segments)
        other_mers = [
            mer for mer in self if mer not in reduce(operator.add, self.segments)]
        components.append(
            pydesc.selection.Set(map(operator.methodcaller('get_pdb_id'), other_mers)))
        return pydesc.selection.SelectionsUnion(components)


class ProteinDescriptor(AbstractDescriptor):
    """AbstractDescriptor subclass.

    Representation of descriptors settled by protein mers.
    """

    def __init__(self, element_obj, list_of_contact_objs):
        """Descriptor contructor.

        Arguments:
        element_obj -- Element instance that becomes central_element.
        list_of_contact_obj -- list of Contact instances that are present in descriptor to be builded.
        """
        AbstractDescriptor.__init__(self, element_obj, list_of_contact_objs)

    @staticmethod
    def build(element_obj, contact_map=None):
        """Static ProteinDescriptor builder.

        Arguments:
        element_obj -- ElementChainable instance.
        contact_map -- instance of contact map. Set on None by default. If so, contact map of structure from which given element was derived from is taken.
        """
        if contact_map is None:
            contact_map = element_obj.derived_from.contact_map
        stc = element_obj.derived_from
        contacts_ = sorted(contact_map.get_monomer_contacts(element_obj.central_monomer))

        def create_contact(input_):
            """Returns Contact or None if failed to create one.

            Argument:
            input_ -- tuple containing: ind_1, ind_2 -- mers inds, value -- contact value.
            """
            (ind_1, ind_2, value) = input_
            try:
                return Contact(Element.build(stc[ind_1]), Element.build(stc[ind_2]))
            except (AttributeError, TypeError, ValueError):
                # TypeError is raised during creation of Contacts based on pairs of
                # mers that cannot be used to create Elements
                return

        contacts_objs = filter(bool, map(create_contact, contacts_))
        return ProteinDescriptor(element_obj, contacts_objs)


class NucleotideDescriptor(AbstractDescriptor):
    """Class that represents nucleotide descriptors."""

    def __init__(self, element_obj, list_of_contact_objs):
        """Descriptor contructor.

        Arguments:
        element_obj -- Element instance that becomes central_element.
        list_of_contact_obj -- list of Contact instances that are present in descriptor to be builded.
        """
        AbstractDescriptor.__init__(self, element_obj, list_of_contact_objs)

    @staticmethod
    def build(element_obj, contact_map=None):
        """Static NucleotideDescriptor builder."""

        def bld_con(der, *args):
            """Builds a contact."""
            return pydesc.structure.Contact(Element.build(der[args[0]]), Element.build(der[args[1]]))

        derf = element_obj.derived_from
        cmer = element_obj.central_monomer

        if isinstance(element_obj.central_monomer, pydesc.mers.Ion):
            raise TypeError

        if not contact_map:
            contact_map = element_obj.derived_from.contact_map

        contacts_objs = []
        contacts_inds = [cmer.ind]
        ions = []
        for i in contact_map.get_monomer_contacts(cmer.ind):
            if type(derf[i[1]]) != pydesc.mers.Ion:
                try:
                    contacts_objs.append(bld_con(derf, *i))
                    contacts_inds.append(i[1])
                except ValueError:
                    continue
            else:
                ions.append(derf[i[1]])
        ion_contacts = {}
        for ion in ions:
            cdist = (cmer.ring_center - ion.rc).calculate_length()
            for c_ind in contact_map.contacts[ion.ind]:
                if c_ind in contacts_inds:
                    continue
                dist = cdist + (derf[c_ind].ring_center - ion.rc).calculate_length()
                try:
                    cnto = bld_con(derf, ion.ind, c_ind)
                except ValueError:
                    continue
                try:
                    shortest_cnt = min(ion_contacts[c_ind], (dist, cnto, ion.ind), key=lambda x: x[0])
                except KeyError:
                    shortest_cnt = (dist, cnto, ion.ind)
                ion_contacts[c_ind] = shortest_cnt

        qualif_ion_cons = [bld_con(derf, i, cmer.ind) for i in set([i[2] for i in ion_contacts.values()])]

        contacts_objs.extend([i[1] for i in ion_contacts.values()] + qualif_ion_cons)

        return pydesc.structure.NucleotideDescriptor(element_obj, filter(bool, contacts_objs))
