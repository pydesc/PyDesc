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
PyDesc Selections.

created: 03.09.2013 - 23.09.2013, Tymoteusz 'hert' Oleniecki
"""
import operator
import re
from abc import ABCMeta
from abc import abstractmethod
from functools import reduce

from pydesc.mers.base import Mer
from pydesc.mers.factories import WrongMerType
from pydesc.structure.topology import PartialStructure
from pydesc.structure.topology import Segment


class Selector:
    """Class able to create new or entangled (sub)structures from selections."""

    def __init__(self, mer_factory):
        """Initialize selector with mer factory."""
        self.mer_factory = mer_factory

    def create_new_structure(self, selection, structure_obj):
        """Creates new user structure instance under selection conditions
        (mers are copied).

        Arguments:
        selection -- selection instance.
        structure_obj -- an instance of any abstract structure subclass.
        distinguish_chains -- True or False; temporarily changes property
        distinguish_chains. Initially set to None, if so default value set
        by constructor is used.

        CHANGES IN structure_obj.converter (its number converter) WILL
        INFLUENCE CREATED User Structure OBJECT.
        """
        set_selection = selection.specify(structure_obj)
        structure = PartialStructure(
            [],  # initialize with empty set of mers
            structure_obj.derived_from.converter,
        )
        list_of_inds = set_selection.get_list_of_inds(structure_obj)
        mf = self.mer_factory
        mers = [mf.copy_mer(structure_obj[ind]) for ind in list_of_inds]
        structure.set_mers(mers)
        return structure


class Selection(metaclass=ABCMeta):
    """A selections of mers of a (sub)structures."""

    def __init__(self):
        """Selections constructor.

        Creates a registry of selected structures as a dictionary,
        which keys are (sub)structures and values are instances of user
        structures.
        """
        pass

    def __add__(self, selection):
        return SelectionsUnion([self, selection])

    def __mul__(self, selection):
        return SelectionsIntersection([self, selection])

    def __sub__(self, selection):
        return SelectionsComplement([self, selection])

    def __repr__(self):
        return "<Selection: %s>" % self.__class__.__name__.lower()

    @abstractmethod
    def specify(self, structure_obj):
        """Returns Set of current selection pdb-id tuples for given structure.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new selection.
        distinguish_chains -- True or False; changes the behaviour of a new
        set in such a way that if the value of a chain is set to False,
        then the chain's character is ignored.
        Initially set to None, if so default value set by constructor is used.
        """
        pass

    def create_structure(self, structure_obj):
        """Creates user structure instance under selection conditions.

        Arguments:
        selection -- selection instance.
        structure_obj -- an instance of any abstract structure subclass.

        CHANGES IN structure_obj.converter (its number converter) WILL
        INFLUENCE CREATED User Structure OBJECT.
        """
        set_selection = self.specify(structure_obj)
        substructure = set_selection.create_structure(structure_obj)
        return substructure

    @staticmethod
    def _finalize_specify(list_of_inds, converter):
        """Creates an instance of selection.Set out of given list of pydesc
        indices.

        Arguments:
        list_of_inds -- a list of PyDesc integers of mers that are to be
        changed into pdb-id tuples.
        converter -- an instance of number converter to be used to convert
        PyDesc integers to pdb-id tuples.
        """
        list_of_pdb_ids = list(map(converter.get_pdb_id, list_of_inds))
        return Set(list_of_pdb_ids)


class Set(Selection):
    """Set of mers PyDesc inds."""

    def __init__(self, list_of_pdb_ids):
        """Set selection constructor.

        Extended superclass method.
        Stores PDB-id tuples of meres to be selected.
        To learn about PDB-id tuples see number converter docstring.

        Arguments:
        list_of_pdb_ids -- list of PDB-id tuples.
        distinguish_chains -- True or False; changes the behaviour of a new
        set in such a way that if the value of a chain is set to False,
        then the chain's character is ignored. Initially set to True.
        """
        Selection.__init__(self)
        self.ids = [id_ for id_ in list_of_pdb_ids]

    def __iter__(self):
        """Returns Set iterator."""
        return iter(self.ids)

    def __repr__(self):
        return "<Selection: set of %s mers>" % str(len(self.ids))

    def get_list_of_inds(self, structure_obj):
        converter = structure_obj.derived_from.converter
        list_of_inds = converter.get_list_of_inds(self.ids)
        return list_of_inds

    def specify(self, structure_obj):
        """Returns Set of current set pdb-id tuples.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new selection.
        """
        list_of_inds = self.get_list_of_inds(structure_obj)
        get_pdb_id_method = structure_obj.derived_from.converter.get_pdb_id
        list_of_pdb_ids = list(map(get_pdb_id_method, list_of_inds))
        return Set(list_of_pdb_ids)

    def create_structure(self, structure_obj):
        """Creates user structure instance under selection conditions.

        Arguments:
        selection -- selection instance.
        structure_obj -- an instance of any abstract structure subclass.

        CHANGES IN structure_obj.converter (its number converter) WILL
        INFLUENCE CREATED User Structure OBJECT.
        """
        inds = self.get_list_of_inds(structure_obj)
        converter = structure_obj.derived_from.converter
        mers = [structure_obj[ind] for ind in inds]
        substructure = PartialStructure(mers, converter)
        return substructure


class Range(Selection):
    """Range of PyDesc inds to select."""

    def __init__(self, pdb_start, pdb_end):
        """Range selection constructor.

        Extended superclass method.
        Stores PDB-id tuples of first and last monomer from range to be
        selected.
        NOTE: last monomer WILL be selected.
        Selections lengths made on different structures may differ due to
        mers with insertion code. It is not possible to select range if
        first or last monomer is not present in given structure.
        To learn about PDB-id tuples see number converter docstring.

        Arguments:
        pdb_start -- a starting residue/nucleotide PDB-id tuple,
        pdb_end -- a ending residue/nucleotide PDB-id tuple.
        """
        Selection.__init__(self)
        self.start = pdb_start
        self.end = pdb_end

    def __repr__(self):
        return "<Selection: range %s - %s>" % (str(self.start), str(self.end))

    def specify(self, structure_obj):
        """Creates user structure instance made of mers from range (
        including last one).

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        """
        converter = structure_obj.derived_from.converter
        start_ind, end_ind = converter.get_list_of_inds([self.start, self.end])
        ids = converter.ind2pdb[start_ind : end_ind + 1]
        return Set(list_of_pdb_ids=ids)

    def create_segment(self, structure_obj):
        """Returns pydesc.structure.Segment."""
        return Segment(structure_obj[self.start], structure_obj[self.end])


class MerAttr(Selection, metaclass=ABCMeta):
    """Class of selections based on different mer attributes."""

    def __init__(self, value):
        super().__init__()
        self.value = value

    @property
    @abstractmethod
    def attr_name(self):
        """Name of attribute determining if mers are to be selected or not."""
        pass

    def specify(self, structure_obj):
        """Return Set selection containing ids specific for given structure."""
        list_of_inds = [
            mer.ind
            for mer in structure_obj
            if getattr(mer, self.attr_name) == self.value
        ]
        return Selection._finalize_specify(
            list_of_inds, structure_obj.derived_from.converter
        )


class ChainSelection(MerAttr):
    """Selection of all mers signed with given chain character."""

    def __init__(self, chain_name):
        MerAttr.__init__(self, chain_name)

    def __repr__(self):
        return "<Selection: chain %s>" % str(self.value)

    @property
    def attr_name(self):
        return "chain"


class MerName(MerAttr):
    """Selection of all mers with given name."""

    def __init__(self, chain_name):
        MerAttr.__init__(self, chain_name)

    @property
    def attr_name(self):
        return "name"

    def __repr__(self):
        return "<Selection: mers called %s>" % self.value.lstrip()


class MerExactType(MerAttr):
    """Select mers of given type (without subclasses).

    Compare to MerSubclasses.
    """

    def __init__(self, cls):
        if not issubclass(cls, Mer):
            raise WrongMerType("Given class has to be subclass of Mer.")
        MerAttr.__init__(self, cls)

    @property
    def attr_name(self):
        return "__class__"

    def __repr__(self):
        return "<Selection of mer exact type: %s>" % self.value.__name__


class MerSubclasses(Selection):
    """Selection of mers of given type."""

    def __init__(self, monomer_subclass):
        """Monomer type selection constructor.

        Extended superclass method.
        Stores a type of monomer - subclass of monomer class.
        NOTE: to learn about monomer subclasses see monomer module
        documentation.

        Argument:
        monomer_subclasses -- a class of mers to be selected.
        """
        Selection.__init__(self)
        self.monomer_subclass = monomer_subclass

    def __repr__(self):
        pattern = re.compile("[A-Z][a-z]*")
        operation = re.findall(pattern, self.monomer_subclass.__name__)
        operation.reverse()
        operation = " ".join(operation).lower()
        return "<Selection: all %s>" % (operation + "s")

    def specify(self, structure_obj):
        """Creates user structure instance made of mers of given type.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        """
        list_of_inds = [
            monomer.ind
            for monomer in structure_obj
            if isinstance(monomer, self.monomer_subclass)
        ]
        return Selection._finalize_specify(
            list_of_inds, structure_obj.derived_from.converter
        )


class Everything(Selection):
    """Selection of all possible mers created by PyDesc."""

    def __init__(self):
        """Overall selection constructor (extended superclass method)."""
        Selection.__init__(self)

    def specify(self, structure_obj):
        """Returns given structure as instance of user structure.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- always set to True.
        """
        list_of_inds = [i.ind for i in structure_obj]
        return self._finalize_specify(
            list_of_inds, structure_obj.derived_from.converter
        )

    @staticmethod
    def create_structure(structure_obj):
        """Create structure in optimal way."""
        return structure_obj


class Nothing(Selection):
    """Empty selection."""

    def __init__(self):
        """Empty selection constructor (extended superclass method)."""
        Selection.__init__(self)

    def specify(self, structure_obj):
        """Returns empty (containing 0 mers) user structure instance.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        """
        return Set([])


class CombinedSelection(Selection, metaclass=ABCMeta):
    """Abstract class, a selection obtained via operation on other
    selections."""

    def __init__(self, list_of_selections):
        """Complex selections constructor.

        Extended method.
        Adds an attribute 'selections', which is a list of selections to be
        used in operations.

        Argument:
        distinguish_chains -- True, False or None; initially set to None.
        If so, all sub-selections use their own distinguish chains setup.
        True/False forces all sub-selections to act like they have property
        'distinguish_chains' set to True/False. See selection docstring to
        learn more.
        """
        Selection.__init__(self)
        self.selections = list_of_selections

    def __iter__(self):
        """Returns iterator over list of sub-selections."""
        return iter(self.selections)

    def iter_recursively(self):
        """Returns recursive iterator over sub-selections."""

        def iter_combined(sel):
            """Returns recursive iterator if possible, otherwise - returns
            given obj.
            """
            try:
                return list(sel.iter_recursively())
            except AttributeError:
                return [sel]

        return iter(reduce(operator.add, list(map(iter_combined, self.selections))))

    def _create_list_of_id_sets(self, structure_obj):
        list_of_id_sets = []
        for selection in self.selections:
            ids = selection.specify(structure_obj).ids
            list_of_id_sets.append(set(ids))
        return list_of_id_sets

    def __repr__(self):
        pattern = re.compile("[A-Z][a-z]*")
        name = self.__class__.__name__.replace("Selections", "")
        operation = " ".join(re.findall(pattern, name)).lower()
        operation = operation.capitalize()
        return "<%s of %i selections>" % (operation, len(self.selections))

    @abstractmethod
    def specify(self, structure_obj):
        """Extended superclass method."""
        pass


class SelectionsUnion(CombinedSelection):
    """An union of two selections.

    E.g.:
    union_1_u_2 = selection_1 + selection_2
    """

    def specify(self, structure_obj):
        """Creates user structure instance made of mers that meet any
        criteria given by sub-selections.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        """
        inds_sets = self._create_list_of_id_sets(structure_obj)
        list_of_ids = set.union(*inds_sets)
        return Set(list_of_ids)


class SelectionsIntersection(CombinedSelection):
    """An intersection of two selections.

    E.g.:
    intersection_1_x_2 = selection_1 * selection_2
    """

    def specify(self, structure_obj):
        """Creates user structure instance made of mers that meet all
        criteria given by sub-selections.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        """
        try:
            list_of_id_sets = self._create_list_of_id_sets(structure_obj)
            list_of_ids = list(set.intersection(*list_of_id_sets))
        except (TypeError, KeyError):
            raise ValueError("Not enough selections to intersect")
        return Set(list_of_ids)


class SelectionsComplement(CombinedSelection):
    """An complement of two selections.

    E.g.:
    complement_1___2 = selection_1 - selection_2
    """

    def specify(self, structure_obj):
        """Creates user structure instance made of mers meeting criteria of
        first selection, but not the second.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        """
        ids_sets = self._create_list_of_id_sets(structure_obj)
        list_of_ids = ids_sets[0] - ids_sets[1]
        return Set(list_of_ids)
