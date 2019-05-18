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
from abc import ABCMeta
from abc import abstractmethod
import operator
import re
from copy import deepcopy
import pydesc.structure
from functools import reduce

from pydesc.mers import MerFactory


class Selection(object):
    """A selections of mers of a (sub)structures."""

    __metaclass__ = ABCMeta

    def __init__(self, distinguish_chains):
        """Selections constructor.

        Creates a registry of selected structures as a dictionary, which keys are (sub)structures and values are
        instances of user structures.

        Argumnet:
        distinguish_chains -- True, False or None; if True chain character is considered while selecting (sub)structure,
        otherwise selection includes all mers with given number and insertion code regardless of chain character.
        None is for complex selections only; in that case all subselections use their own distinguish chains setup.
        """
        self._distinguish_chains = distinguish_chains  # get and set provided by property

    def __add__(self, selection):
        return SelectionsUnion([self, selection])

    def __mul__(self, selection):
        return SelectionsIntersection([self, selection])

    def __sub__(self, selection):
        return SelectionsRelativeComplement([self, selection])

    def __repr__(self):
        return "<Selection: %s>" % self.__class__.__name__.lower()

    def create_structure(self, structure_obj, distinguish_chains=None):
        """Creates user structure instance under selection conditions.

        Arguments:
        structure_obj -- an instance of any abstract structure subclass.
        distinguish_chains -- True or False; temporarily changes property distinguish_chains.
        Initially set to None, if so default value set by constructor is used.

        CHANGES IN structure_obj.converter (its numberconverter) WILL INFLUENCE CREATED User Structure OBJECT.
        """
        if distinguish_chains is None:
            distinguish_chains = self._distinguish_chains
        return self.specify(structure_obj, distinguish_chains).create_structure(structure_obj, distinguish_chains)

    def create_new_structure(self, structure_obj, distinguish_chains=None):
        """Creates new user structure instance under selection conditions (mers are copied).

        Arguments:
        structure_obj -- an instance of any abstract structure subclass.
        distinguish_chains -- True or False; temporarily changes property distinguish_chains.
        Initially set to None, if so default value set by constructor is used.

        CHANGES IN structure_obj.converter (its numberconverter) WILL INFLUENCE CREATED User Structure OBJECT.
        """
        if distinguish_chains is None:
            distinguish_chains = self._distinguish_chains
        return self.specify(structure_obj, distinguish_chains).create_new_structure(structure_obj, distinguish_chains)

    @property
    def distinguish_chains(self):
        """Property that takes boolean values; describes selection ability to distinguish chain character.

        If this property is set to True chain character is considered while selecting (sub)structure, otherwise
        selection includes all mers with given number and insertion code regardless of chain character.
        """
        return self._distinguish_chains

    @abstractmethod
    def specify(self, structure_obj, distinguish_chains=None):
        """Returns Set of current selection pdb-id tuples for given structure.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new selection.
        distinguish_chains -- True or False; changes the behaviour of a new set in such a way that if the value of a
        chain is set to False, then the chain's character is ignored.
        Initially set to None, if so default value set by constructor is used.
        """
        pass

    @staticmethod
    def _finalize_specify(list_of_inds, converter):
        """Creates an instance of selection.Set out of given list of pydesc indices.

        Arguments:
        list_of_inds -- a list of PyDesc integers of mers that are to be changed into pdb-id tuples.
        converter -- an instance of number converter to be used to convert PyDesc integers to pdb-id tuples.
        """
        list_of_pdb_ids = map(converter.get_pdb_id, list_of_inds)
        return Set(list_of_pdb_ids)


# =========== Simple selections
class Set(Selection):
    """Set of mers PyDesc inds."""

    def __init__(self, list_of_pdb_ids, distinguish_chains=True):
        """Set selection constructor.

        Extended superclass method.
        Stores PDB-id tuples of monomeres to be selected.
        To learn about PDB-id tuples see number converter docstring.

        Arguments:
        list_of_pdb_ids -- list of PDB-id tuples.
        distinguish_chains -- True or False; changes the behaviour of a new set in such a way that if the value of a chain is set to False, then the chain's character is ignored. Initially set to True.
        """
        Selection.__init__(self, distinguish_chains)
        self.ids = [id_ for id_ in list_of_pdb_ids]

    def __iter__(self):
        """Returns Set iterator."""
        return iter(self.ids)

    def __repr__(self):
        return "<Selection: set of %s mers>" % str(len(self.ids))

    def create_structure(self, structure_obj, distinguish_chains=None):
        """Creates user structure based on given set of PDB-ids.

        Arguments:
        structure_obj -- an instance of any abstract structure subclass.
        distinguish_chains -- True or False; temporarily changes property distinguish_chains. Initially set to None, if so default value set by constructor is used.

        NOTE: CHANGES IN structure_obj.converter (its number_converter) WILL INFLUENCE CREATED User Structure OBJECT.
        """
        inds = filter(bool, structure_obj.derived_from.converter.get_list_of_inds(self.ids))
        substructure = pydesc.structure.PartialStructure([structure_obj[ind] for ind in inds],
                                                         structure_obj.derived_from.converter)
        return substructure

    def create_new_structure(self, structure_obj, distinguish_chains=None):
        """Creates new user structure based on given set of PDB-ids with copied mers.

        Arguments:
        structure_obj -- an instance of any abstract structure subclass.
        distinguish_chains -- True or False; temporarily changes property distinguish_chains. Initially set to None, if so default value set by constructor is used.

        NOTE: CHANGES IN structure_obj.converter (its number_converter) WILL INFLUENCE CREATED User Structure OBJECT.
        """
        structure = pydesc.structure.PartialStructure([], structure_obj.derived_from.converter)
        inds = filter(bool, structure_obj.derived_from.converter.get_list_of_inds(self.ids))
        mf = MerFactory()
        mers = [mf.copy_mer(structure_obj[ind]) for ind in inds]
        structure.set_mers(mers)
        return structure

    def specify(self, structure_obj, distinguish_chains=None):
        """Returns Set of current set pdb-id tuples.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new selection.
        distinguish_chains -- True or False; changes the behaviour of a new set in such a way that if the value of a chain is set to False, then the chain's character is ignored. Initially set to None, if so default value set by constructor is used.
        """
        if distinguish_chains is None:
            distinguish_chains = self._distinguish_chains
        list_of_inds = filter(bool, structure_obj.derived_from.converter.get_list_of_inds(self.ids, distinguish_chains))
        list_of_pdb_ids = map(structure_obj.derived_from.converter.get_pdb_id, list_of_inds)
        return Set(list_of_pdb_ids, distinguish_chains)


class Range(Selection):
    """Range of PyDesc inds to select."""

    def __init__(self, pdb_start, pdb_end, distinguish_chains=True):
        """Range selection constructor.

        Extended superclass method.
        Stores PDB-id tuples of first and last monomer from range to be selected.
        NOTE: last monomer WILL be selected.
        Selections lenghts made on different structures may differ due to mers with insertion code.
        It is not possible to select range if first or last monomer is not present in given structure.
        To learn about PDB-id tuples see number converter docstring.

        Arguments:
        pdb_start -- a starting residue/nucleotide PDB-id tuple,
        pdb_end -- a ending residue/nucleotide PDB-id tuple.
        distinguish_chains -- True or False; changes the behaviour of a new set in such a way that if the value of a chain is set to False, then the chain's character is ignored. Initially set to True.
        """
        Selection.__init__(self, distinguish_chains)
        # ~ if pdb_start == pdb_end:
        # ~ raise ValueError("Cannot create range selection using one mer (same start and end)")
        self.start = pdb_start
        self.end = pdb_end

    def __repr__(self):
        return "<Selection: range %s - %s>" % (str(self.start), str(self.end))

    def specify(self, structure_obj, distinguish_chains=None):
        """Creates user structure instance made of mers from range.

        Extended superclass method.
        Returns an empty structure if it is impossible to identify first or last monomer in given structure.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- True or False; changes the behaviour of a new set in such a way that if the value of a chain is set to False, then the chain's character is ignored. Initially set to None, if so default value set by constructor is used.
        """
        if distinguish_chains is None:
            distinguish_chains = self._distinguish_chains
        start = pydesc.numberconverter.PDB_id.create_from_string(str(self.start))
        end = pydesc.numberconverter.PDB_id.create_from_string(str(self.end))
        inds = structure_obj.derived_from.converter.get_list_of_inds([start, end], distinguish_chains)
        start_ind, end_ind = inds[::2], inds[1::2]
        list_of_inds = []
        for start, end in zip(start_ind, end_ind):
            if None in (start, end):
                continue
            temp_inds = [monomer.ind for monomer in structure_obj[start:end]]
            list_of_inds.extend(temp_inds)
        return Selection._finalize_specify(self, list_of_inds, structure_obj.derived_from.converter)

    def create_segment(self, structure_obj, distinguish_chains=None):
        """Returns pydesc.structure.Segment."""
        return pydesc.structure.Segment(structure_obj[self.start], structure_obj[self.end])


class ChainSelection(Selection):
    """Selection of all mers signed with given chain character."""

    def __init__(self, chain_char):
        """Chain selection constructor.

        Extended superclass method.
        Stores a character of chain to be selected.

        Argument:
        chain_name -- a character of chain to be selected.
        """
        Selection.__init__(self, True)
        self.chain = chain_char

    def __repr__(self):
        return "<Selection: chain %s>" % str(self.chain)

    def specify(self, structure_obj, distinguish_chains=True):
        """Creates user structure instance made of mers belonging to one chain according to PDB file.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- always set to True.
        """
        chain = [chain for chain in structure_obj.chains if self.chain == chain.chain_char]
        list_of_inds = map(operator.attrgetter('ind'), *chain)
        return Selection._finalize_specify(self, list_of_inds, structure_obj.derived_from.converter)


class MonomerName(Selection):
    """Selection of all mers with given name."""

    def __init__(self, monomer_name):
        """Monomer name selection constructor.

        Extended superclass method.
        Stores a name of monomer.
        NOTE: nucleotide names start with two white characters.

        Argument:
        monomer_name -- a name of nucleotides to be selected.
        """
        Selection.__init__(self, False)
        self.monomer_name = monomer_name

    def __repr__(self):
        return "<Selection: mers called %s>" % self.monomer_name.lstrip()

    def specify(self, structure_obj, distinguish_chains=False):
        """Creates user structure instance made of mers with given name.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- always set to False.
        """
        list_of_inds = [monomer.ind for monomer in structure_obj if monomer.name == self.monomer_name]
        return Selection._finalize_specify(self, list_of_inds, structure_obj.derived_from.converter)


class MonomerType(Selection):
    """Selection of mers of given type."""

    def __init__(self, monomer_subclass):
        """Monomer type selection constructor.

        Extended superclass method.
        Stores a type of monomer - subclass of monomer class.
        NOTE: to learn about monomer subclasses see monomer module documentation.

        Argument:
        monomer_subclasses -- a class of mers to be selected.
        """
        Selection.__init__(self, False)
        self.monomer_subclass = monomer_subclass

    def __repr__(self):
        pattern = re.compile('[A-Z]{1}[a-z]*')
        operation = re.findall(pattern, self.monomer_subclass.__name__)
        operation.reverse()
        operation = " ".join(operation).lower()
        return "<Selection: all %s>" % (operation + "s")

    def specify(self, structure_obj, distinguish_chains=False):
        """Creates user structure instance made of mers of given type.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- always set to False.
        """
        list_of_inds = [monomer.ind for monomer in structure_obj if isinstance(monomer, self.monomer_subclass)]
        return Selection._finalize_specify(self, list_of_inds, structure_obj.derived_from.converter)


class Everything(Selection):
    """Selection of all possible mers created by PyDesc."""

    def __init__(self):
        """Overall selection constructor (extended superclass method)."""
        Selection.__init__(self, True)

    def specify(self, structure_obj, distinguish_chains=False):
        """Returns given structure as instance of user structure.

        Extended superclass method.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- always set to True.
        """
        list_of_inds = map(operator.attrgetter('ind'), structure_obj)
        return self._finalize_specify(list_of_inds, structure_obj.derived_from.converter)

    def create_structure(self, structure_obj, distinguish_chains=None):
        return structure_obj


class Nothing(Selection):
    """Empty selection."""

    def __init__(self):
        """Empty selection constructor (extended superclass method)."""
        Selection.__init__(self, True)

    def create_structure(self, structure_obj, distinguish_chains=False):
        """Returns empty UserStructure.

        Arguments:
        structure_obj -- instance od pydesc.structure.AbstractStructure.
        distinguish_chains -- boolean value set always to False. See pydesc.selection.Selection.create_structure for more informations.
        """
        substructure = pydesc.structure.PartialStructure([], structure_obj.derived_from.converter)
        return substructure

    def specify(self, structure_obj, distinguish_chains=True):
        """Returns empty (containing 0 mers) user structure instance.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- always set to True.
        """
        return Set([])


# =============
class CombinedSelection(Selection):
    """Abstract class, a selection obtained via operation on other selections."""

    __metaclass__ = ABCMeta

    def __init__(self, list_of_selections, distinguish_chains=None):
        """Complex selections constructor.

        Extended method.
        Adds an attribute 'selections', which is a list of selections to be used in operations.

        Argumnet:
        distinguish_chains -- True, False or None; initially set to None. If so, all subselections use their own distinguish chains setup. True/False forces all subselctions to act like they have property 'distinguish_chains' set to True/False. See selection docstring to learn more.
        """
        Selection.__init__(self, distinguish_chains)
        self.selections = list_of_selections

    def __iter__(self):
        """Returns iterator over list of subselections."""
        return iter(self.selections)

    def iter_recursively(self):
        """Returns recurive iterator over subselections."""

        def iter_combined(sel):
            """Returns recursiv iterator if possible, otherwise - returns given obj"""
            try:
                return list(sel.iter_recursively())
            except AttributeError:
                return [sel]

        return iter(reduce(operator.add, map(iter_combined, self.selections)))

    def __repr__(self):
        pattern = re.compile('[A-Z]{1}[a-z]*')
        operation = " ".join(re.findall(pattern, self.__class__.__name__.replace("Selections", ""))).lower()
        return "<%s of %s selections>" % (operation.capitalize(), str(len(self.selections)))

    @abstractmethod
    def specify(self, structure_obj, distinguish_chains=None):
        """Extended superclass method."""
        pass
        # PyLint satisfied


# ============= Combined selections
class SelectionsUnion(CombinedSelection):
    """An union of two selections.

    You can easily create unions adding two selections, e.g.:
    union_1_and_2 = selection_1 + selection_2
    """

    def specify(self, structure_obj, distinguish_chains=None):
        """Creates user structure instance made of mers that meet any criteria given by subselections.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- True or False; changes the behaviour of a new set in such a way that if the value of a chain is set to False, then the chain's character is ignored. Initially set to None, if so default values set by selections' constructors are used.
        """
        if distinguish_chains is None:
            distinguish_chains = self._distinguish_chains
        list_of_pdb_ids = reduce(operator.add,
                                 [list(selection.specify(structure_obj, distinguish_chains)) for selection in
                                  self.selections])
        list_of_pdb_ids = dict((list_of_pdb_ids.index(pdb_id), pdb_id) for pdb_id in list_of_pdb_ids)
        # removing repeated mers by assigning them to the same key in dictionary
        list_of_pdb_ids = [list_of_pdb_ids[key] for key in list_of_pdb_ids.keys()]
        return Set(list_of_pdb_ids, distinguish_chains)


class SelectionsIntersection(CombinedSelection):
    """An intersection of two selections.

    You can easily create intersections multiplying two selections, e.g.:
    intersection_1_x_2 = selection_1 * selection_2
    """

    def specify(self, structure_obj, distinguish_chains=None):
        """Creates user structure instance made of mers that meet all criteria given by subselections.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- True or False; changes the behaviour of a new set in such a way that if the value of a chain is set to False, then the chain's character is ignored. Initially set to None, if so default values set by selections' constructors are used.
        """
        if distinguish_chains is None:
            distinguish_chains = self._distinguish_chains
        try:
            list_of_pdb_ids = list(set.intersection(
                *map(set, map(operator.methodcaller('specify', structure_obj, distinguish_chains), self.selections))))
            # only mers present in all structures connected wiht selections are in list above; they are stored in order imposed by first selection
        except (TypeError, KeyError):
            raise ValueError("Not enough selections to intersect")
        return Set(list_of_pdb_ids, distinguish_chains)


class SelectionsRelativeComplement(CombinedSelection):
    """An relative complement of two selections.

    You can easily create relative complements subtracting two selections, e.g.:
    relComp_1_exc_2 = selection_1 - selection_2
    """

    def specify(self, structure_obj, distinguish_chains=None):
        """Creates user structure instance made of mers that meet first given selection criteria, but do not meet criteria given by secound selection.

        Arguments:
        structure_obj -- a structure to be base for new structure.
        distinguish_chains -- True or False; changes the behaviour of a new set in such a way that if the value of a chain is set to False, then the chain's character is ignored. Initially set to None, if so default values set by selections' constructors are used.
        """
        if distinguish_chains is None:
            distinguish_chains = self._distinguish_chains
        list_of_pdb_ids = [pdb_id for pdb_id in self.selections[0].specify(structure_obj, distinguish_chains) if
                           pdb_id not in self.selections[1].specify(structure_obj, distinguish_chains)]
        # adding mers to list, if they are not present in structure created under second selection contidions
        return Set(list_of_pdb_ids, distinguish_chains)
