"""
Classes that represents molecular structures and their substructures.

created: 10.07.2013 - , Tymoteusz 'hert' Oleniecki
"""

import os.path

import mdtraj

import pydesc.dbhandler
import pydesc.geometry
from pydesc import warnexcept
from pydesc.config import ConfigManager
from pydesc.chemistry.factories import BioPythonMerFactory
from pydesc.chemistry.factories import MDTrajMerFactory
from pydesc.numberconverter import NumberConverterFactory
from pydesc.structure.topology import Chain
from pydesc.structure.topology import Structure
from pydesc.structure.trajectory import Trajectory
from pydesc.structure.trajectory import Trajectory

# pylint: disable=no-member
ConfigManager.new_branch("element")
ConfigManager.element.set_default("element_chainable_length", 5)
ConfigManager.new_branch("structure")
ConfigManager.structure.set_default("dssp_path", "dssp")


# pylint: enable=no-member


class StructureLoader:
    """Loads structures from the databases using a given designation."""

    def __init__(
        self,
        handler=pydesc.dbhandler.MetaHandler(),
        parser=pydesc.dbhandler.MetaParser(QUIET=True),
        mer_factory=BioPythonMerFactory(),
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

    def create_chain(self, pdb_chain, structure):
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
            mer_dct, warns = self.mer_factory.create(
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

    def create_structure(self, model, path, converter):
        """Return structure for given model

        Argument:
            model -- Bio.PDB model.
            path -- path to pdb file storing the structure.
            converter -- id converter for structure to be created.
        """
        name = model.get_full_id()[0]
        structure = Structure(name, path, converter)

        chains = [self.create_chain(pdb_chain, structure) for pdb_chain in model]

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
        To choose specific BioUnit type "<bio unit_ ode>/<number>", e.g.
        "Unit://1no5/1"
        in case of MetaHandler or
        "1no5/1"
        in case of UnitHandler.
        """
        warnexcept.set_filters()
        path, open_files = self._get_files_and_path(code, path)
        code = self._get_code(code, path)
        models = self._get_models(open_files, code)
        converter = NumberConverterFactory().from_pdb_models(models)
        structures = [self.create_structure(model, path, converter) for model in models]

        return structures


class TrajectoryLoader:
    """Trajectory loading class."""

    def __init__(self, mer_proxy_factory=None):
        """Initialize loader with mer_proxy_factory creating residues with proxy
        atoms. If None is given (default) -- assumes factory able to copy MDTraj mer
        factory."""
        if mer_proxy_factory is None:
            mer_proxy_factory = MDTrajMerFactory()
        self.mer_factory = mer_proxy_factory

    def load_trajectory(self, trajectory_path, structure_obj):
        """Create Trajectory instance loaded from given path for topology in given
        structure instance."""
        md_trajectory = mdtraj.load(trajectory_path, top=structure_obj.path)
        trajectory_obj = self.create_trajectory(md_trajectory, structure_obj)
        return trajectory_obj

    def create_trajectory(self, mdtraj_obj, structure_obj):
        """Create Trajectory instance from given MDTraj.Trajectory reading
        topology information from given pydesc structure instance."""
        trajectory_obj = Trajectory(
            structure_obj.name, structure_obj.path, structure_obj.converter, mdtraj_obj
        )
        chains = [
            self.create_chain(chain, trajectory_obj) for chain in structure_obj.chains
        ]
        trajectory_obj.finalize(chains)
        return trajectory_obj

    def create_chain(self, chain, trajectory_obj):
        """Copy given chain into chain containing proxy atoms to given trajectory
        object."""
        mers = [self.mer_factory.create(mer, trajectory_obj) for mer in chain]
        return Chain(trajectory_obj, chain.name, mers)
