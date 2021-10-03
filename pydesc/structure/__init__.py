from pathlib import Path

import mdtraj

from pydesc import warnexcept
from pydesc.chemistry.factories import BioPythonAtomSetFactory
from pydesc.chemistry.factories import MDTrajAtomSetFactory
from pydesc.config import ConfigManager
from pydesc.numberconverter import NumberConverterFactory
from pydesc.parsers import MetaParser
from pydesc.structure.topology import Chain
from pydesc.structure.topology import Structure
from pydesc.structure.trajectory import Trajectory

# pylint: disable=no-member
ConfigManager.new_branch("structure")
ConfigManager.structure.set_default("dssp_path", "dssp")


# pylint: enable=no-member


class StructureLoader:
    """Loads structures from the databases using a given designation.

    Args:
        parser: file parser returning BioPython's structures.
        atom_set_factory: factory able to produce AtomSet instances from BioPython
            residues.

    """

    def __init__(
        self, parser=MetaParser(QUIET=True), atom_set_factory=BioPythonAtomSetFactory(),
    ):
        self.parser = parser
        self.atom_set_factory = atom_set_factory

    def load_structures(self, file_handlers, common_converter=False):
        """Returns a list of Structure instances and the NumberConverter
        instance.

        Args:
            file_handlers: sequence of file-like objects to load structures from.
            common_converter(bool): determines if all models should have common id
                converter. If set to True, inds for all residues with the same
                PDB id will also be the same (cross-model, cross-file).
                Meant to facilitate work with different conformations of the same
                molecule or set of mutants with short inserts.
                False by default.

        Returns:
            list: Structure instances.

        """

        def get_converter():
            if not common_converter:
                return nc_factory.from_pdb_models([model])
            return converter

        self._assert_is_list(file_handlers)
        warnexcept.set_filters()
        models, paths = self._get_models_n_paths(file_handlers)
        nc_factory = NumberConverterFactory()
        if common_converter:
            converter = nc_factory.from_pdb_models(models)
        structures = []
        for model, path in zip(models, paths):
            converter = get_converter()
            structure = self.create_structure(model, path, converter)
            structures.append(structure)

        return structures

    @staticmethod
    def _assert_is_list(handlers):
        try:
            n_files = len(handlers)
        except TypeError:
            msg = (
                f"StructureLoader takes list of file handlers as an argument."
                f"Make sure file-like objects are in sequence, even if there is "
                f"only one object."
            )
            raise ValueError(msg)
        else:
            if n_files < 1:
                msg = f"At least 1 file handler required, {n_files} passed."
                raise ValueError(msg)

    def _get_models_n_paths(self, files):
        """Return list of models read from given PDB files."""
        models = []
        paths = []
        for handler in files:
            path = self._get_path(handler)
            name = self._get_name(path)
            with handler:
                pdb_structure = self.parser.get_structure(name, handler)
            new_models = [pdb_model for pdb_model in pdb_structure]
            models.extend(new_models)
            paths.extend([path] * len(new_models))
        return models, paths

    @staticmethod
    def _get_path(file_handler):
        file_path = Path(file_handler.name)
        return file_path

    @staticmethod
    def _get_name(file_path):
        """Return structure name."""
        name = file_path.stem
        return name

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
            raise ValueError("Impossible to load some sets of atoms.")

        mers = []
        hits = dict((klass, 0) for klass in self.atom_set_factory.chainable)

        for pdb_residue in pdb_chain:
            mer_dct, warns = self.atom_set_factory.create(
                pdb_residue=pdb_residue, structure_obj=structure, warn_in_place=False
            )
            if mer_dct is None:
                continue

            mers.append((mer_dct, warns))
            for mer_class in hits:
                hits[mer_class] += int(mer_class in mer_dct)

        try:
            winner = max(hits, key=lambda hit: hits[hit])
        except ValueError:
            winner = None
        others = self.atom_set_factory.other
        chain_mers = []
        for mer_dct, warns in mers:
            accepted_mer = pick_mer(mer_dct, winner, others)
            chain_mers.append(accepted_mer)
            warns.raise_all(type(accepted_mer))

        return Chain(structure, pdb_chain.get_id(), chain_mers)


class TrajectoryLoader:
    """Trajectory loading class."""

    def __init__(self, mer_proxy_factory=None):
        """Initialize loader with mer_proxy_factory creating residues with proxy
        atoms. If None is given (default) -- assumes factory able to copy MDTraj mer
        factory."""
        if mer_proxy_factory is None:
            mer_proxy_factory = MDTrajAtomSetFactory()
        self.mer_factory = mer_proxy_factory

    def load_trajectory(self, trajectory_path, structure_obj):
        """Create Trajectory instance loaded from given path for topology in given
        structure instance."""
        md_trajectory = mdtraj.load(trajectory_path, top=str(structure_obj.path))
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
