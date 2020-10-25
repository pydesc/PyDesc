import pytest

from pydesc.api.trajectory import freeze_frame
from pydesc.api.trajectory import from_frames
from pydesc.chemistry.factories import CopyingFactor
from pydesc.chemistry.full_atom import Nucleotide
from pydesc.chemistry.full_atom import Residue
from pydesc.selection import AtomSetSubclass
from pydesc.selection import Selector
from pydesc.selection import Set
from pydesc.structure import StructureLoader
from pydesc.structure import TrajectoryLoader
from pydesc.structure.topology import Structure
from pydesc.structure.trajectory import Trajectory
from pydesc.warnexcept import FrameNotAvailable


@pytest.mark.parametrize("ext", ["dcd", "xtc"])
@pytest.mark.parametrize("stc", ["mdm2", "trf1-dna"])
class TestTrajectory:
    trajs = {}
    stcs = {}

    @classmethod
    def load_stc(cls, stc):
        try:
            return cls.stcs[stc]
        except KeyError:
            pass
        structure_loader = StructureLoader()
        structures = structure_loader.load_structures(
            path="tests/data/test_trajectories/topologies/%s.pdb" % stc
        )
        cls.stcs[stc] = structures[0]
        return structures[0]

    @classmethod
    def load_traj(cls, ext, stc):
        try:
            return cls.trajs[(ext, stc)]
        except KeyError:
            pass
        structure = cls.load_stc(stc)
        trajectory_loader = TrajectoryLoader()
        pth = "tests/data/test_trajectories/%s/%s_5frames.%s"
        trajectory = trajectory_loader.load_trajectory(pth % (ext, stc, ext), structure)
        cls.trajs[(ext, stc)] = trajectory
        return trajectory

    def test_loader(self, ext, stc):
        trajectory = self.load_traj(ext, stc)
        assert isinstance(trajectory, Trajectory)
        assert trajectory.get_n_frames() == 5

    def test_residue(self, ext, stc):
        trajectory = self.load_traj(ext, stc)
        for mer in trajectory:
            if isinstance(mer, Residue):
                break
        else:
            pytest.xfail()
            return
        coords1 = mer.atoms["CA"].vector
        pseudo_coords1 = mer.cbx

        trajectory.set_frame(1)

        pseudo_coords2 = mer.cbx
        coords2 = mer.atoms["CA"].vector
        assert "cbx" in mer.pseudoatoms

        trajectory.set_frame(0)

        coords3 = mer.atoms["CA"].vector

        assert tuple(coords1) != tuple(coords2)
        assert tuple(coords1) == tuple(coords3)
        assert tuple(pseudo_coords1) != tuple(pseudo_coords2)
        assert mer.pseudoatoms == {}

    def test_nucleotide(self, ext, stc):
        trajectory = self.load_traj(ext, stc)
        for mer in trajectory:
            if isinstance(mer, Nucleotide):
                break
        else:
            pytest.xfail()
            return
        coords1 = mer.atoms["P"].vector
        pseudo_coords1 = mer.prc

        trajectory.set_frame(1)

        pseudo_coords2 = mer.prc
        coords2 = mer.atoms["P"].vector
        assert "prc" in mer.pseudoatoms

        trajectory.set_frame(0)

        coords3 = mer.atoms["P"].vector

        assert tuple(coords1) != tuple(coords2)
        assert tuple(coords1) == tuple(coords3)
        assert tuple(pseudo_coords1) != tuple(pseudo_coords2)
        assert mer.pseudoatoms == {}

    def test_frames(self, ext, stc):
        trajectory = self.load_traj(ext, stc)

        for i in range(5):
            trajectory.set_frame(i)

        with pytest.raises(FrameNotAvailable):
            trajectory.set_frame(42)

    def test_selection(self, ext, stc):
        trajectory = self.load_traj(ext, stc)
        subclass = type(trajectory[0])
        selection = AtomSetSubclass(subclass)
        new_stc = selection.create_structure(trajectory)
        assert len(new_stc) > 0
        assert trajectory[0] in new_stc

        trajectory.set_frame(2)

        frame2_atom_coords = tuple(trajectory[0])[0].vector
        selection_atom_coords = tuple(new_stc[0])[0].vector
        assert tuple(frame2_atom_coords) == tuple(selection_atom_coords)

    def test_freezing_substructure(self, ext, stc):
        trajectory = self.load_traj(ext, stc)
        get_id = trajectory.converter.get_pdb_id

        pdb_inds = [get_id(mer.ind) for mer in trajectory[:5]]
        set_sel = Set(pdb_inds)

        picker = Selector(CopyingFactor())
        picker.create_new_structure(set_sel, trajectory)

    def test_freezing_structure(self, ext, stc):
        trajectory = self.load_traj(ext, stc)

        trajectory.set_frame(2)

        new_structure = freeze_frame(trajectory)

        assert isinstance(new_structure, Structure)
        for traj_chain, new_chain in zip(trajectory.chains, new_structure.chains):
            assert traj_chain.name == new_chain.name
            assert len(traj_chain) == len(new_chain)
        assert len(new_structure) == len(trajectory)

        mer0_atoms = new_structure[0].atoms
        assert len(mer0_atoms) >= 1

        for atom in mer0_atoms:
            new_vector = new_structure[0].atoms[atom].vector
            old_vector = trajectory[0].atoms[atom].vector
            assert (new_vector == old_vector).all()

        trajectory.set_frame(1)

        for atom in mer0_atoms:
            new_vector = new_structure[0].atoms[atom].vector
            old_vector = trajectory[0].atoms[atom].vector
            assert (new_vector != old_vector).all()

    def test_assembling_frozen_frames(self, ext, stc):
        trajectory = self.load_traj(ext, stc)
        frames = []
        for index in range(5):
            trajectory.set_frame(index)
            frame = freeze_frame(trajectory)
            frames.append(frame)

        new_trajectory = from_frames(frames)

        assert trajectory.get_n_frames() == new_trajectory.get_n_frames()

        assert len(new_trajectory) == len(trajectory)
        assert len(new_trajectory.chains) == len(trajectory.chains)
        old_chains = trajectory.chains
        new_chains = new_trajectory.chains
        for old_chain, new_chain in zip(old_chains, new_chains):
            assert len(new_chain) == len(old_chain)

        for index in range(5):
            trajectory.set_frame(index)
            new_trajectory.set_frame(index)
            assert (trajectory.md_matrix == new_trajectory.md_matrix).all()
