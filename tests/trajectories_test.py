import pytest

from pydesc.mers.factories import BioPythonMerFactory
from pydesc.mers.full_atom import Nucleotide
from pydesc.mers.full_atom import Residue
from pydesc.selection import MerSubclasses
from pydesc.selection import Selector
from pydesc.selection import Set
from pydesc.structure import StructureLoader
from pydesc.structure import TrajectoryLoader
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
        coords1 = mer.CA.vector
        pseudo_coords1 = mer.cbx

        trajectory.set_frame(1)

        pseudo_coords2 = mer.cbx
        coords2 = mer.CA.vector
        assert "cbx" in mer.pseudoatoms

        trajectory.set_frame(0)

        coords3 = mer.CA.vector

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
        coords1 = mer.P.vector
        pseudo_coords1 = mer.prc

        trajectory.set_frame(1)

        pseudo_coords2 = mer.prc
        coords2 = mer.P.vector
        assert "prc" in mer.pseudoatoms

        trajectory.set_frame(0)

        coords3 = mer.P.vector

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
        selection = MerSubclasses(subclass)
        new_stc = selection.create_structure(trajectory)
        assert len(new_stc) > 0
        assert trajectory[0] in new_stc

        trajectory.set_frame(2)

        frame2_atom_coords = tuple(trajectory[0])[0].vector
        selection_atom_coords = tuple(new_stc[0])[0].vector
        assert tuple(frame2_atom_coords) == tuple(selection_atom_coords)

    def test_freezing_substructure(self, ext, stc):
        trajectory = self.load_traj(ext, stc)

        pdb_inds = [mer.get_pdb_id() for mer in trajectory[:5]]
        set_sel = Set(pdb_inds)

        picker = Selector(BioPythonMerFactory)
        picker.create_new_structure(set_sel, trajectory)
