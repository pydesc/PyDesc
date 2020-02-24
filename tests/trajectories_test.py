import pytest

from pydesc.structure import TrajectoryLoader
from pydesc.structure import StructureLoader

from pydesc.warnexcept import FrameNotAvailable


def test_trajectory_loader_dcd():
    structure_loader = StructureLoader()
    structures = structure_loader.load_structures(
        path="tests/data/test_trajectories/dcd/mdm2.pdb"
    )

    trajectory_loader = TrajectoryLoader()
    trajectory = trajectory_loader.load_trajectory(
        "tests/data/test_trajectories/dcd/mdm2_5frames.dcd", structures[0]
    )

    assert trajectory.get_n_frames() == 5

    coords1 = trajectory[3].CA.vector
    pseudo_coords1 = trajectory[3].cbx
    trajectory.set_frame(1)
    pseudo_coords2 = trajectory[3].cbx
    coords2 = trajectory[3].CA.vector
    trajectory.set_frame(0)
    coords3 = trajectory[3].CA.vector

    assert tuple(coords1) != tuple(coords2)
    assert tuple(coords1) == tuple(coords3)
    assert tuple(pseudo_coords1) != tuple(pseudo_coords2)
    assert trajectory[3].pseudoatoms == {}

    with pytest.raises(FrameNotAvailable):
        trajectory.set_frame(42)
