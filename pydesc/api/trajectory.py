import mdtraj
import numpy

from pydesc.mers.factories import CopyingFactor
from pydesc.selection import Everything
from pydesc.selection import Selector
from pydesc.structure.topology import Chain
from pydesc.structure.topology import Structure
from pydesc.structure.trajectory import Trajectory


def freeze_frame(trajectory):
    """Return current frame of given trajectory as new structure."""
    factory = CopyingFactor()
    selection = Everything()
    selector = Selector(factory)
    frame = Structure(trajectory.name, trajectory.path, trajectory.converter)
    chains = []
    for chain in trajectory.chains:
        chain_mers = selector.create_new_structure(selection, chain)
        new_chain = Chain(frame, chain.chain_name, chain_mers)
        chains.append(new_chain)
    frame.finalize(chains)
    return frame


def from_frames(frames, topology=None):
    """Create trajectory from list of frames as separate structures.

    Note that all frames have to contain all mers and atoms that topology contains.
    If the come from different files -- first frame will be treated as topology.

    """
    if topology is None:
        topology = frames[0]
    name = topology.name
    path = topology.path
    converter = topology.converter

    md_traj = mdtraj.load_frame(path, 0)
    trajectory = Trajectory(name, path, converter, md_traj)
    n_frames = len(frames)
    n_atoms = md_traj.n_atoms
    coords = numpy.zeros((n_frames, n_atoms, 3))

    for n, frame in enumerate(frames):
        for mer in frame:
            for atom in mer:
                atom_index = trajectory.serial_map[atom.serial_number]
                coords[n, atom_index] = atom.vector
    trajectory.md_matrix = coords

    return trajectory
