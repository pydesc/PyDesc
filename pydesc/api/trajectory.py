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
"""Package providing convenience functions for dealing with trajectories."""

import mdtraj
import numpy

from pydesc.chemistry.factories import CopyingFactor
from pydesc.selection import Everything
from pydesc.selection import Selector
from pydesc.structure import TrajectoryLoader
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
    path_str = str(topology.path)
    converter = topology.converter

    md_traj = mdtraj.load_frame(path_str, 0)
    trajectory = Trajectory(name, path_str, converter, md_traj)
    n_frames = len(frames)
    n_atoms = md_traj.n_atoms
    coords = numpy.zeros((n_frames, n_atoms, 3))

    for n, frame in enumerate(frames):
        for mer in frame:
            for atom in mer:
                atom_index = trajectory.serial_map[atom.serial_number]
                coords[n, atom_index] = atom.vector
    trajectory.md_matrix = coords
    loader = TrajectoryLoader()
    chains = [loader.create_chain(chain, trajectory) for chain in topology.chains]
    trajectory.finalize(chains)

    return trajectory
