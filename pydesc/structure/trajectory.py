from pydesc.structure.topology import Structure
from pydesc.warnexcept import FrameNotAvailable


class Trajectory(Structure):
    """Extension of Structure meant to consist of objects containing AtomProxy
    instances to given MDTraj trajectory."""

    def __init__(self, name, path, converter_obj, md_trajectory):
        """Initialize Trajectory, extended Structure, with superclass arguments plus
        MDTraj trajectory."""
        super().__init__(name, path, converter_obj)
        self._frame = 0
        self.md_matrix = md_trajectory.xyz
        atoms = md_trajectory.topology.atoms
        self.serial_map = {atom.serial: atom.index for atom in atoms}

    def get_atom_coords(self, atom):
        """Get coords of given proxy atom."""
        atom_index = self.serial_map[atom.serial_number]
        return self.md_matrix[self._frame, atom_index]

    def set_frame(self, n):
        """Change trajectory frame to given one."""
        if not (0 <= n <= self.md_matrix.shape[0]):
            raise FrameNotAvailable("Frame %i out of range." % n)
        if n != self._frame:
            for mer in self:
                mer.reset_dynamic_cache()
        self._frame = n

    def get_frame(self):
        """Return currently set frame."""
        return self._frame

    def get_n_frames(self):
        """Get number of frames in trajectory."""
        return self.md_matrix.shape[0]

    def __repr__(self):
        return "<Trajectory %s>" % self.name
