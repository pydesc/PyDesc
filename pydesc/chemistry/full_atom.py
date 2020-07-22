import numpy
import scipy.linalg

import pydesc.geometry
from pydesc.chemistry.base import AtomSet
from pydesc.chemistry.base import Mer
from pydesc.chemistry.base import Ligand
from pydesc.chemistry.base import Pseudoatom
from pydesc.chemistry.base import register_dynamic_feature
from pydesc.chemistry.base import register_pseudoatom
from pydesc.warnexcept import IncompleteParticle
from pydesc.warnexcept import NoConfiguration
from pydesc.warnexcept import warn
from pydesc.warnexcept import WrongAtomSetType

norm = scipy.linalg.get_blas_funcs("nrm2")


class Residue(Mer):
    """Representation of a residue."""

    @staticmethod
    def calculate_angles_static(structure_obj):
        """Calculates all torsion angles of residues in given structure.

        Argument:
        structure_obj -- instance of AbstractStructure subclass.

        Fills 'angles' property in all residues in given (sub)structure.
        Calculates them using numpy, much faster than non-static Residue
        method.
        """
        residues = [mer for mer in structure_obj if isinstance(mer, Residue)]
        nres = len(residues)
        if nres == 0:
            return
        n, ca, c = numpy.transpose(
            numpy.array([[a.vector for a in r.iter_bb_atoms()] for r in residues]),
            (1, 0, 2),
        )[[0, 1, 2]]
        pc = numpy.empty((nres, 3), dtype=numpy.float32)
        nn = numpy.empty((nres, 3), dtype=numpy.float32)

        pc[1:] = c[:-1]
        nn[:-1] = n[1:]

        no_prev = numpy.fromiter((r.previous_mer is None for r in residues), dtype=bool)
        no_next = numpy.fromiter((r.next_mer is None for r in residues), dtype=bool)

        pc[no_prev] = n[no_prev]
        nn[no_next] = c[no_next]

        cca = c - ca
        cnn = c - nn
        nca = n - ca
        npc = n - pc

        pl1 = numpy.cross(cca, cnn)  # vectors perpendicular to plane 1
        pl2 = numpy.cross(nca, cca)  # vectors perpendicular to plane 2
        pl3 = numpy.cross(npc, nca)  # vectors perpendicular to plane 3

        with numpy.errstate(divide="ignore", invalid="ignore"):
            pl1, pl2, pl3 = (
                pl / numpy.sqrt(numpy.einsum("ij,ij->i", pl, pl)).reshape(-1, 1)
                for pl in (pl1, pl2, pl3)
            )

        angs = []
        for planes, direction in (((pl1, pl2), -cca), ((pl2, pl3), nca)):
            cos = numpy.einsum("ij,ij->i", *planes)
            cpr = numpy.cross(*planes)
            sin = numpy.sqrt(numpy.einsum("ij,ij->i", cpr, cpr))
            sign = numpy.sign(numpy.einsum("ij,ij->i", direction, cpr))

            t2 = numpy.arctan2(sin, cos) * sign
            t1 = numpy.nan_to_num(t2)

            angs.append(t1)

        for res, (psi, phi) in zip(residues, list(zip(*angs))):
            res.dynamic_properties["angles"] = (psi, phi)

    def __init__(self, ind, name, chain, atoms):
        """Residue initializer.

        Arguments:
        ind -- mers index.
        name -- mers name.
        chain -- mers chain name (str).
        atoms -- dict mapping atoms names to Atom instances.

        Raises Warning if a given pdb_residue does not contain proper atoms
        or if its atoms occur in wrong distances.
        Extended MonomerChainable method.
        See also config file docstring.

        Config parameters in branch ConfigManager.mers.residue:
        min_c_ca_dist
        max_c_ca_dist
        min_ca_n_dist
        max_ca_n_dist
        min_c_o_dist
        max_c_o_dist
        old_cbx_calculation -- True or False
        """
        super().__init__(ind, name, chain, atoms)

    @register_pseudoatom
    def backbone_average(self):
        """Calculates coordinates of average ca pseudoatom and adds it to
        current residue pseudoatoms.

        Average ca is calculated as moving average for configurable number
        of residues around current residue.
        """
        steps = self.get_config("moving_average")
        if not steps % 2 == 1:
            raise ValueError("Wrong Number of steps for moving average.")
        average_ca = numpy.array(self.ca.vector)
        next_mer = last_mer = self
        cnt = 1
        try:
            for _ in range(steps // 2):
                next_mer = next_mer.next_mer
                last_mer = last_mer.previous_mer
                average_ca += next_mer.ca.vector + last_mer.ca.vector
                cnt += 2
        except AttributeError:
            # AttributeError is raised by mers at the beginning and at the
            # end of chain they have no next/previous mers
            pass

        return Pseudoatom(numpy_vec=(average_ca / cnt), name="bb_average")

    @register_dynamic_feature
    def angles(self):
        """Calculates torsion angles of residue and fills 'angles' property."""
        ang_psi, ang_phi = 0.0, 0.0

        prm = self.previous_mer
        nxm = self.next_mer

        atoms = [self.atoms["N"], self.atoms["CA"], self.atoms["C"]]

        pl2 = pydesc.geometry.Plane.build(*atoms)

        if prm is not None:
            pl3 = pydesc.geometry.Plane.build(*([prm.atoms["C"]] + atoms[:2]))
            ang_phi = pl2.dihedral_angle(pl3)

        if nxm is not None:
            pl1 = pydesc.geometry.Plane.build(*(atoms[1:] + [nxm.atoms["N"]]))
            ang_psi = pl1.dihedral_angle(pl2)

        return ang_psi, ang_phi

    @property
    def ca(self):
        """Property that returns current residue alpha carbon Atom object."""
        return self.atoms["CA"]

    @register_pseudoatom
    def cbx(self):
        """Adds Pseudoatom containing coordinates of the point that lies 1A
        farther from carbon alpha, than does carbon beta; or carbon alpha
        coordinates for GLY.
        """
        if self.name == "GLY":
            n_2_ca = self.atoms["N"] - self.atoms["CA"]
            c_2_ca = self.atoms["C"] - self.atoms["CA"]
            average_ca_cb_distance = 1.53
            cbx = (n_2_ca + c_2_ca).get_unit_vector() * (average_ca_cb_distance + 1)
            cbx = self.ca.vector + cbx
        else:
            try:
                ca = self.atoms["CA"].vector
                cb = self.atoms["CB"].vector
            except KeyError:
                msg = "AtomSet lacks CA or CB, cannot calculate cbx."
                raise IncompleteParticle(msg)
            vec = cb - ca
            nrm = norm(vec)
            vec = vec * ((nrm + 1) / nrm)
            cbx = ca + vec

        return Pseudoatom(numpy_vec=cbx, name="cbx")

    @register_pseudoatom
    def rc(self):
        """Return pseudoatom storing side chain geometric center."""
        return super().rc


class Nucleotide(Mer):
    """Representation of a nucleotide."""

    def __init__(self, ind, name, chain, atoms):
        """Nucleotide constructor.

        Arguments:
        pdb_residue -- BioPython Bio.PDB.Residue instance based on which
        Nucleotide is being created.
        structure -- the Structure instance to which Nucleotide belongs.

        Raises Warning if given pdb_residue does not contain proper atoms or if
        its atoms occur in wrong distances.
        Extended MonomerChainable method.
        See also config file docstring.

        Config parameters in branch ConfigManager.mers.nucleotide:
        min_o5'_p_dist
        max_o5'_p_dist
        min_c5'_o5'_dist
        max_c5'_o5'_dist
        min_c4'_c5'_dist
        max_c4'_c5'_dist
        min_c3'_c4'_dist
        max_c3'_c4'_dist
        min_c3'_o3'_dist
        max_c3'_o3'_dist
        """
        super().__init__(ind, name, chain, atoms)

        rats = self.get_config("ring_atoms")

        def flag(name, atom):
            if name in rats:
                atom.ring_flag = True
                return True
            atom.ring_flag = False

        self.ring_atoms = {
            name: atom for name, atom in list(self.atoms.items()) if flag(name, atom)
        }

    def calculate_proximate_ring_center(self):
        """Adds pseudoatom representing center of the base ring being closer to
        glycosidic bond."""
        vec = numpy.array([0.0, 0.0, 0.0])
        for at in ("C4", "C5", "N7", "C8", "N9"):
            vec += self.atoms[at].vector
        vec /= 5.0
        return Pseudoatom(numpy_vec=vec, name="prc")

    @register_pseudoatom
    def prc(self):
        """Get ring center of base ring closest to sugar."""
        try:
            return self.calculate_proximate_ring_center()
        except KeyError:
            return self.ring_center

    @register_dynamic_feature
    def ring_center(self):
        """Adds pseudoatom representing base ring center."""
        try:
            vec = (self.ring_atoms["N1"].vector + self.ring_atoms["C4"].vector) * 0.5
        except KeyError:
            raise IncompleteParticle("Lacking N1 or C4, unable to create Nucleotide.")
        return Pseudoatom(numpy_vec=vec, name="ring_center")

    @register_dynamic_feature
    def ring_plane(self):
        """Adds pydesc.geometry.Plane object representing base to current
        nucleotide pseudoatom dictionary."""
        at1, at2, at3 = (
            self.ring_atoms["C2"],
            self.ring_atoms["C4"],
            self.ring_atoms["C6"],
        )
        return pydesc.geometry.Plane.build(at1, at2, at3)

    @register_pseudoatom
    def nx(self):
        """Adds pseudoatom representing extended by 1.4A vector along glycosidic
        bond."""
        at1 = self.atoms["C1'"]
        try:
            at2 = self.N9
        except AttributeError:
            at2 = self.N1
        vec = (at2 - at1).vector
        nrm = norm(vec)
        nvec = vec * ((nrm + 1.4) / nrm)

        nx = at1.vector + nvec
        return Pseudoatom(numpy_vec=nx, name="nx")


class MonoatomicIon(Ligand):
    """Representation of an monoatomic ion ligand."""

    def __init__(self, ind, name, chain, atoms):
        super(MonoatomicIon, self).__init__(ind, name, chain, atoms)
        if len(self.atoms) != 1:
            msg = "Multiple atoms passed to initialization."
            raise WrongAtomSetType(msg)

    def get_radius(self):
        """Return ion radius."""
        name = max(self.atoms)
        try:
            return self.get_config("radii")[name]
        except KeyError:
            warn(NoConfiguration("No radius for %s ions." % name))
            return 2.5


class Compound(Ligand):
    """Representation of any chemical compound beside biopolymers."""

    def __init__(self, ind, name, chain, atoms):
        super(Compound, self).__init__(ind, name, chain, atoms)


AtomSet.reset_config_cache()
