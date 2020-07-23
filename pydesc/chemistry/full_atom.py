# Copyright 2019 Tymoteusz Oleniecki
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

"""Classes storing mers and ligands data in full-atom representation."""

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
    """Representation of a residue.

    Implicit checks are done during initialization. See settings in branch
    "ConfigManager.chemistry.residue".

    Args:
        ind: residue id in pydesc structure.
        name: residue name.
        chain: residue chain name.
        atoms(dict): map of names to Atom instances.

    """

    @staticmethod
    def calculate_angles_static(structure_obj):
        """Calculate torsion angles of residues in given structure.

        Only subclasses of Residue are taken into account.

        Sets "angles" property in residues.
        Purpose of this method is to speed up calculation of angles -- its faster
        than calling property in every residue.

        Args:
            structure_obj: sequence of AtomSet instances.

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

        no_prev = numpy.fromiter((r.prev_mer is None for r in residues), dtype=bool)
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
            res.dynamic_features["angles"] = (psi, phi)

    def __init__(self, ind, name, chain, atoms):
        super().__init__(ind, name, chain, atoms)

    @register_pseudoatom
    def backbone_average(self):
        """Pseudoatom; moving average CA.

        Moving average frame is configurable (see ConfigManager).
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
                last_mer = last_mer.prev_mer
                average_ca += next_mer.ca.vector + last_mer.ca.vector
                cnt += 2
        except AttributeError:
            # AttributeError is raised by mers at the beginning and at the
            # end of chain they have no next/previous mers
            pass

        return Pseudoatom(numpy_vec=(average_ca / cnt), name="bb_average")

    @register_dynamic_feature
    def angles(self):
        """Torsion angles psi and phi."""
        ang_psi, ang_phi = 0.0, 0.0

        prm = self.prev_mer
        nxm = self.next_mer

        atoms = [self.atoms["N"], self.atoms["CA"], self.atoms["C"]]

        pl2 = pydesc.geometry.Plane.build(*atoms)

        if prm is not None:
            pl3 = pydesc.geometry.Plane.build(*([prm.atoms["C"]] + atoms[:2]))
            ang_phi = pl2.dihedral_angle(pl3)

        if nxm is not None:
            pl1 = pydesc.geometry.Plane.build(*(atoms[1:] + [nxm.atoms["N"]]))
            ang_psi = pl1.dihedral_angle(pl2)

        # TODO: sign is sometimes wrong -- check why
        return ang_psi, ang_phi

    @property
    def ca(self):
        """Carbon alpha."""
        return self.atoms["CA"]

    @register_pseudoatom
    def cbx(self):
        """Pseudoatom; CA->CB vector extended by 1A (fixed for GLY)."""
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

    def get_adjusted_length(self):
        """Get distance between backbone_average pseudoatoms of this and
        the next mer or None if distance cannot be computed."""
        try:
            return abs(self.backbone_average - self.next_mer.backbone_average)
        except AttributeError:
            return None


class Nucleotide(Mer):
    """Representation of a nucleotide.

    Implicit checks are done during initialization. See settings in branch
    "ConfigManager.chemistry.nucleotide".

    Args:
        ind: residue id in pydesc structure.
        name: residue name.
        chain: residue chain name.
        atoms(dict): map of names to Atom instances.

    """

    def __init__(self, ind, name, chain, atoms):
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
        """Calculate and return pseudoatom representing center of the base ring being
        closer to glycosidic bond."""
        vec = numpy.array([0.0, 0.0, 0.0])
        for at in ("C4", "C5", "N7", "C8", "N9"):
            vec += self.atoms[at].vector
        vec /= 5.0
        return Pseudoatom(numpy_vec=vec, name="prc")

    @register_pseudoatom
    def prc(self):
        """Pseudoatom; center of base ring closest to sugar."""
        try:
            return self.calculate_proximate_ring_center()
        except KeyError:
            return self.ring_center

    @register_pseudoatom
    def ring_center(self):
        """Pseudoatom; (larger) base ring center."""
        try:
            vec = (self.ring_atoms["N1"].vector + self.ring_atoms["C4"].vector) * 0.5
        except KeyError:
            raise IncompleteParticle("Lacking N1 or C4, cannot calculate ring center.")
        return Pseudoatom(numpy_vec=vec, name="ring_center")

    @register_dynamic_feature
    def ring_plane(self):
        """Dynamic feature; base plane."""
        at1, at2, at3 = (
            self.ring_atoms["C2"],
            self.ring_atoms["C4"],
            self.ring_atoms["C6"],
        )
        return pydesc.geometry.Plane.build(at1, at2, at3)

    @register_pseudoatom
    def nx(self):
        """Pseudoatom; vector along glycosidic bond extended by 1.4A ."""
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
