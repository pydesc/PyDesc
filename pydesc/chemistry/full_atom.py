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
from pydesc.chemistry.base import Ligand
from pydesc.chemistry.base import Mer
from pydesc.chemistry.base import Pseudoatom
from pydesc.chemistry.base import register_dynamic_feature
from pydesc.chemistry.base import register_pseudoatom
from pydesc.config import ConfigManager
from pydesc.warnexcept import IncompleteParticle
from pydesc.warnexcept import NoConfiguration
from pydesc.warnexcept import WrongAtomSetType
from pydesc.warnexcept import warn

norm = scipy.linalg.get_blas_funcs("nrm2")

ConfigManager.chemistry.new_branch("residue")
ConfigManager.chemistry.residue.set_default("bb_bond_threshold", 2.0)
ConfigManager.chemistry.residue.set_default("backbone_atoms", ("N", "CA", "C"))
ConfigManager.chemistry.residue.set_default("indicators", ("CA", "cbx"))
ConfigManager.chemistry.residue.set_default("adjusted_segment_length", 18.0)
ConfigManager.chemistry.residue.set_default("check_distances", False)
ConfigManager.chemistry.residue.set_default(
    "crucial_atom_distances", (("C", "CA", 1.35, 1.71), ("CA", "N", 1.35, 1.75))
)
ConfigManager.chemistry.residue.set_default(
    "code",
    {
        "ILE": "I",
        "GLN": "Q",
        "GLX": "Z",
        "GLY": "G",
        "GLU": "E",
        "CYS": "C",
        "HIS": "H",
        "SER": "S",
        "LYS": "K",
        "PRO": "P",
        "ASX": "B",
        "ASN": "N",
        "VAL": "V",
        "THR": "T",
        "ASP": "D",
        "TRP": "W",
        "PHE": "F",
        "ALA": "A",
        "MET": "M",
        "LEU": "L",
        "ARG": "R",
        "TYR": "Y",
    },
)
ConfigManager.chemistry.residue.set_default(
    "additional_code",
    {
        "DNP": "A",
        "ABI": "A",
        "ALM": "A",
        "MAA": "A",
        "TIH": "A",
        "FLA": "A",
        "DAL": "A",
        "CSD": "A",
        "BNN": "A",
        "HAC": "A",
        "PRR": "A",
        "AYA": "A",
        "CHG": "A",
        "DHA": "A",
        "TPQ": "A",
        "SEG": "A",
        "DIV": "V",
        "MVA": "V",
        "DVA": "V",
        "BUG": "L",
        "DLE": "L",
        "CLE": "L",
        "NLN": "L",
        "NLE": "L",
        "NLP": "L",
        "MLE": "L",
        "LEF": "L",
        "DIL": "I",
        "IIL": "I",
        "DPR": "P",
        "HYP": "P",
        "MSE": "M",
        "OMT": "M",
        "CXM": "M",
        "FME": "M",
        "MME": "M",
        "DAH": "F",
        "PHI": "F",
        "DPN": "F",
        "HPQ": "F",
        "PHL": "F",
        "LTR": "W",
        "TPL": "W",
        "DTR": "W",
        "TRO": "W",
        "HTR": "W",
        "MSA": "G",
        "SAR": "G",
        "MPQ": "G",
        "GLZ": "G",
        "GSC": "G",
        "GL3": "G",
        "NMC": "G",
        "DSN": "S",
        "SEL": "S",
        "SEP": "S",
        "SET": "S",
        "SAC": "S",
        "SVA": "S",
        "MIS": "S",
        "OAS": "S",
        "TPO": "T",
        "ALO": "T",
        "DTH": "T",
        "BMT": "T",
        "BCS": "C",
        "SOC": "C",
        "C5C": "C",
        "C6C": "C",
        "SCS": "C",
        "PEC": "C",
        "DCY": "C",
        "EFC": "C",
        "SCY": "C",
        "SMC": "C",
        "CSX": "C",
        "BUC": "C",
        "CSO": "C",
        "PR3": "C",
        "CCS": "C",
        "CEA": "C",
        "CME": "C",
        "CSP": "C",
        "CSS": "C",
        "CSW": "C",
        "CY1": "C",
        "CY3": "C",
        "CYG": "C",
        "CYM": "C",
        "CYQ": "C",
        "SCH": "C",
        "SHC": "C",
        "OCS": "C",
        "CAS": "C",
        "TYQ": "Y",
        "TYS": "Y",
        "TYB": "Y",
        "STY": "Y",
        "DTY": "Y",
        "IYR": "Y",
        "PAQ": "Y",
        "TYY": "Y",
        "PTR": "Y",
        "TYI": "Y",
        "MEN": "N",
        "DGN": "Q",
        "MGN": "Q",
        "2AS": "D",
        "ASB": "D",
        "DAS": "D",
        "ASK": "D",
        "ASL": "D",
        "ASQ": "D",
        "BHD": "D",
        "ASA": "D",
        "DSP": "D",
        "5HP": "E",
        "CGU": "E",
        "DGL": "E",
        "GMA": "E",
        "GGL": "E",
        "PCA": "E",
        "DLY": "K",
        "LYM": "K",
        "LLY": "K",
        "LYZ": "K",
        "KCX": "K",
        "LLP": "K",
        "TRG": "K",
        "SHR": "K",
        "ALY": "K",
        "ARM": "R",
        "ACL": "R",
        "HAR": "R",
        "HMR": "R",
        "AGM": "R",
        "DAR": "R",
        "HIC": "H",
        "3AH": "H",
        "NEM": "H",
        "NEP": "H",
        "DHI": "H",
        "MHS": "H",
        "HIP": "H",
    },
)
ConfigManager.chemistry.new_branch("nucleotide")
ConfigManager.chemistry.nucleotide.set_default("bb_bond_threshold", 2.0)
ConfigManager.chemistry.nucleotide.set_default(
    "code",
    {
        "G": "G",
        "C": "C",
        "U": "U",
        "A": "A",
        "DG": "G",
        "DA": "A",
        "DT": "T",
        "DC": "C",
    },
)
ConfigManager.chemistry.nucleotide.set_default(
    "backbone_atoms", ("P", "O5'", "C5'", "C4'", "C3'", "O3'")
)
ConfigManager.chemistry.nucleotide.set_default(
    "ring_atoms", ("N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9")
)
ConfigManager.chemistry.nucleotide.set_default("check_distances", False)
ConfigManager.chemistry.nucleotide.set_default(
    "crucial_atom_distances",
    (
        ("P", "O5'", 1.54, 1.66),
        ("O5'", "C5'", 1.34, 1.54),
        ("C5'", "C4'", 1.44, 1.56),
        ("C4'", "C3'", 1.46, 1.58),
        ("C3'", "O3'", 1.37, 1.49),
    ),
)
ConfigManager.chemistry.nucleotide.set_default(
    "indicators", ("C3'", "P", "ring_center")
)
ConfigManager.chemistry.new_branch("monoatomicion")
ConfigManager.chemistry.monoatomicion.set_default("indicators", ("gc",))
ConfigManager.chemistry.monoatomicion.set_default(
    "radii",
    {
        "BE": 0.59,
        "BA": 1.49,
        "BI": 1.17,
        "BK": 1.1,
        "BR": 1.82,
        "RU": 0.82,
        "RE": 0.77,
        "TM": 1.17,
        "RA": 1.62,
        "RB": 1.66,
        "RH": 0.805,
        "P": 0.58,
        "GE": 0.87,
        "GD": 1.078,
        "GA": 0.76,
        "OS": 0.77,
        "C": 0.3,
        "HO": 1.041,
        "HF": 0.85,
        "HG": 1.33,
        "PR": 1.13,
        "PT": 0.94,
        "PU": 1.14,
        "PB": 1.33,
        "PA": 1.16,
        "PD": 1.0,
        "PO": 1.08,
        "PM": 1.11,
        "ZN": 0.88,
        "K": 1.52,
        "O": 1.26,
        "S": 1.7,
        "W": 0.8,
        "EU": 1.31,
        "ZR": 0.86,
        "ER": 1.03,
        "MG": 0.86,
        "MO": 0.83,
        "MN": 0.97,
        "AU": 1.51,
        "FR": 1.94,
        "FE": 0.92,
        "NI": 0.83,
        "NA": 1.16,
        "NB": 0.86,
        "ND": 1.43,
        "ES": 0.928,
        "NP": 1.24,
        "B": 0.41,
        "CO": 0.885,
        "CM": 1.11,
        "CL": 1.67,
        "CA": 1.14,
        "CF": 1.09,
        "CE": 1.15,
        "N": 1.32,
        "V": 0.93,
        "CS": 1.81,
        "CR": 0.94,
        "CU": 0.91,
        "SR": 1.32,
        "SI": 0.54,
        "SN": 0.83,
        "SM": 1.36,
        "SC": 0.885,
        "SB": 0.9,
        "SE": 1.84,
        "YB": 1.16,
        "DY": 1.21,
        "LA": 1.172,
        "F": 1.19,
        "LI": 0.9,
        "TL": 1.64,
        "LU": 1.001,
        "TH": 1.08,
        "TI": 1.0,
        "TE": 2.07,
        "TB": 1.063,
        "TC": 0.785,
        "TA": 0.86,
        "AC": 1.26,
        "AG": 1.29,
        "I": 2.06,
        "IR": 0.82,
        "AM": 1.4,
        "AL": 0.675,
        "AS": 0.72,
        "U": 1.165,
        "AT": 0.76,
        "IN": 0.94,
        "Y": 1.04,
        "CD": 1.09,
        "XE": 0.62,
    },
)
ConfigManager.chemistry.new_branch("compound")
ConfigManager.chemistry.compound.set_default("indicators", ("gc",))


def calculate_residues_angles_vectorized(residues):
    """Calculate torsion angles for given residues.

    Sets "angles" property in residues.
    Purpose of this method is to speed up calculation of angles -- its faster
    than calling property in every residue.

    Args:
        residues: sequence of Residue instances.

    """
    n_res = len(residues)
    if n_res == 0:
        return
    n, ca, c = numpy.transpose(
        numpy.array([[a.vector for a in r.iter_bb_atoms()] for r in residues]),
        (1, 0, 2),
    )[[0, 1, 2]]
    pc = numpy.empty((n_res, 3), dtype=numpy.float32)
    nn = numpy.empty((n_res, 3), dtype=numpy.float32)

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
        t1 = numpy.nan_to_num(t2).astype(float)

        angs.append(t1)

    for res, (psi, phi) in zip(residues, list(zip(*angs))):
        res.dynamic_features["angles"] = (psi, phi)


class FullAtomMer(Mer):
    """Abstract superclass for full-atom Residues and Nucleotides."""

    @register_pseudoatom
    def rc(self):
        """Geometrical center of side chain (non-backbone atoms)."""
        non_backbone_coordinates = [a.vector for a in self.iter_nbb_atoms()]
        with numpy.errstate(divide="raise", invalid="raise"):
            vector = numpy.mean(non_backbone_coordinates, 0)
        return Pseudoatom(numpy_vec=vector, name="rc")


class FullAtomLigand(Ligand):
    """Abstract superclass for full-atom Ligands."""

    @property
    def rc(self):
        """The same as "gc"."""
        return self.gc


class Residue(FullAtomMer):
    """Representation of a residue.

    Implicit checks are done during initialization. See settings in branch
    "ConfigManager.chemistry.residue".

    Args:
        ind: residue id in pydesc structure.
        name: residue name.
        chain: residue chain name.
        atoms(dict): map of names to Atom instances.

    """

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
            direction = (self.atoms["CA"] - self.atoms["C"]).vector
            pl3 = pydesc.geometry.Plane.build(*([prm.atoms["C"]] + atoms[:2]))
            ang_phi = pl2.dihedral_angle(pl3, direction)

        if nxm is not None:
            direction = (self.atoms["N"] - self.atoms["CA"]).vector
            pl1 = pydesc.geometry.Plane.build(*(atoms[1:] + [nxm.atoms["N"]]))
            ang_psi = pl1.dihedral_angle(pl2, direction)

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

    @property
    def simple_secondary_structure(self):
        """Secondary structure in simple 3-letter code.

        H -- helix
        E -- extended strand
        C -- coil
        """
        temp = ConfigManager.structure_mon.simple_secondary_structure_code
        return temp[self._ss]


class Nucleotide(FullAtomMer):
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
            at2 = self.atoms["N9"]
        except KeyError:
            at2 = self.atoms["N1"]
        vec = (at2 - at1).vector
        nrm = norm(vec)
        nvec = vec * ((nrm + 1.4) / nrm)

        nx = at1.vector + nvec
        return Pseudoatom(numpy_vec=nx, name="nx")


class MonoatomicIon(FullAtomLigand):
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


class Compound(FullAtomLigand):
    """Representation of any chemical compound beside biopolymers."""

    def __init__(self, ind, name, chain, atoms):
        super(Compound, self).__init__(ind, name, chain, atoms)


AtomSet.reset_config_cache()
