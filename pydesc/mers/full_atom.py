import numpy
import scipy.linalg

import pydesc.geometry
from pydesc.warnexcept import IncompleteParticle
from pydesc.warnexcept import NoConfiguration
from pydesc.warnexcept import warn
from pydesc.warnexcept import WrongMerType
from .base import Mer
from .base import MerChainable
from .base import MerOther
from .base import Pseudoatom

norm = scipy.linalg.get_blas_funcs("nrm2")


class Residue(MerChainable):
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

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Residue constructor.

        Arguments:
        pdb_residue -- BioPython Bio.PDB.Residue instance based on which the
        Residue is being created.
        structure_obj -- the Structure instance which the Residue belongs to.

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
        MerChainable.__init__(self, structure_obj, ind, name, chain, atoms)
        self.calculate_cbx()

    def finalize(self):
        """Method called by structures to calculate and set attributes that
        need structural
        information to be calculated.
        """
        super(Residue, self).finalize()
        self.calculate_backbone_average()

    def calculate_backbone_average(self):
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

        self.pseudoatoms["backbone_average"] = Pseudoatom(numpy_vec=(average_ca / cnt))

    @property
    def angles(self):
        """Property that returns torsion angles (in order: psi and phi) of
        residue."""
        try:
            return self.dynamic_properties["angles"]
        except KeyError:
            self.calculate_angles()
            return self.dynamic_properties["angles"]

    def calculate_angles(self):
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

        self.dynamic_properties["angles"] = (ang_psi, ang_phi)

    @property
    def ca(self):
        """Property that returns current residue alpha carbon Atom object."""
        return self.atoms["CA"]

    def calculate_cbx(self):
        """Adds Pseudoatom containing coordinates of the point that lies 1A
        farther from carbon alpha, than does carbon beta; or carbon alpha
        coordinates for GLY.
        """
        if self.get_config("legacy_cbx_calculation"):
            self.calculate_cbx_legacy()
            return
        if self.name == "GLY":
            n_2_ca = self.atoms["N"] - self.atoms["CA"]
            c_2_ca = self.atoms["C"] - self.atoms["CA"]
            average_ca_cb_distance = 1.53
            cbx = (n_2_ca + c_2_ca).get_unit_vector() * (average_ca_cb_distance + 1)
        else:
            try:
                ca = self.atoms["CA"].vector
                cb = self.atoms["CB"].vector
            except KeyError:
                raise IncompleteParticle(
                    "Mer lacks CA or CB, cannot calculate residue's cbx."
                )
            vec = cb - ca
            nrm = norm(vec)
            vec = vec * ((nrm + 1) / nrm)
            cbx = ca + vec

        self.pseudoatoms["cbx"] = Pseudoatom(numpy_vec=cbx, name="cbx")

    def calculate_cbx_legacy(self):
        """Creates pydesc.geometry.Coord instance containing coordinates of
        cbx calculated in legacy mode and assigns it to residue cbx property.

        Needs numpy to proceed. Uses Kabsch algorithm to superpose patternal
        set of C, CA, N and CBX, eeven for GLY.
        """
        pattern = [
            [1.26462, -0.673997, -3.024425],
            [0, 0, -2.5],
            [0, 0, 0],
            [-1.23670, -0.656232, -3.010602],
        ]
        # positions of atoms/points C, C alfa, C beta extended by 1 A and N,
        # respectively
        try:
            coords = [self.atoms[i] for i in ("C", "CA", "CB", "N")]
        except KeyError:
            coords = (
                self.atoms["C"],
                self.atoms["CA"],
                (self.atoms["CA"] - self.atoms["C"])
                + (self.atoms["CA"] - self.atoms["N"]),
                self.atoms["N"],
            )
        bb_coords = [coord_obj.get_coord() for coord_obj in coords]
        #

        # ===============
        # foregin code starts here
        # ===============
        # all subsequent comments untile notice were made by author
        # assertions replaced by ValueErrors
        # pylint:disable=invalid-name, no-member

        # check for consistency
        if len(bb_coords) != len(pattern):
            raise ValueError("Wrong lenght of backbone: mer %s" % str(self))
        L = len(bb_coords)

        # must alway center the two proteins to avoid
        # affine transformations.  Center the two proteins
        # to their selections.
        COM1 = numpy.sum(bb_coords, axis=0) / float(L)
        COM2 = numpy.sum(pattern, axis=0) / float(L)
        bb_coords = bb_coords - COM1
        pattern = pattern - COM2

        # This beautiful step provides the answer. V and Wt are the orthonormal
        # bases that when multiplied by each other give us the rotation
        # matrix, U. S, (Sigma, from SVD) provides us with the error!  Isn't
        # SVD great!
        V, S, Wt = numpy.linalg.svd(numpy.dot(numpy.transpose(pattern), bb_coords))

        # we already have our solution, in the aaults from SVD.
        # we just need to check for reflections and then produce
        # the rotation.  V and Wt are orthonormal, so their det's
        # are +/-1.
        reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        # U is simply V*Wt
        U = numpy.dot(V, Wt)

        # rotate and translate the molecule
        pattern = numpy.dot((pattern), U) + COM1
        pattern = pattern.tolist()

        # pylint: enable=invalid-name, no-member
        # =============
        # end of foreign code
        # =============
        self.pseudoatoms["cbx"] = Pseudoatom(*pattern[2], name="cbx")


class Nucleotide(MerChainable):  # TODO: Improve ConfigManager access

    """Representation of a nucleotide."""

    def __init__(self, structure_obj, ind, name, chain, atoms):
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
        MerChainable.__init__(self, structure_obj, ind, name, chain, atoms)

        rats = self.get_config("ring_atoms")

        def flag(name, atom):
            if name in rats:
                atom.ring_flag = True
                return True
            atom.ring_flag = False

        self.ring_atoms = {
            name: atom for name, atom in list(self.atoms.items()) if flag(name, atom)
        }

        self.calculate_ring_center()
        self.calculate_proximate_ring_center()
        self.ring_plane = None
        self.calculate_ring_plane()
        self.calculate_nx()
        self.ion_neighbours = []

    def calculate_ring_center(self):
        """Adds pseudoatom representing base ring center."""
        try:
            vec = (self.ring_atoms["N1"].vector + self.ring_atoms["C4"].vector) * 0.5
        except KeyError:
            raise IncompleteParticle("Lacking N1 or C4, unable to create Nucleotide.")
        self.pseudoatoms["ring_center"] = Pseudoatom(numpy_vec=vec, name="ring_center")

    def calculate_ring_plane(self):
        """Adds pydesc.geometry.Plane object representing base to current
        nucleotide pseudoatom dictionary."""
        at1, at2, at3 = (
            self.ring_atoms["C2"],
            self.ring_atoms["C4"],
            self.ring_atoms["C6"],
        )
        self.ring_plane = pydesc.geometry.Plane.build(
            at1, at2, at3
        )  # pylint:disable=attribute-defined-outside-init
        # current method is called by init

    def calculate_proximate_ring_center(self):
        """Adds pseudoatom representing center of the base ring being closer to
        glycosidic bond."""
        try:
            vec = numpy.array([0.0, 0.0, 0.0])
            for at in ("C4", "C5", "N7", "C8", "N9"):
                vec += self.atoms[at].vector
            vec /= 5.0
            self.pseudoatoms["prc"] = Pseudoatom(numpy_vec=vec, name="prc")
        except KeyError:
            pass

    @property
    def prc(self):
        """Get ring center of base ring closest to sugar."""
        try:
            return self.pseudoatoms["prc"]
        except KeyError:
            return self.pseudoatoms["ring_center"]

    def calculate_nx(self):
        """Adds pseudoatom representing extended by 1.4A vector along
        glycosidic bond."""
        at1 = self.atoms["C1'"]
        try:
            at2 = self.N9
        except AttributeError:
            at2 = self.N1
        vec = (at2 - at1).vector
        nrm = norm(vec)
        nvec = vec * ((nrm + 1.4) / nrm)

        nx = at1.vector + nvec
        self.pseudoatoms["nx"] = Pseudoatom(numpy_vec=nx, name="nx")


# pylint:enable=no-self-use


class Ion(MerOther):
    """Representation of an ion ligand."""

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Ion constructor.

        Sets basic attributes.

        Arguments:
        pdb_residue -- Bio.PDB.Residue instance representing ion.
        structure_obj -- instance of parental PyDesc structure.
        """
        super(Ion, self).__init__(structure_obj, ind, name, chain, atoms)
        if len(self.atoms) != 1:
            raise WrongMerType(
                "Failed to create Ion, given BioPython residue consists of "
                "to many atoms."
            )

    def get_radius(self):
        """Return ion radius."""
        name = max(self.atoms)
        try:
            return self.get_config("radii")[name]
        except KeyError:
            warn(NoConfiguration("No radius for %s ions." % name))
            return 2.5


class Ligand(MerOther):
    """Representation of any ligand except ions."""

    def __init__(self, structure_obj, ind, name, chain, atoms):
        """Ligand constructor.

        Sets basic attributes.

        Arguments:
        pdb_residue -- Bio.PDB.Residue instance representing ligands other
        than ions.
        structure_obj -- instance of parental PyDesc structure.
        """
        super(Ligand, self).__init__(structure_obj, ind, name, chain, atoms)


Mer.reset_config_cache()
