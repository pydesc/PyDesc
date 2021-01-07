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

"""Classes storing mers and ligands data in backbone-trace representation."""

import scipy.linalg

from pydesc.chemistry.base import Mer
from pydesc.chemistry.base import Pseudoatom
from pydesc.chemistry.base import register_pseudoatom
from pydesc.config import ConfigManager

norm = scipy.linalg.get_blas_funcs("nrm2")

ConfigManager.chemistry.new_branch("catrace")
ConfigManager.chemistry.catrace.set_default("bb_bond_threshold", 5.0)
ConfigManager.chemistry.catrace.set_default("backbone_atoms", ("CA",))
ConfigManager.chemistry.catrace.set_default("indicators", ("CA", "mpp"))
ConfigManager.chemistry.catrace.set_default("mpp_length", 2.53)
ConfigManager.chemistry.new_branch("ptrace")
ConfigManager.chemistry.ptrace.set_default("bb_bond_threshold", 7.5)
ConfigManager.chemistry.ptrace.set_default("backbone_atoms", ("P",))
ConfigManager.chemistry.ptrace.set_default("indicators", ("P",))


class CATrace(Mer):
    """Carbon alpha trace representation of residue."""

    @register_pseudoatom
    def mpp(self):
        """Most protruding (from backbone) point."""
        try:
            cbx_coord = self._calculate_mpp()
        except AttributeError:
            cbx_coord = self._calculate_terminus_mpp()
        return Pseudoatom(numpy_vec=cbx_coord, name="cbx")

    def _calculate_mpp(self):
        """Calculate mpp coord for CA in the middle of chain."""
        av_value = self.get_config("mpp_length")
        self_ca = self.atoms["CA"]
        previous_ca = self.prev_mer.atoms["CA"]
        next_ca = self.next_mer.atoms["CA"]
        prev_to_this = self_ca - previous_ca
        next_to_this = self_ca - next_ca
        direction = (prev_to_this + next_to_this).get_unit_vector()
        cbx_coord = self_ca.vector + direction * av_value
        return cbx_coord

    def _calculate_terminus_mpp(self):
        """Calculate mpp coord for terminal CAs."""
        av_value = self.get_config("mpp_length")
        adjacent_mer = self.next_mer or self.prev_mer
        self_ca = self.atoms["CA"]
        if adjacent_mer is None:
            return self_ca.vector.copy() + [av_value, 0, 0]
        direction = (adjacent_mer.atoms["CA"] - self_ca).get_unit_vector()
        cbx_coord = self_ca.vector + direction * -av_value
        return cbx_coord


class PTrace(Mer):
    """Phosphorus trace representation of nucleotides."""

    pass
