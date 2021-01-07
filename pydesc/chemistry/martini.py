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

"""Classes storing mers and ligands data in Martini representation."""

from pydesc.chemistry.base import Mer
from pydesc.config import ConfigManager

ConfigManager.chemistry.new_branch("martiniresidue")
ConfigManager.chemistry.martiniresidue.set_default("backbone_atoms", ("BB",))
ConfigManager.chemistry.martiniresidue.set_default("bb_bond_threshold", 5.0)


class MartiniResidue(Mer):
    """Residue in representation suitable for Martini force field."""

    @property
    def last_sc(self):
        """Side-chain atom of greatest index."""
        self_len = len(self.atoms)
        try:
            return self.atoms[f"SC{self_len - 1}"]
        except KeyError:
            return self.atoms["BB"]
