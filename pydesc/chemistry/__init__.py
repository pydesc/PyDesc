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

"""Classes representing building blocks of structures: mer and compounds."""

from pydesc.config import ConfigManager

# pylint: disable=no-member
ConfigManager.new_branch("chemistry")
ConfigManager.chemistry.set_default("backbone_atoms", ())
ConfigManager.chemistry.set_default("bb_bond_threshold", 2.0)
ConfigManager.chemistry.set_default("solvent", ["HOH"])
ConfigManager.chemistry.set_default("moving_average", 3)
ConfigManager.chemistry.new_branch("mer")
ConfigManager.chemistry.mer.set_default("check_distances", False)

ConfigManager.new_branch("structure_mon")
ConfigManager.structure_mon.set_default(
    "simple_secondary_structure_code",
    {
        "H": "H",
        "B": "E",
        "E": "E",
        "G": "H",
        "I": "H",
        "T": "C",
        "S": "C",
        "-": "C",
        "=": "=",
    },
)
# TODO: should that be residue setting?

# pylint: enable=no-member
