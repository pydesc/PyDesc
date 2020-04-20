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

"""
Classes that represents mers present in representations of (sub)structures.

created: 11.07.2013 - 31.07.2013, Tymoteusz 'hert' Oleniecki
"""

from pydesc.config import ConfigManager

# pylint: disable=no-member
ConfigManager.new_branch("mers")
ConfigManager.mers.set_default("monomer_acceptable_distance", 2.0)
ConfigManager.mers.set_default("solvent", ["HOH"])
ConfigManager.mers.new_branch("nucleotide")
ConfigManager.mers.new_branch("residue")
ConfigManager.mers.new_branch("monomerchainable")
ConfigManager.mers.new_branch("ion")
ConfigManager.mers.new_branch("ligand")
ConfigManager.mers.set_default("backbone_atoms", ())
ConfigManager.mers.monomerchainable.set_default("check_distances", False)
ConfigManager.mers.residue.set_default(
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
ConfigManager.mers.residue.set_default(
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
ConfigManager.mers.residue.set_default("backbone_atoms", ("N", "CA", "C"))
ConfigManager.mers.residue.set_default("check_distances", False)
ConfigManager.mers.residue.set_default(
    "crucial_atom_distances", (("C", "CA", 1.35, 1.71), ("CA", "N", 1.35, 1.75))
)
ConfigManager.mers.residue.set_default("indicators", ("CA", "cbx"))
ConfigManager.mers.residue.set_default("adjusted_segment_length", 18.0)
ConfigManager.mers.nucleotide.set_default(
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
ConfigManager.mers.nucleotide.set_default(
    "backbone_atoms", ("P", "O5'", "C5'", "C4'", "C3'", "O3'")
)
ConfigManager.mers.nucleotide.set_default(
    "ring_atoms", ("N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9")
)
ConfigManager.mers.nucleotide.set_default("check_distances", False)
ConfigManager.mers.nucleotide.set_default(
    "crucial_atom_distances",
    (
        ("P", "O5'", 1.54, 1.66),
        ("O5'", "C5'", 1.34, 1.54),
        ("C5'", "C4'", 1.44, 1.56),
        ("C4'", "C3'", 1.46, 1.58),
        ("C3'", "O3'", 1.37, 1.49),
    ),
)
ConfigManager.mers.nucleotide.set_default("indicators", ("C3'", "P", "ring_center"))
ConfigManager.mers.set_default("moving_average", 3)
ConfigManager.mers.ion.set_default("indicators", ("rc",))
ConfigManager.mers.ion.set_default(
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

ConfigManager.mers.ligand.set_default("indicators", ("rc",))
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

# pylint: enable=no-member
