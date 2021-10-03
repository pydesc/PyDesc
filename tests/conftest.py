import os.path
from itertools import chain
from pathlib import Path

import pytest

from pydesc.chemistry.full_atom import MonoatomicIon
from pydesc.chemistry.full_atom import Nucleotide
from pydesc.chemistry.full_atom import Residue
from pydesc.config import ConfigManager

ConfigManager.warnings.class_filters.UnknownParticleName = "ignore"
ConfigManager.warnings.class_filters.Info = "ignore"

BASE_DIR = Path(__file__).parent
TEST_STRUCTURES_DIR = BASE_DIR / "data" / "test_structures"
TEST_TRAJECTORIES_DIR = BASE_DIR / "data" / "test_trajectories"

CIF_FILES = {}
PDB_FILES = {}

for path in TEST_STRUCTURES_DIR.glob("**/*[.pdb|.cif]"):
    if path.is_dir():
        continue
    type_name = path.parent.stem
    if path.suffix == ".pdb":
        PDB_FILES.setdefault(type_name, []).append(path)
    elif path.suffix == ".cif":
        CIF_FILES.setdefault(type_name, []).append(path)

ONE_OF_EACH_KIND = [PDB_FILES[key][-1] for key in PDB_FILES] + [
    CIF_FILES[key][-1] for key in CIF_FILES
]

PURE_TYPES_2_MERS_DICT = {
    "dna_only": Nucleotide,
    "dna_only_nmr": Nucleotide,
    "rna_only": Nucleotide,
    "rna_only_nmr": Nucleotide,
    "prots_only": Residue,
    "prots_only_nmr": Residue,
    "PorCA_only": MonoatomicIon,
}

NMR_TYPES = ("prots_only_nmr", "rna_only_nmr", "dna_only_nmr")

DIR_DICT = {v: k for k, v in list(PURE_TYPES_2_MERS_DICT.items())}


def get_keys_iterator(*dcts, keys=None):
    if keys is None:
        keys = list(set([key for dct in dcts for key in dct]))
    return list(chain(*[dct.get(key, []) for dct in dcts for key in keys]))


def get_last_entries(*dcts, keys=None):
    if keys is None:
        keys = [key for dct in dcts for key in dct]
    return [dct[key][-1] for dct in dcts for key in keys if key in dct]


@pytest.fixture(scope="session")
def structures_dir():
    return TEST_STRUCTURES_DIR


@pytest.fixture(scope="session")
def cmaps_dir():
    return BASE_DIR / "data" / "validated_cmaps"


@pytest.fixture(scope="session")
def alignments_dir():
    return BASE_DIR / "data" / "test_alignments"


@pytest.fixture(scope="session")
def trajectories_path():
    return os.path.join(TEST_TRAJECTORIES_DIR, "{type}", "{traj}.{type}")


@pytest.fixture(scope="session")
def topologies_path():
    return os.path.join(TEST_TRAJECTORIES_DIR, "topologies", "{topo}.pdb")


@pytest.fixture
def pure_types_2_mers():
    return PURE_TYPES_2_MERS_DICT


@pytest.fixture(
    scope="session", params=ONE_OF_EACH_KIND,
)
def any_structure_file(request):
    return request.param


@pytest.fixture(
    params=[PDB_FILES[kind][-1] for kind in NMR_TYPES])
def nmr_structure_of_each_kind(request):
    return request.param

@pytest.fixture(scope="session")
def mixed_structures_path():
    return TEST_STRUCTURES_DIR / "mixed"


@pytest.fixture(scope="session")
def ligands_structures_path():
    return TEST_STRUCTURES_DIR / "ligands"


def pytest_addoption(parser):
    parser.addoption(
        "--all-structures",
        action="store_true",
        help="run on all structure",
        default=False,
    )


PDB_FILE = "pdb_file"
CIF_FILE = "cif_file"
PURE_FILE = "pure_file"
PROTEINS = "protein_file"
RNA = "rna_file"
DNA = "dna_file"
NUCLEOTIDE = "nuclei_file"
TRACE = "trace_file"


def pytest_generate_tests(metafunc):
    def set_parameter_if_all(parameter_name, if_true, if_false):
        if parameter_name in metafunc.fixturenames:
            if metafunc.config.getoption("--all-structures"):
                all_files = if_true
            else:
                all_files = if_false
            metafunc.parametrize(parameter_name, all_files)

    set_parameter_if_all(
        CIF_FILE, get_keys_iterator(CIF_FILES), get_last_entries(CIF_FILES)
    )

    set_parameter_if_all(
        PDB_FILE, get_keys_iterator(PDB_FILES), get_last_entries(PDB_FILES)
    )

    pure_keys = [key for key in PDB_FILES if "only" in key]
    set_parameter_if_all(
        PURE_FILE,
        get_keys_iterator(PDB_FILES, CIF_FILES, keys=pure_keys),
        get_last_entries(PDB_FILES, CIF_FILES, keys=pure_keys),
    )

    set_parameter_if_all(
        PROTEINS,
        get_keys_iterator(PDB_FILES, CIF_FILES, keys=["prots_only"]),
        get_last_entries(PDB_FILES, CIF_FILES, keys=["prots_only"]),
    )

    rna_long = get_keys_iterator(PDB_FILES, CIF_FILES, keys=["rna_only"])
    rna_short = get_last_entries(PDB_FILES, CIF_FILES, keys=["rna_only"])
    set_parameter_if_all(RNA, rna_long, rna_short)

    dna_long = get_keys_iterator(PDB_FILES, CIF_FILES, keys=["dna_only"])
    dna_short = get_last_entries(PDB_FILES, CIF_FILES, keys=["dna_only"])
    set_parameter_if_all(DNA, dna_long, dna_short)

    set_parameter_if_all(NUCLEOTIDE, rna_long + dna_long, rna_short + dna_short)

    set_parameter_if_all(
        TRACE,
        get_keys_iterator(PDB_FILES, CIF_FILES, keys="PorCA_only"),
        get_last_entries(PDB_FILES, CIF_FILES, keys="PorCA_only"),
    )
