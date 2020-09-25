import os.path

import pytest

from pydesc.config import ConfigManager
from pydesc.chemistry.full_atom import MonoatomicIon
from pydesc.chemistry.full_atom import Nucleotide
from pydesc.chemistry.full_atom import Residue

ConfigManager.warnings.class_filters.UnknownParticleName = "ignore"
ConfigManager.warnings.class_filters.Info = "ignore"

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_STRUCTURES_DIR = os.path.join(BASE_DIR, "data", "test_structures")
TEST_TRAJECTORIES_DIR = os.path.join(BASE_DIR, "data", "test_trajectories")

CIF_FILES_WITH_TYPE = []

FILE_TYPES = []
PDB_FILES_WITH_TYPE = []
PDB_FILES_WITH_TYPE_SHORT = []
PDB_FILES_DICT = {}
PDB_FILES_DICT_SHORT = {}

for type_ in os.listdir(TEST_STRUCTURES_DIR):
    FILE_TYPES.append(type_)

    files_list = os.listdir(os.path.join(TEST_STRUCTURES_DIR, type_))

    PDB_FILES_WITH_TYPE_SHORT.append(os.path.join(type_, files_list[-1]))
    PDB_FILES_DICT[type_] = []
    PDB_FILES_DICT_SHORT[type_] = [files_list[-1]]

    for structure_file in files_list:
        if structure_file.endswith(".pdb"):
            PDB_FILES_WITH_TYPE.append(os.path.join(type_, structure_file))
            PDB_FILES_DICT[type_].append(structure_file)
        elif structure_file.endswith(".cif"):
            CIF_FILES_WITH_TYPE.append((type_, structure_file))

PDB_FILES_WITH_PURE_TYPE = [
    file_name for file_name in PDB_FILES_WITH_TYPE if "only" in file_name
]
PDB_FILES_WITH_PURE_TYPE_SHORT = [
    file_name for file_name in PDB_FILES_WITH_TYPE_SHORT if "only" in file_name
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

DIR_DICT = {v: k for k, v in list(PURE_TYPES_2_MERS_DICT.items())}


@pytest.fixture(scope="session")
def structures_dir():
    return TEST_STRUCTURES_DIR


@pytest.fixture
def cmaps_dir():
    return os.path.join(BASE_DIR, "data", "validated_cmaps")


@pytest.fixture
def alignments_dir():
    return os.path.join(BASE_DIR, "data", "test_alignments")


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
    scope="session",
    params=[os.path.join(TEST_STRUCTURES_DIR, i) for i in PDB_FILES_WITH_TYPE_SHORT],
)
def structure_file_w_type_short(request):
    return request.param


@pytest.fixture(scope="session")
def mixed_structures_path():
    return os.path.join(TEST_STRUCTURES_DIR, "mixed")


@pytest.fixture(scope="session")
def ligands_structures_path():
    return os.path.join(TEST_STRUCTURES_DIR, "ligands")


def pytest_addoption(parser):
    parser.addoption(
        "--all-structures",
        action="store_true",
        help="run on all structure",
        default=False,
    )


FILES_WITH_TYPE = "structure_file_w_type"
PURE_FILES_WITH_TYPE = "structure_file_w_pure_type"
PROTEINS = "protein_file"
RNA = "rna_file"
DNA = "dna_file"
NUCLEOTIDE = "nuclei_file"


def pytest_generate_tests(metafunc):
    def set_parameter_if_all(parameter_name, if_true, if_false):
        if parameter_name in metafunc.fixturenames:
            if metafunc.config.getoption("--all-structures"):
                all_files = if_true
            else:
                all_files = if_false
            metafunc.parametrize(parameter_name, all_files)

    set_parameter_if_all(
        FILES_WITH_TYPE,
        [os.path.join(TEST_STRUCTURES_DIR, i) for i in PDB_FILES_WITH_TYPE],
        [os.path.join(TEST_STRUCTURES_DIR, i) for i in PDB_FILES_WITH_TYPE_SHORT],
    )
    set_parameter_if_all(
        PURE_FILES_WITH_TYPE,
        [os.path.join(TEST_STRUCTURES_DIR, i) for i in PDB_FILES_WITH_PURE_TYPE],
        [os.path.join(TEST_STRUCTURES_DIR, i) for i in PDB_FILES_WITH_PURE_TYPE_SHORT],
    )
    set_parameter_if_all(
        PROTEINS,
        [
            os.path.join(TEST_STRUCTURES_DIR, "prots_only", i)
            for i in PDB_FILES_DICT["prots_only"]
        ],
        [
            os.path.join(TEST_STRUCTURES_DIR, "prots_only", i)
            for i in PDB_FILES_DICT_SHORT["prots_only"]
        ],
    )
    set_parameter_if_all(
        RNA,
        [
            os.path.join(TEST_STRUCTURES_DIR, "rna_only", i)
            for i in PDB_FILES_DICT["rna_only"]
        ],
        [
            os.path.join(TEST_STRUCTURES_DIR, "rna_only", i)
            for i in PDB_FILES_DICT_SHORT["rna_only"]
        ],
    )
    set_parameter_if_all(
        DNA,
        [
            os.path.join(TEST_STRUCTURES_DIR, "dna_only", i)
            for i in PDB_FILES_DICT["dna_only"]
        ],
        [
            os.path.join(TEST_STRUCTURES_DIR, "dna_only", i)
            for i in PDB_FILES_DICT_SHORT["dna_only"]
        ],
    )
    set_parameter_if_all(
        NUCLEOTIDE,
        [i for i in PDB_FILES_WITH_PURE_TYPE if "dna" in i or "rna" in i],
        [i for i in PDB_FILES_WITH_PURE_TYPE_SHORT if "dna" in i or "rna" in i],
    )
