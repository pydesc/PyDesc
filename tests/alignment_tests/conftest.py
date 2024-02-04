import pytest
from common import get_arange_array
from common import get_n_mocked_structures
from common import get_trivial_array

from pydesc.alignment.base import Alignment


@pytest.fixture
def mocked_structures2():
    return get_n_mocked_structures(2)


@pytest.fixture
def mocked_structures3():
    return get_n_mocked_structures(3)


@pytest.fixture
def triple_alignment(mocked_structures3):
    arr = get_trivial_array(10, 3)
    alignment = Alignment(mocked_structures3, arr)
    return alignment


@pytest.fixture
def triple_nontrivial_alignment(mocked_structures3):
    arr = get_arange_array(10, 3)
    return Alignment(mocked_structures3, arr)


@pytest.fixture
def pairwise_nontrivial_alignment(mocked_structures2):
    arr = get_arange_array(10, 2)
    return Alignment(mocked_structures2, arr)
