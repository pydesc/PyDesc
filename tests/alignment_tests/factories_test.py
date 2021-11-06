from unittest.mock import Mock

import pytest

from pydesc.alignment.base import DASH
from pydesc.alignment.factory import AlignmentFactory


@pytest.fixture
def structure():
    return Mock()


def test_from_list_of_inds(structure):
    factory = AlignmentFactory()
    structures = [structure] * 3
    inds_lists = ((0, 1, 2, 4), (2, 3, 4, 5), (1, None, 42, 24))
    alignment = factory.create_from_list_of_inds(structures, inds_lists)
    assert alignment.inds.shape == (3, 4)
    assert alignment.inds[1, 3] is DASH
