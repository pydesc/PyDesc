import os.path
import pytest
import numpy
from unittest.mock import MagicMock

from scipy.sparse import dok_matrix

from pydesc.contacts.geometrical import PointsDistanceCriterion
from pydesc.contacts.base import ContactsAlternative, ContactsConjunction, \
    ContactsExclusiveDisjunction, NotCriterion
from pydesc.structure import StructureLoader


@pytest.fixture
def structure(structures_dir):
    sl = StructureLoader()
    path_str = os.path.join(structures_dir, "rna_only", "1KIS.pdb")
    stc, = sl.load_structures(path=path_str)
    return stc


@pytest.fixture
def mocked_criteria():
    crit1 = MagicMock()
    crit2 = MagicMock()
    crit3 = MagicMock()
    for crit in crit1, crit2, crit3:
        crit.calculate_contacts.return_value = dok_matrix((3, 3), dtype=numpy.uint8)
    crit1.calculate_contacts.return_value[:, 0] = 2
    crit1.calculate_contacts.return_value[:, 1] = 1
    crit2.calculate_contacts.return_value[0, :] = 2
    crit2.calculate_contacts.return_value[1, :] = 1
    crit3.calculate_contacts.return_value.setdiag(2)
    crit3.calculate_contacts.return_value[0, 1] = 1
    crit3.calculate_contacts.return_value[1, 2] = 1

    return crit1, crit2, crit3


def test_point_distance(structure):
    mer1 = structure[18]  # B:19
    mer2 = structure[31]  # B:32
    mer3 = structure[4]  # A:5
    # B:19.rc is 10.5 or so far from B:32.rc
    # B:19.rc to A:5 is way more

    crt = PointsDistanceCriterion("rc", 12, 0)

    res = crt.calculate_contacts(structure)

    assert res[mer1.ind, mer2.ind] == 2
    assert res[mer1.ind, mer3.ind] == 0


def test_alternative(mocked_criteria):
    structure = MagicMock()
    crit = ContactsAlternative(*mocked_criteria)

    res = crit.calculate_contacts(structure)

    expected_result = numpy.array(
        [2, 2, 2,
         2, 2, 1,
         2, 1, 2]
    ).reshape((3, 3))

    assert (res.toarray() == expected_result).all()


def test_conjunction(mocked_criteria):
    structure = MagicMock()
    crit = ContactsConjunction(*mocked_criteria)

    res = crit.calculate_contacts(structure)

    expected_result = numpy.array(
        [2, 1, 0,
         0, 1, 0,
         0, 0, 0]
    ).reshape((3, 3))

    assert (res.toarray() == expected_result).all()


def test_exclusive_disjunction(mocked_criteria):
    structure = MagicMock()
    crit = ContactsExclusiveDisjunction(*mocked_criteria)

    res = crit.calculate_contacts(structure)

    expected_result = numpy.array(
        [0, 1, 2,
         1, 1, 1,
         2, 1, 2]
    ).reshape((3, 3))

    assert (res.toarray() == expected_result).all()


def test_not(mocked_criteria):
    structure = MagicMock()
    crit = mocked_criteria[0]
    not_crit = NotCriterion(crit)

    res = not_crit.calculate_contacts(structure)

    tmp_res = crit.calculate_contacts(structure)
    expected_res = dok_matrix((3, 3))
    expected_res[tmp_res == 2] = 0
    expected_res[tmp_res == 1] = 1
    expected_res[tmp_res == 0] = 2

    assert (res.toarray() == expected_res.toarray()).all()
