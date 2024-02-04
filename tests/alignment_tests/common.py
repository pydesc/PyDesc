from unittest.mock import MagicMock

import numpy


def get_n_mocked_structures(n, *, start=0):
    lst = [MagicMock(name=f"stc{n}") for n in range(start, start + n)]
    for n, mock in enumerate(lst):  # has to be done this way because of MagicMock
        mock.name = f"stc{n}"
    return lst


def get_trivial_array(l, n):
    arr, _ = numpy.indices((l, n))
    return arr.astype(object)


def get_arange_array(l, n):
    return numpy.arange(n * l).reshape(n, l).T.astype(object)
