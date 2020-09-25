"""Alignment loader classes."""
import re
import csv
from abc import ABCMeta
from abc import abstractmethod

import numpy

from pydesc.alignment.base import MultipleColumnsAlignment
from pydesc.alignment.base import PairAlignment
from pydesc.numberconverter import PDBid


def get_column_alignment_class(array):
    _, n_structures = array.shape
    if n_structures == 2:
        return PairAlignment
    elif n_structures > 2:
        return MultipleColumnsAlignment
    raise ValueError("Not enough columns to create alignment.")


class AbstractLoader(metaclass=ABCMeta):

    def __init__(self, path):
        self.path = path
        self.data = {}
        self.metadata = {}

    @abstractmethod
    def read_metadata(self):
        pass

    @abstractmethod
    def create_alignment(self, structures):
        pass


class CSVLoader(AbstractLoader):

    def __init__(self, path, delimiter="\t"):
        super().__init__(path)
        self.delimiter = delimiter

    def read_metadata(self):
        if 'structures' not in self.data:
            self._read_file()
        return self.metadata

    def _read_file(self):
        with open(self.path) as csv_file:
            reader = csv.reader(csv_file, delimiter=self.delimiter)
            structures_labels = next(reader)
            rows = [i for i in reader if i]
        self.metadata['labels'] = structures_labels
        self.data['rows'] = rows

    @staticmethod
    def _parse_pdb_id(id_str):
        match = re.match("([^:]+):([0-9]*)([^0-9,:]?)", id_str)
        try:
            chain = match.group(1)
        except AttributeError:
            return None
        no = int(match.group(2))
        i_code = match.group(3) or None
        return PDBid((chain, no, i_code))

    @staticmethod
    def _get_ind(converter, pdb_id):
        if pdb_id is None:
            return None
        return converter.get_ind(pdb_id)

    def create_alignment(self, structures):
        if 'rows' not in self.data:
            self._read_file()
        converters = [structure.converter for structure in structures]
        length = len(self.data['rows'])
        array = numpy.empty((length, len(structures)))
        for i, row in enumerate(self.data['rows']):
            row_pdb_ids = [self._parse_pdb_id(id_str) for id_str in row]
            iterator = zip(row_pdb_ids, converters)
            row_inds = [self._get_ind(converter, pdb_id) for pdb_id, converter in
                        iterator]
            array[i] = row_inds
        alignment_class = get_column_alignment_class(array)
        return alignment_class(structures, array)
