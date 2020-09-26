"""Alignment loader classes."""
import csv
import re
from abc import ABCMeta
from abc import abstractmethod

import numpy

from pydesc.alignment.base import MultipleColumnsAlignment
from pydesc.alignment.base import PairAlignment, JoinedPairAlignments
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
        self.structure_labels = None

    @abstractmethod
    def _read_file(self):
        pass

    def read_metadata(self):
        if self.structure_labels is None:
            self._read_file()
        metadata = dict(self.metadata)
        metadata['labels'] = self.structure_labels
        return metadata

    @abstractmethod
    def create_alignment(self, structures):
        pass


class CSVLoader(AbstractLoader):

    def __init__(self, path, delimiter="\t"):
        super().__init__(path)
        self.delimiter = delimiter

    def _read_file(self):
        with open(self.path) as csv_file:
            reader = csv.reader(csv_file, delimiter=self.delimiter)
            structures_labels = next(reader)
            rows = [i for i in reader if i]
        self.structure_labels = structures_labels
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


class PALLoader(AbstractLoader):

    def _read_file(self):
        with open(self.path) as file:
            n_structures = int(file.readline())
            labels = [file.readline().strip() for _ in range(n_structures)]
            self.structure_labels = labels
            while True:
                header = file.readline()
                if not header:
                    break
                if not header.startswith(">"):
                    # not a header
                    continue
                n_lines = int(file.readline().replace("@", ""))
                self.data[header] = [file.readline() for _ in range(n_lines)]

    @staticmethod
    def _get_structures(header, label_map):
        match = re.match(">(.*) (.*)\n", header)
        label1, label2 = match.group(1), match.group(2)
        structure1 = label_map[label1]
        structure2 = label_map[label2]
        return structure1, structure2

    @staticmethod
    def _parse_pdb_id(pdb_str):
        match = re.match("([^0-9]*)([0-9]*)(.?)", pdb_str.strip())
        chain = match.group(1) or None
        no = int(match.group(2))
        i_code = match.group(3) or None
        return PDBid((chain, no, i_code))

    def _parse_ranges(self, ranges):
        id_pattern = "[^0-9]*[0-9]*.?"
        range_pattern = f"({id_pattern}) -- ({id_pattern})"
        match = re.match(f"{range_pattern} <--> {range_pattern}", ranges)
        ids = [match.group(i) for i in range(1, 5)]
        ids = [self._parse_pdb_id(id_str) for id_str in ids]
        return ids

    @staticmethod
    def _fill_chain_if_none(pdb_id, new_chain_name):
        if pdb_id.chain is None:
            pdb_id = PDBid((new_chain_name, pdb_id.ind, pdb_id.icode))
            return pdb_id
        return pdb_id

    def _fill_chains(self, ranges_ids, structures):
        stc1chain = None
        stc2chain = None
        stc1, stc2 = structures
        if len(stc1.chains) == 1:
            stc1chain = stc1.chains[0].chain_name
        if len(stc2.chains) == 1:
            stc2chain = stc2.chains[0].chain_name
        new_ids = []
        for id1s, id1e, id2s, id2e in ranges_ids:
            id1s = self._fill_chain_if_none(id1s, stc1chain)
            id1e = self._fill_chain_if_none(id1e, stc1chain)
            id2s = self._fill_chain_if_none(id2s, stc2chain)
            id2e = self._fill_chain_if_none(id2e, stc2chain)
            new_ids.append((id1s, id1e, id2s, id2e))
        return new_ids

    @staticmethod
    def _unwrap_range(structure, id1, id2):
        ind1 = structure.converter.get_ind(id1)
        ind2 = structure.converter.get_ind(id2)
        inds = [mer.ind for mer in structure[ind1: ind2]]
        return inds

    def _unwrap_ranges(self, ranges_ids, structures):
        stc1, stc2 = structures
        inds1 = []
        inds2 = []
        for id1s, id1e, id2s, id2e in ranges_ids:
            inds1.extend(self._unwrap_range(stc1, id1s, id1e))
            inds2.extend(self._unwrap_range(stc2, id2s, id2e))
        return inds1, inds2

    def create_alignment(self, structures):
        if not self.data:
            self._read_file()
        label2structure = dict(zip(self.structure_labels, structures))
        pair_alignments = []
        for header, ranges in self.data.items():
            structures_subset = self._get_structures(header, label2structure)
            ranges_ids = [self._parse_ranges(line) for line in ranges]
            ranges_ids = self._fill_chains(ranges_ids, structures_subset)
            inds1, inds2 = self._unwrap_ranges(ranges_ids, structures_subset)
            inds_rows = numpy.array(tuple(zip(inds1, inds2)))
            pair_alignment = PairAlignment(structures_subset, inds_rows)
            pair_alignments.append(pair_alignment)

        if len(pair_alignments) == 1:
            return pair_alignments[0]
        return JoinedPairAlignments(pair_alignments)
