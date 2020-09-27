"""Alignment loader classes."""
import csv
import re
from abc import ABCMeta
from abc import abstractmethod

import numpy

from pydesc.alignment.base import JoinedPairAlignments
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
        self._data = {}
        self._metadata = {}
        self._structure_labels = None
        self._read = False

    @property
    def data(self):
        if not self._read:
            self._read_file()
        return self._data

    @property
    def metadata(self):
        if not self._read:
            self._read_file()
        return self._metadata

    @property
    def structure_labels(self):
        if self._structure_labels is None:
            self._read_file()
        return self._structure_labels

    @abstractmethod
    def _read_file(self):
        self._read = True

    def read_metadata(self):
        metadata = dict(self.metadata)
        metadata["labels"] = self.structure_labels
        return metadata

    @abstractmethod
    def load_partial_alignment(self, structures_map):
        pass

    def load_alignment(self, structures):
        if len(structures) != len(self.structure_labels):
            raise ValueError(
                "Number of given structures must match number of "
                "structure labels in alignment file. Consider using "
                "'load_partial_alignment' method instead."
            )
        structure_map = dict(zip(self.structure_labels, structures))
        return self.load_partial_alignment(structure_map)


class CSVLoader(AbstractLoader):
    def __init__(self, path, delimiter="\t"):
        super().__init__(path)
        self.delimiter = delimiter

    def _read_file(self):
        with open(self.path) as csv_file:
            reader = csv.reader(csv_file, delimiter=self.delimiter)
            structures_labels = next(reader)
            rows = [numpy.array(i) for i in reader if i]
        self._structure_labels = structures_labels
        self._data["rows"] = rows
        super()._read_file()

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

    def load_partial_alignment(self, structures_map):
        converters = [
            structures_map[label].converter for label in self.structure_labels
        ]
        length = len(self.data["rows"])
        array = numpy.empty((length, len(structures_map)))
        array_indices = [
            i
            for i, label in enumerate(self.structure_labels)
            if label in structures_map
        ]
        for i, row in enumerate(self.data["rows"]):
            row = row[array_indices]
            row_pdb_ids = [self._parse_pdb_id(id_str) for id_str in row]
            iterator = zip(row_pdb_ids, converters)
            row_inds = [
                self._get_ind(converter, pdb_id) for pdb_id, converter in iterator
            ]
            array[i] = row_inds
        alignment_class = get_column_alignment_class(array)
        return alignment_class(structures_map, array)


class PALLoader(AbstractLoader):
    def _read_file(self):
        with open(self.path) as file:
            n_structures = int(file.readline())
            labels = [file.readline().strip() for _ in range(n_structures)]
            self._structure_labels = labels
            while True:
                header = file.readline()
                if not header:
                    break
                if not header.startswith(">"):
                    # not a header
                    continue
                n_lines = int(file.readline().replace("@", ""))
                self._data[header] = [file.readline() for _ in range(n_lines)]
        super()._read_file()

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
    def _fit_chain_by_ind(converter, pdb_id):
        try:
            (matching_id,) = [pid for pid in converter.ind2pdb if pid[1:] == pdb_id[1:]]
        except ValueError:
            icode = pdb_id.icode or ""
            mer_id = f"{pdb_id.ind}{icode}"
            raise ValueError(f"Mer {mer_id} matches none or more than one chain.")
        return matching_id

    def _fill_chain_if_none(self, pdb_id, new_chain_name, converter):
        if pdb_id.chain is None:
            if new_chain_name is None:
                return self._fit_chain_by_ind(converter, pdb_id)
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
            id1s = self._fill_chain_if_none(id1s, stc1chain, stc1.converter)
            id1e = self._fill_chain_if_none(id1e, stc1chain, stc1.converter)
            id2s = self._fill_chain_if_none(id2s, stc2chain, stc2.converter)
            id2e = self._fill_chain_if_none(id2e, stc2chain, stc2.converter)
            new_ids.append((id1s, id1e, id2s, id2e))
        return new_ids

    @staticmethod
    def _unwrap_range(structure, id1, id2):
        ind1 = structure.converter.get_ind(id1)
        ind2 = structure.converter.get_ind(id2)
        inds = [mer.ind for mer in structure[ind1:ind2]]
        return inds

    def _unwrap_ranges(self, ranges_ids, structures):
        stc1, stc2 = structures
        inds1 = []
        inds2 = []
        for id1s, id1e, id2s, id2e in ranges_ids:
            inds1.extend(self._unwrap_range(stc1, id1s, id1e))
            inds2.extend(self._unwrap_range(stc2, id2s, id2e))
        return inds1, inds2

    def load_partial_alignment(self, structures_map):
        pair_alignments = []
        for header, ranges in self.data.items():
            try:
                structures_subset = self._get_structures(header, structures_map)
            except KeyError:
                continue
            ranges_ids = [self._parse_ranges(line) for line in ranges]
            ranges_ids = self._fill_chains(ranges_ids, structures_subset)
            inds1, inds2 = self._unwrap_ranges(ranges_ids, structures_subset)
            inds_rows = numpy.array(tuple(zip(inds1, inds2)), dtype=numpy.uint32)
            pair_alignment = PairAlignment(structures_subset, inds_rows)
            pair_alignments.append(pair_alignment)

        if len(pair_alignments) == 1:
            return pair_alignments[0]
        return JoinedPairAlignments(pair_alignments)


class FASTALoader(AbstractLoader):
    def __init__(self, path, miss_match_characters=(".", "-")):
        super().__init__(path)
        self.mmc = miss_match_characters
        self._metadata["#"] = {}
        self._metadata["ranges"] = {}

    @property
    def ranges(self):
        return self.metadata["ranges"]

    def _read_file(self):
        structure_labels = []
        with open(self.path) as fasta_file:
            current_label = None
            while True:
                line = fasta_file.readline()
                if not line:
                    break  # eof
                elif line.startswith("#"):
                    continue  # comment
                elif not line.strip():
                    continue  # blanks
                elif line.startswith(">"):
                    label, comment, ranges = self._parse_label(line)
                    structure_labels.append(label)
                    current_label = label
                    self._metadata["ranges"][label] = ranges
                    self._metadata["#"][label] = comment
                    continue
                sequence = self._data.get(current_label, "")
                self._data[current_label] = sequence + line.strip()

        self._structure_labels = structure_labels
        super()._read_file()

    @staticmethod
    def _parse_label(label_line):
        match = re.match(r">(\w*)[\ \t]*([\w\ ]*)(\[.*\])?", label_line)
        label = match.group(1)
        comment = match.group(2)
        ranges = match.group(3) or "[]"
        return label, comment, ranges

    def load_partial_alignment(self, structures_map):
        mer_iterators = {}
        for label, structure in structures_map.items():
            ranges = self._parse_ranges(self.ranges[label])
            mer_iterator = self._iter_inds(structure, ranges)
            mer_iterators[label] = mer_iterator
        lengths = {len(self.data[label]) for label in structures_map}
        if len(lengths) != 1:
            raise ValueError(
                "Sequences of given structures in fasta alignment file " "are uneven."
            )
        length = max(lengths)
        array = numpy.empty((length, len(structures_map)))
        structures = []
        for i, label in enumerate(self.structure_labels):
            if label not in structures_map:
                continue
            sequence = self.data[label]
            mer_i = mer_iterators[label]
            column = [self._get_ind(char, mer_i) for char in sequence]
            array[:, i] = numpy.array(column)
            structures.append(structures_map[label])
        array = self._remove_single_structure_rows(array)
        alignment_class = get_column_alignment_class(array)
        return alignment_class(structures, array)

    @staticmethod
    def _parse_ranges(ranges_str):
        range_pattern = r"\w+:\d+[^\d\]\[]?-\d+[^\d\]\[]?"
        ranges = re.findall(range_pattern, ranges_str)
        range_pattern = r"(\w+):(\d+)(\D?)-(\d+)(\D?)"
        results = []
        for range_str in ranges:
            match = re.match(range_pattern, range_str)
            chain = match.group(1)
            id1 = match.group(2)
            icode1 = match.group(3) or None
            id2 = match.group(4)
            icode2 = match.group(5) or None
            pdb_id1 = PDBid((chain, id1, icode1))
            pdb_id2 = PDBid((chain, id2, icode2))
            results.append((pdb_id1, pdb_id2))
        return results

    @staticmethod
    def _iter_inds(structure, ranges):
        if not ranges:
            first = structure.converter.ind2pdb[0]
            last = structure.converter.ind2pdb[-1]
            ranges = ((first, last),)
        converter = structure.converter
        for start_id, end_id in ranges:
            start = converter.get_ind(start_id)
            end = converter.get_ind(end_id)
            for mer in structure[start:end]:
                yield mer.ind

    def _get_ind(self, char, iterator):
        if char in self.mmc:
            return numpy.nan
        if re.match("[a-z]", char):
            _ = next(iterator)
            return numpy.nan
        return next(iterator)

    @staticmethod
    def _remove_single_structure_rows(array):
        _, n_structures = array.shape
        n_nans = numpy.count_nonzero(numpy.isnan(array), axis=1)
        array = array[n_nans < (n_structures - 1)]
        return array
