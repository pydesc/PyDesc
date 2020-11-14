# Copyright 2020 Tymoteusz 'vdhert' Oleniecki
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
"""Alignment loader classes."""

import csv
import re
from abc import ABC
from abc import abstractmethod

import numpy

from pydesc.alignment.base import DASH
from pydesc.alignment.base import JoinedPairAlignments
from pydesc.alignment.base import MultipleAlignment
from pydesc.alignment.base import PairAlignment
from pydesc.alignment.base import drop_single_mer_rows
from pydesc.numberconverter import PDBid


def get_column_alignment_class(array):
    """Given array of aligned inds, return class that best suits it or raise
    ValueError for invalid array."""
    _, n_structures = array.shape
    if n_structures == 2:
        return PairAlignment
    elif n_structures > 2:
        return MultipleAlignment
    raise ValueError("Not enough columns to create alignment.")


class AbstractLoader(ABC):
    """Common loader abstract superclass."""

    def __init__(self, path):
        self.path = path
        self._data = {}
        self._metadata = {}
        self._structure_labels = None
        self._read = False

    @property
    def data(self):
        """Return alignment data read from alignment file."""
        if not self._read:
            self._read_file()
        return self._data

    @property
    def metadata(self):
        """Return metadata read from alignment file."""
        if not self._read:
            self._read_file()
        return self._metadata

    @property
    def structure_labels(self):
        """Return structure labels read from alignment file."""
        if self._structure_labels is None:
            self._read_file()
        return self._structure_labels

    @abstractmethod
    def _read_file(self):
        self._read = True

    def read_metadata(self):
        """Read meta data from alignment file, especially structure labels.

        Metadata is returned as dict with guaranteed key "labels". Other data depends of
        alignment file type.

        """
        metadata = dict(self.metadata)
        metadata["labels"] = self.structure_labels
        return metadata

    @abstractmethod
    def load_alignment_mapping(self, structures_map):
        """Load alignment mapping labels to structures as in given map.

        Passing map lacking some labels present in alignment file will result in
        loading partial alignment.

        Args:
            structures_map(dict): map with labels(strings) as keys and structures as
            values.

        """
        pass

    def _prepare_map_for_structures(self, structures):
        structure_map = dict(zip(self.structure_labels, structures))
        return structure_map

    def load_alignment(self, structures):
        """Load alignment assuming labels correspond with given structures (order
        matters).

        If given list of structures is shorter, only as many columns will be loaded
        as many structures were given.

        """
        structure_map = self._prepare_map_for_structures(structures)
        alignment = self.load_alignment_mapping(structure_map)
        return alignment


class CSVLoader(AbstractLoader):
    """Loader able to read alignment in CSV (or similar) format."""

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
        if match is None:
            return None
        chain = match.group(1)
        no = int(match.group(2))
        i_code = match.group(3) or None
        return PDBid((chain, no, i_code))

    @staticmethod
    def _get_ind(converter, pdb_id):
        if pdb_id is None:
            return DASH
        return converter.get_ind(pdb_id)

    def load_alignment_mapping(self, structures_map):
        valid_labels = [i for i in self.structure_labels if i in structures_map]
        converters = [structures_map[label].converter for label in valid_labels]
        length = len(self.data["rows"])
        array = numpy.empty((length, len(structures_map)), dtype=object)
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
    def __init__(self, path):
        super(PALLoader, self).__init__(path)
        self.version = None
        self.id_parser = None

    @staticmethod
    def _get_id_parser(version):
        dct = {
            "1.0": PALIdParserV1(),
            "2.0": PALIdParserV2(),
        }
        return dct[version]

    def _read_file(self):
        with open(self.path) as file:
            first_line = file.readline()
            version_match = re.match("v:([0-9]+.[0-9]+)", first_line)
            if version_match is not None:
                self.version = version_match.group(1)
                n_structures = int(file.readline())
            else:
                n_structures = int(first_line)
                self.version = "1.0"
            self.id_parser = self._get_id_parser(self.version)
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

    def load_alignment_mapping(self, structures_map):
        joined_pairs = self.load_joined_pairs_mapping(structures_map)
        alignment = joined_pairs.to_columns()
        return alignment

    def load_joined_pairs(self, structures):
        """Load alignment as container of pair alignments.

        That is native approach for PAL format. Might be useful if single file stores
        disjointed or inconsistent pair alignment, not multiple alignment separated
        to particular pair alignments.

        """
        structure_map = self._prepare_map_for_structures(structures)
        alignment = self.load_joined_pairs_mapping(structure_map)
        return alignment

    def load_joined_pairs_mapping(self, structures_map):
        """Method corresponding to load_alignment_mapping, but returning container
        for pair alignments."""
        pair_alignments = []
        for header, ranges in self.data.items():
            try:
                structures_subset = self._get_structures(header, structures_map)
            except KeyError:
                continue
            ranges_ids = [self._parse_ranges(line) for line in ranges]
            ranges_ids = self._fill_chains(ranges_ids, structures_subset)
            inds1, inds2 = self._unwrap_ranges(ranges_ids, structures_subset)
            inds_rows = numpy.array(tuple(zip(inds1, inds2)), dtype=object)
            pair_alignment = PairAlignment(structures_subset, inds_rows)
            pair_alignments.append(pair_alignment)

        if len(pair_alignments) == 1:
            return pair_alignments[0]
        return JoinedPairAlignments(pair_alignments)

    @staticmethod
    def _get_structures(header, label_map):
        match = re.match(">(.*) (.*)\n", header)
        label1, label2 = match.group(1), match.group(2)
        structure1 = label_map[label1]
        structure2 = label_map[label2]
        return structure1, structure2

    def _parse_ranges(self, ranges):
        range_pattern = f"(\\S*) -- (\\S*)"
        match = re.match(f"{range_pattern} <--> {range_pattern}", ranges)
        ids = [match.group(i) for i in range(1, 5)]
        ids = [self.id_parser.parse_pdb_id(id_str) for id_str in ids]
        return ids

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

    def _fill_chain_if_none(self, pdb_id, new_chain_name, converter):
        if pdb_id.chain is None:
            if new_chain_name is None:
                return self._fit_chain_by_ind(converter, pdb_id)
            pdb_id = PDBid((new_chain_name, pdb_id.ind, pdb_id.icode))
            return pdb_id
        return pdb_id

    @staticmethod
    def _fit_chain_by_ind(converter, pdb_id):
        try:
            (matching_id,) = [pid for pid in converter.ind2pdb if pid[1:] == pdb_id[1:]]
        except ValueError:
            icode = pdb_id.icode or ""
            mer_id = f"{pdb_id.ind}{icode}"
            raise ValueError(f"Mer {mer_id} matches none or more than one chain.")
        return matching_id

    def _unwrap_ranges(self, ranges_ids, structures):
        stc1, stc2 = structures
        inds1 = []
        inds2 = []
        for id1s, id1e, id2s, id2e in ranges_ids:
            inds1.extend(self._unwrap_range(stc1, id1s, id1e))
            inds2.extend(self._unwrap_range(stc2, id2s, id2e))
        return inds1, inds2

    @staticmethod
    def _unwrap_range(structure, id1, id2):
        ind1 = structure.converter.get_ind(id1)
        ind2 = structure.converter.get_ind(id2)
        inds = [mer.ind for mer in structure[ind1:ind2]]
        return inds


class PALIdParser(ABC):
    """Abstract PDB id parser for PAL loader."""

    @abstractmethod
    def parse_pdb_id(self, pdb_str):
        """Parse given pdb id string."""
        pass


class PALIdParserV1(PALIdParser):
    """PDB id parser for PAL loader in version 1.0."""

    def __init__(self):
        self.pid_pattern = "([^0-9]*)([0-9]*)(.?)"

    def parse_pdb_id(self, pdb_str):
        match = re.match(self.pid_pattern, pdb_str.strip())
        chain = match.group(1) or None
        no = int(match.group(2))
        i_code = match.group(3) or None
        return PDBid((chain, no, i_code))


class PALIdParserV2(PALIdParser):
    """PDB id parser for PAL loader in version 2.0."""

    def __init__(self):
        self.long_pid_pattern = "([^:]*):([0-9]*)(.?)"
        self.short_pid_pattern = "([0-9]*)(.?)"

    def parse_pdb_id(self, pdb_str):
        if ":" in pdb_str:
            return self._parse_long(pdb_str)
        return self._parse_short(pdb_str)

    def _parse_short(self, pdb_str):
        match = re.match(self.short_pid_pattern, pdb_str.strip())
        no = int(match.group(1))
        i_code = match.group(2) or None
        return PDBid((None, no, i_code))

    def _parse_long(self, pdb_str):
        match = re.match(self.long_pid_pattern, pdb_str.strip())
        chain = match.group(1) or None
        no = int(match.group(2))
        i_code = match.group(3) or None
        return PDBid((chain, no, i_code))


class FASTALoader(AbstractLoader):
    def __init__(self, path, miss_match_characters=(".", "-")):
        super().__init__(path)
        self.mmc = miss_match_characters
        self._metadata["#"] = {}
        self._metadata["ranges"] = {}

    @property
    def ranges(self):
        """Ranges stored in alignment file."""
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
        match = re.match(r">(\w*)[ \t]*([\w ]*)(\[.*])?", label_line)
        label = match.group(1)
        comment = match.group(2)
        ranges = match.group(3) or "[]"
        return label, comment, ranges

    def load_alignment_mapping(self, structures_map):
        mer_iterators = {}
        for label, structure in structures_map.items():
            ranges = self._parse_ranges(self.ranges[label])
            mer_iterator = self._iter_inds(structure, ranges)
            mer_iterators[label] = mer_iterator
        lengths = {len(self.data[label]) for label in structures_map}
        if len(lengths) != 1:
            msg = "Sequences of given structures in fasta alignment file " "are uneven."
            raise ValueError(msg)
        length = max(lengths)
        array = numpy.empty((length, len(structures_map)), dtype=object)
        structures = []
        ind = 0
        for label in self.structure_labels:
            if label not in structures_map:
                continue
            sequence = self.data[label]
            mer_i = mer_iterators[label]
            column = [self._get_ind(char, mer_i) for char in sequence]
            array[:, ind] = numpy.array(column)
            structures.append(structures_map[label])
            ind += 1
        alignment_class = get_column_alignment_class(array)
        array = drop_single_mer_rows(array)
        alignment = alignment_class(structures, array)
        return alignment

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
            return DASH
        if re.match("[a-z]", char):
            _ = next(iterator)
            return DASH
        return next(iterator)
