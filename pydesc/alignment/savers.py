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
"""Classes writing alignments to files."""

import csv
from abc import ABC
from abc import abstractmethod
from collections import OrderedDict
from itertools import groupby
from itertools import zip_longest

import numpy

from pydesc.alignment.base import DASH


class AbstractSaver(ABC):
    """Alignment saver interface."""

    @abstractmethod
    def save(self, stream, alignment, *, names=None):
        """Write given alignment to given file-like object.

        Args:
            stream(file-like): object with write method.
            alignment: Pair- or MultipleAlignment.
            names(dict): map structure->label(string) telling saver how to substitute
                structures with names. Must be passed as keyword argument.

        """
        pass


class CSVSaver(AbstractSaver):
    """Class saving alignments as CSV (or similar) files.

    Args:
        delimiter(str): tab by default; string to be put between columns.
        dash(str): string to be put as empty cell ("-" by default).

    """

    def __init__(self, delimiter="\t", dash="-", sequence_dict=None):
        self.delimiter = delimiter
        self.dash = dash
        self.seq_dct = sequence_dict

    def save(self, stream, alignment, *, names=None):
        if names is None:
            names = {}
        structures = alignment.structures
        writer = csv.writer(stream, delimiter=self.delimiter)
        names = [names.get(structure, structure.name) for structure in structures]
        writer.writerow(names)
        for row in alignment.iter_rows():
            text_row = self._format_row(row, structures)
            writer.writerow(text_row)

    def _format_row(self, row, structures):
        text_row = []
        for ind, structure in zip(row, structures):
            if ind is DASH:
                text = self.dash
            else:
                pdb_id = structure.converter.get_pdb_id(ind)
                letter_code = self._get_letter(structure[ind])
                text = self._format_id(pdb_id, letter_code)
            text_row.append(text)
        return text_row

    def _get_letter(self, mer):
        if self.seq_dct is None:
            return mer.seq
        return self.seq_dct[mer.name]

    @staticmethod
    def _format_id(pdb_id, letter):
        icode = pdb_id.icode or ""
        ind = f"{pdb_id.ind}{icode}"
        return f"{pdb_id.chain}:{ind}:{letter}"


class FASTASaver(AbstractSaver):
    """Class saving alignments as FASTA files.

    Args:
        dash(str): string to be put as empty cell ("-" by default).
        wrap_at(int): max length for single line of sequence.

    """

    def __init__(self, dash="-", wrap_at=None, sequence_dict=None):
        self.dash = dash
        self.wrap_at = wrap_at
        self.seq_dct = sequence_dict

    def save(self, stream, alignment, *, names=None):
        sequences = {}
        all_segments = {}
        structures = alignment.structures
        single_mer_rows = numpy.count_nonzero(alignment.inds == DASH, axis=1) == 1
        for structure, column in zip(structures, alignment.iter_columns()):
            inds = [ind for ind in column if ind is not DASH]
            segments = {}
            for k, group in groupby(inds, key=lambda ind: structure[ind].chain):
                segments.update(self._get_segments(tuple(group)))
            coded_mers = [
                self._get_mer_code(structure, ind, lower)
                for ind, lower in zip(column, single_mer_rows)
            ]
            sequence = "".join(coded_mers)
            sequences[structure] = sequence
            all_segments[structure] = segments

        for structure in structures:
            formatted_segments = []
            for pair in sorted(all_segments[structure].items()):
                segment = self._format_pdb_segment(structure, pair)
                formatted_segments.append(segment)
            segments = ", ".join(formatted_segments)
            header = f">{structure.name} [{segments}]\n"
            stream.write(header)
            self._write_sequence(stream, sequences[structure])

    @staticmethod
    def _get_segments(inds):
        segments = {}
        itr = iter(inds)
        try:
            starting_ind = next(itr)
        except StopIteration:
            return segments
        for ind1, ind2 in zip(iter(inds), itr):
            if ind1 + 1 == ind2:
                continue
            segments[starting_ind] = ind1
            starting_ind = ind2
        segments[starting_ind] = inds[-1]
        return segments

    @staticmethod
    def _format_pdb_segment(structure, pair):
        start, end = pair
        converter = structure.converter
        pdb_start = converter.get_pdb_id(start)
        pdb_end = converter.get_pdb_id(end)
        start_string = pdb_start.format(chain=True)
        end_string = pdb_end.format(chain=False)
        string = f"{start_string}-{end_string}"
        return string

    def _get_mer_code(self, structure, ind, lower_case):
        if ind is DASH:
            return self.dash
        letter = self._get_letter(structure[ind])
        if lower_case:
            return letter.lower()
        return letter.upper()

    def _get_letter(self, mer):
        if self.seq_dct is None:
            return mer.seq
        return self.seq_dct[mer.name]

    def _write_sequence(self, stream, sequence):
        if self.wrap_at is None:
            return stream.write(f"{sequence}\n")
        iterators = [iter(sequence)] * self.wrap_at
        for chunk in zip_longest(*iterators, fillvalue=" "):
            chunk = "".join(chunk)
            stream.write(f"{chunk}\n")


class PALSaver(AbstractSaver):
    """Class saving alignments as PAL(2.0) files."""

    def __init__(self):
        self.version = "2.0"

    def save(self, stream, alignment, *, names=None):
        names = self._prepare_names(alignment, names)
        chains_edges = self._get_chain_edges(alignment.structures)
        self._write_global_header(stream, alignment, names)
        joined_pairs = alignment.to_joined_pairs()
        for pair_alignment in joined_pairs.pair_alignments:
            self._write_pair_header(stream, pair_alignment, names)
            segments_map = self._prepare_segments(pair_alignment, chains_edges)
            self._write_segments(stream, pair_alignment, segments_map)
            stream.write("\n")

    @staticmethod
    def _prepare_names(alignment, names):
        all_names = {structure: structure.name for structure in alignment.structures}
        if names is not None:
            all_names.update(names)
        return all_names

    @staticmethod
    def _get_chain_edges(structures):
        chain_map = {}
        for structure in structures:
            for chain in structure.chains:
                end_ind = chain[-1].ind
                chain_map.setdefault(structure, []).append(end_ind)
        return chain_map

    def _write_global_header(self, stream, alignment, names):
        stream.write(f"v:{self.version}\n")
        stream.write(f"{len(alignment.structures)}\n")
        for structure in alignment.structures:
            stream.write(f"{names[structure]}\n")
        stream.write("\n")

    @staticmethod
    def _write_pair_header(stream, pair_alignment, names):
        structure1, structure2 = pair_alignment.structures
        name1 = names[structure1]
        name2 = names[structure2]
        stream.write(f">{name1} {name2}\n")

    @staticmethod
    def _prepare_segments(pair_alignment, chains_edges):
        left_structure, right_structure = pair_alignment.structures
        segments_map = OrderedDict()
        staring_index = 0
        length, _ = pair_alignment.inds.shape
        if length == 0:
            return segments_map
        for row_index in range(length - 1):
            nth_row, mth_row = pair_alignment.inds[[row_index, row_index + 1]]
            if numpy.all(nth_row + 1 == mth_row):
                right_is_edge = nth_row[0] in chains_edges[left_structure]
                left_is_edge = nth_row[1] in chains_edges[right_structure]
                if (not right_is_edge) and (not left_is_edge):
                    continue
            segments_map[staring_index] = row_index
            staring_index = row_index + 1
        segments_map[staring_index] = length - 1
        return segments_map

    def _write_segments(self, stream, pair_alignment, segments_map):
        left_structure, right_structure = pair_alignment.structures
        n_segments = len(segments_map)
        stream.write(f"@{n_segments}\n")
        for star_row, end_row in segments_map.items():
            left_start, right_start = pair_alignment.inds[star_row]
            left_end, right_end = pair_alignment.inds[end_row]
            left_segment = self._format_segment(left_structure, left_start, left_end)
            right_segment = self._format_segment(
                right_structure, right_start, right_end
            )
            line = f"{left_segment} <--> {right_segment}  (0,1)\n"
            stream.write(line)

    @staticmethod
    def _format_segment(structure, start_ind, end_ind, chain=True):
        converter = structure.converter
        start_pdb = converter.get_pdb_id(start_ind).format(chain)
        end_pdb = converter.get_pdb_id(end_ind).format(chain)
        return f"{start_pdb} -- {end_pdb}"
