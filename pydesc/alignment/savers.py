import csv
from abc import ABC
from abc import abstractmethod
from collections import OrderedDict
from itertools import groupby
from itertools import zip_longest

import numpy

from pydesc.alignment.base import DASH


def get_segments(inds):
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


def format_pdb_segment(structure, pair):
    start, end = pair
    converter = structure.converter
    pdb_start = converter.get_pdb_id(start)
    pdb_end = converter.get_pdb_id(end)
    start_string = pdb_start.format(chain=True)
    end_string = pdb_end.format(chain=False)
    string = f"{start_string}-{end_string}"
    return string


class AbstractSaver(ABC):
    @abstractmethod
    def save(self, stream, alignment, *, names=None):
        pass


class CSVSaver(AbstractSaver):
    def __init__(self, delimiter="\t", dash="-"):
        self.delimiter = delimiter
        self.dash = dash

    def save(self, stream, alignment, *, names=None):
        if names is None:
            names = {}
        structures = alignment.structures
        writer = csv.writer(stream, delimiter=self.delimiter)
        names = [names.get(structure, structure.name) for structure in structures]
        writer.writerow(names)
        for row in alignment.iter_rows():
            text_row = self.format_row(row, structures)
            writer.writerow(text_row)

    def format_row(self, row, structures):
        text_row = []
        for ind, structure in zip(row, structures):
            if ind is DASH:
                text = self.dash
            else:
                pdb_id = structure.converter.get_pdb_id(ind)
                letter_code = structure[ind].seq
                text = self.format_id(pdb_id, letter_code)
            text_row.append(text)
        return text_row

    @staticmethod
    def format_id(pdb_id, letter):
        icode = pdb_id.icode or ""
        ind = f"{pdb_id.ind}{icode}"
        return f"{pdb_id.chain}:{ind}:{letter}"


class FASTASaver(AbstractSaver):
    def __init__(self, dash="-", wrap_at=None):
        self.dash = dash
        self.wrap_at = wrap_at

    def save(self, stream, alignment, *, names=None):
        sequences = {}
        all_segments = {}
        structures = alignment.structures
        single_mer_rows = numpy.count_nonzero(alignment.inds == DASH, axis=1) == 1
        for structure, column in zip(structures, alignment.iter_columns()):
            inds = [ind for ind in column if ind is not DASH]
            segments = {}
            for k, group in groupby(inds, key=lambda ind: structure[ind].chain):
                segments.update(get_segments(tuple(group)))
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
                segment = format_pdb_segment(structure, pair)
                formatted_segments.append(segment)
            segments = ", ".join(formatted_segments)
            header = f">{structure.name} [{segments}]\n"
            stream.write(header)
            self._write_sequence(stream, sequences[structure])

    def _get_mer_code(self, structure, ind, lower_case):
        if ind is DASH:
            return self.dash
        letter = structure[ind].seq
        if lower_case:
            return letter.lower()
        return letter.upper()

    def _write_sequence(self, stream, sequence):
        if self.wrap_at is None:
            return stream.write(f"{sequence}\n")
        iterators = [iter(sequence)] * self.wrap_at
        for chunk in zip_longest(*iterators, fillvalue=" "):
            chunk = "".join(chunk)
            stream.write(f"{chunk}\n")


class PALSaver(AbstractSaver):
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
