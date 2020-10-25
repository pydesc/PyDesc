import csv
import numpy
from itertools import zip_longest, groupby
from abc import ABC
from abc import abstractmethod

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
