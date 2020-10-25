import csv
from abc import ABC
from abc import abstractmethod

from pydesc.alignment.base import DASH


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
