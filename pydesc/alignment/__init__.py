# Copyright 2017 Tymoteusz Oleniecki
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

import pydesc.config
import pydesc.dbhandler
import pydesc.selection
import pydesc.structure
import pydesc.warnexcept

import re
import numpy
import operator
import tempfile
import contextmeta
import xml.dom.minidom

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline

from abc import ABCMeta
from itertools import chain
from StringIO import StringIO

from pydesc.cydesc.overfit import Overfit

# pylint: disable=no-member
pydesc.config.ConfigManager.new_branch("alignments")
pydesc.config.ConfigManager.alignments.set_default("multiple_alignment_mode", "strict")
# pylint: enable=no-member


def get_segments(mer_sequence, sort=True):
    """Returns list of tuples containing starting and ending segment mers.

    Argument:
    mer_sequence -- list or tuple of pydesc mers.

    One-mer segments are also returned (tuples containt the same mer twice).
    """
    if sort:
        mer_sequence = sorted(mer_sequence, key=lambda mer: mer.ind)
    segs = []
    try:
        prev = start = mer_sequence[0]
    except IndexError:
        return []
    for mer in mer_sequence[1:] + [None]:
        if mer is not None and mer.previous_monomer == prev:
            prev = mer
            continue
        segs.append((start, prev))
        start = prev = mer
    return segs


def split_to_int_and_icode(string):
    """Returns tuple of pdb integer and pdb insertion code for given pdb id.

    Argument:
    string -- loaded from alignment file pdb number of residue - containing digital int and insertion code (sign or letter).

    If there is no insertion code - None is returned as insertion code.
    """
    if string.isdigit():
        pdb_number = int(string)
        insertion_code = None
    else:
        print repr(string)
        pdb_number = int(string[:-1])
        insertion_code = string[-1]
    return pdb_number, insertion_code


def blast_sequences(structure_1, structure_2, nucleotide=False, **kwargs):  # pylint:disable=too-many-locals
    """Return PairAlignment based on sequence alignment performed by BLAST.

    Arguments:
    structure_1, structure_2 -- pydesc.AbstractStructure subclasses instances to be aligned.
    nucleotides -- boolean. By default set to False. If so, uses blastp, otherwise blastn is used.
    kwargs -- all other key-word arguments to be passed to blast.

    To run this method NCBI+ has to be installed.
    Returns tuple containing two elements: NCBI BLAST result parsed by BioPython and pydesc PairAlignment object.
    See BioPython help to read more about BLAST result.
    """
    s1mers = pydesc.selection.MonomerType(pydesc.monomer.MonomerChainable).create_structure(structure_1)
    s2mers = pydesc.selection.MonomerType(pydesc.monomer.MonomerChainable).create_structure(structure_2)
    seq1 = SeqRecord(Seq(s1mers.get_sequence()), id=str(structure_1))
    seq2 = SeqRecord(Seq(s2mers.get_sequence()), id=str(structure_2))

    blast = NcbiblastnCommandline if nucleotide else NcbiblastpCommandline

    if 0 in (len(seq1), len(seq2)):
        pydesc.warnexcept.warn(Warning('One of following structures you are trying to BLAST has no sequence -- %s, %s' % (str(structure_1), str(structure_2))))

    with tempfile.NamedTemporaryFile() as file1:
        with tempfile.NamedTemporaryFile() as file2:
            SeqIO.write(seq1, file1.name, "fasta")
            SeqIO.write(seq2, file2.name, "fasta")
            blastres = blast(query=file1.name, subject=file2.name, outfmt=5, **kwargs)()[0]
    try:
        res = NCBIXML.read(StringIO(blastres)).alignments[0].hsps[0]
    except IndexError:
        raise IndexError('No alignment for found with BLAST.')

    s1i = iter(list(s1mers)[res.query_start - 1:res.query_end])
    s2i = iter(list(s2mers)[res.sbjct_start - 1:res.sbjct_end])
    aligned_mers = []
    for match, chr1, chr2 in zip(res.match, res.query, res.sbjct):
        if match != ' ':
            aligned_mers.append((s1i.next(), s2i.next()))
            continue
        for char, item in zip((chr1, chr2), (s1i, s2i)):
            if char != '-':
                dummy = item.next()

    return (res, Alignment.build_from_list_of_mers((structure_1, structure_2), aligned_mers))

def find_sequence(sequence, structure, nucleotide=False, **kwargs):
    """Returns substructures that fit to given pattern.

    Arguments:
    sequence -- string; sequence pattern to be found in structure.
    structure -- pydesc.AbstractStructure subclasses instances to search with given pattern.
    nucleotides -- boolean. By default set to False. If so, uses blastp, otherwise blastn is used.
    kwargs -- all other key-word arguments to be passed to blast.

    To run this method NCBI+ has to be installed.
    Returns list of tuples containing two elements: substructure and NCBI BLAST result parsed by BioPython.
    See BioPython help to read more about BLAST result.
    """
    seqmers = pydesc.selection.MonomerType(pydesc.monomer.MonomerChainable).create_structure(structure)
    seq2 = SeqRecord(Seq(seqmers.get_sequence()), id=str(structure))

    blast = NcbiblastnCommandline if nucleotide else NcbiblastpCommandline

    with tempfile.NamedTemporaryFile() as file1:
        with tempfile.NamedTemporaryFile() as file2:
            SeqIO.write(SeqRecord(Seq(sequence), id='tempseq'), file1.name, "fasta")
            SeqIO.write(seq2, file2.name, "fasta")
            blastres = blast(query=file1.name, subject=file2.name, outfmt=5, **kwargs)()[0]

    try:
        res = NCBIXML.read(StringIO(blastres)).alignments[0].hsps
    except IndexError:
        raise IndexError('No alignment for %s found with BLAST.' % structure.name)

    def mk_thx(fit):
        mrs = list(seqmers)[fit.sbjct_start - 1: fit.sbjct_end]
        inds = map(len, filter(bool, fit.sbjct.split('-')))
        i1 = 0
        rngs = []
        for i2 in inds:
            i2 = i2 + i1
            rngs.append(pydesc.selection.Range(mrs[i1].pid, mrs[i2 - 1].pid).create_structure(structure))
            i1 = i2
        return rngs
    #~ mk_thx = lambda x: pydesc.selection.Range(list(seqmers)[x.sbjct_start - 1].pid, list(seqmers)[x.sbjct_end - 1].pid).create_structure(structure)

    return [(mk_thx(mtch), mtch) for mtch in res]


class AlignmentLoader(object):

    """Class responsible for loading alignments from files."""

    def __init__(self, handler=pydesc.dbhandler.MetaHandler(), reload_=False, model_no=None):
        """Alignment loader constructor.

        Arguments:
        handler -- instance of pydesc.dbhandler. By default set to MetaHandler.
        reload_ -- True or False; determines if structures are to be reloaded when same structures are aligned in differen alignments. By default set to False.
        model_no -- dictionary with structures codes as keys and their model numbers as values. Determines which models are to be loaded by structure loader.
        Initiall set to None. If so - empty dictionary is used.
        """
        self.structure_loader = pydesc.structure.StructureLoader()
        self.handler = handler
        self.reload_structures = reload_
        self.structures = {}
        if model_no is None:
            model_no = {}
        self.__models = model_no

    def load_alignment(self, file_name):
        """Returns alignment object loaded from file.

        Argument:
        file_name -- path to alignment file.

        Method is able to load csv, fasta, xml and pal files.
        """
        extension = file_name.split("/")[-1].split('.')[1]
        if extension == 'xml':
            file_data = self.load_xml(file_name)
        elif extension == 'csv':
            file_data = self.load_csv(file_name)
        elif extension == "pal":
            file_data = self.load_pal(file_name)
        else:
            file_data = self.load_fasta(file_name)
        return file_data

    # ??? koniecznie napisac, ze struktury nie sa ladowane drugi raz i ze NMR tylko pierwszy jest brany
    def load_structure(self, structure_code, path=None):
        """Returns pydesc structure objects of given pdb code.

        Argument:
        structure_code -- string; structure code. See pydesc.structure.StructureLoader for more information.
        path -- path to pdb file. Optional.

        Returns structure with given code. If structure was already loaded - new structure is not loaded, old structure is returned.
        """
        # ew co z NMR?
        # zaznaczyc tez, na czym polega rownosc struktur, bo zaladowanie
        # podobnej struktury spowoduje niemoznosc liczenia domkniecia, przejsc
        # etc.!!!
        if not self.reload_structures:
            try:
                return self.structures[structure_code]
            except KeyError:
                pass
        try:
            model = self.__models[structure_code]
        except KeyError:
            model = 0
        if path is None:
            self.structures[structure_code] = self.structure_loader.load_structure(structure_code)[model]
        else:
            self.structures[structure_code] = self.structure_loader.load_structure(structure_code, path=path)[model]
        return self.structures[structure_code]

    def load_xml(self, file_name):
        """Returns proper alignment object loaded from xml file.

        Argument:
        file_name -- file path to xml file.

        Returns PairAlignment or MultipleAlignment object, depending on number of aligned structures.
        """

        def convert(mer, structure_obj):
            """Converts mer string to pydesc ind.

            Arguments:
            mer -- string; mer name.
            structure_obj -- pydesc AbstrctStructure instance.
            """
            try:
                pdb_number, dummy_code, chain = mer.split(":")
                pdb_int, icode = split_to_int_and_icode(pdb_number.strip())
                id_tuple = (chain, pdb_int, icode)
                ind = structure_obj.converter.get_ind(id_tuple)
                return structure_obj[ind]
            except Exception:   # ??? jakiego typu bledy przy nieuliniowionych merach?
                return "-"

        alignment_doc = xml.dom.minidom.parse(file_name).childNodes[0]
        structures = [self.load_structure(member.childNodes[0].nodeValue[:4]) for member in alignment_doc.getElementsByTagName('member')]
        aligned_mers = [[mer.childNodes[0].nodeValue for mer in row.getElementsByTagName('meq')] for row in alignment_doc.getElementsByTagName('row')]
        aligned_mers = [map(convert, alignment_tuple, structures) for alignment_tuple in aligned_mers]
        return Alignment.build_from_list_of_mers(structures, aligned_mers)

    def load_csv(self, file_name):
        """Returns proper alignment object. loaded from csv file.

        Argument:
        file_name -- file path.

        Returns PairAlignment or MultipleAlignment object, depending on number of aligned structures.
        """

        with open(file_name, "r") as file_obj:
            csv_lines = file_obj.read().splitlines()
        structures = [self.load_structure(structure_obj) for structure_obj in csv_lines[0].split()]

        def get_mers((structure_index, pdb_tuple)):
            """Returns mer ind for given structure index and pdb_id-tuple."""
            if type(pdb_tuple) is str:
                return pdb_tuple
            return structures[structure_index][structures[structure_index].converter.get_ind(pdb_tuple)]

        def convert_csv_mer_to_pdb_tuple(csv_mer):
            """Converts csv mers to pdb tuples."""
            if csv_mer.count(":") != 2:
                return csv_mer
            mer_parts = csv_mer.split(":")
            pdb_int, icode = split_to_int_and_icode(mer_parts[1])
            return (mer_parts[0], pdb_int, icode)

        aligned_mers = []
        for line in csv_lines[1:]:
            line = line.split()
            alignment_tuple = map(convert_csv_mer_to_pdb_tuple, line)
            alignment_tuple = tuple(map(get_mers, enumerate(alignment_tuple)))
            aligned_mers.append(alignment_tuple)
        return Alignment.build_from_list_of_mers(structures, aligned_mers)

    def load_pal(self, file_name):  # pylint:disable=too-many-locals
        """Returns proper alignment object loaded from pal file.

        Argument:
        file_name -- file path.

        Returns PairAlignment or MultipleAlignment object, depending on number of aligned structures.
        """

        def load(stcn):
            try:
                return self.load_structure(stcn[:-1])
            except pydesc.dbhandler.InvalidID:
                return self.load_structure(stcn)

        with open(file_name) as file_:
            lines = map(str.strip, file_.readlines())

        stc_n = int(lines.pop(0))
        structures = {i.name: i for i in [load(lines.pop(0)) for i in range(stc_n)]}
        alignments = []
        while True:
            try:
                line = lines.pop(0)
            except IndexError:
                break
            try:
                (r1s, r1e), (r2s, r2e) = [map(str.strip, i.split("--")) for i in line[:line.find("(")].split("<-->")]
                try:
                    aligned_mers.extend(zip(s1[r1s: r1e], s2[r2s: r2e]))
                except KeyError:
                    # aligned_mers.extend(zip(s1[s1.name[5].upper() + r1s: s1.name[5].upper() + r1e], s2[s2.name[5].upper() + r2s: s2.name[5].upper() + r2e]))
                    aligned_mers.extend(zip(s1[ch1 + r1s: ch1 + r1e], s2[ch2 + r2s: ch2 + r2e]))
            except ValueError:
                # print line
                if line.startswith(">"):
                    try:
                        alignments.append(PairAlignment((s1, s2), aligned_mers))
                    except NameError:
                        pass

                    s1name, s2name=line[1:].split()

                    # print s1name, s2name
                    try:
                        s1 = structures[s1name]
                    except KeyError:
                        s1 = structures[s1name[:-1]]
                        ch1 = s1name[-1]

                    try:
                        s2 = structures[s2name]
                    except KeyError:
                        s2 = structures[s2name[:-1]]
                        ch2 = s2name[-1]

                    aligned_mers = []
                # elif line.strip() == '':
                    # alignments.append(PairAlignment((s1, s2), aligned_mers))

        alignments.append(PairAlignment((s1, s2), aligned_mers))

        return Alignment.build_from_pair_alignments(alignments)

    def load_fasta(self, file_name):
        """Returns proper alignment object loaded from fasta file.

        Argument:
        file_name -- file path.

        Returns PairAlignment or MultipleAlignment object, depending on number of aligned structures.
        """

        letters = pydesc.config.ConfigManager.monomer.residue.residue_additional_code.values()

        with open(file_name, "r") as file_:
            fasta_blocks = file_.read().split(">")[1:]

        headers = []
        seqs = []
        for i in fasta_blocks:
            splt = i.split('\n')
            headers.append(splt[0])
            seqs.append("".join(splt[1:]))
        structures = [self.load_structure(i.split()[0]) for i in headers]

        def cvrt(rng):
            chn=rng.split(":")[0]
            return [chn + i for i in rng.split(":")[1].split('-')]

        def mk_rng(patt, stream):
            try:
                return reduce(operator.add, [pydesc.selection.Range(*cvrt(i)) for i in re.findall(patt, stream)[0][1:-1].split(',')])
            except IndexError:
                return pydesc.selection.Everything()

        rng_pattern = re.compile('\[.:.*\]')
        iters = [iter(mk_rng(rng_pattern, hea).create_structure(stc)) for hea, stc in zip(headers, structures)]

        def drop_mer((it, le), dbg=False):
            if not le.isalpha():
                return '-'
            mer = it.next()
            if le not in letters:
                return mer  # for unusual mers we assume next one is right one
            while mer.seq != le:
                if dbg:
                    print 'skipping', mer, 'because of', le
                mer = it.next()
            return mer

        try:
            mrs = [tuple(map(drop_mer, zip(iters, col))) for col in zip(*seqs)]
        except StopIteration:
            rlen = [len(mk_rng(rng_pattern, hea).create_structure(stc)) for hea, stc in zip(headers, structures)]
            alen = [len(filter(lambda x: str.isalpha(x) and x in letters, i)) for i in seqs]
            bad_ent = filter(lambda (rl, al, stc): rl < al, zip(rlen, alen, structures))
            if bad_ent == []:
                for stc, hea, seq in zip(structures, headers, seqs):
                    it = iter(mk_rng(rng_pattern, hea).create_structure(stc))
                    try:
                        [drop_mer((it, s), True) for s in seq]
                    except StopIteration:
                        bad_ent.append(stc)
            raise ValueError('Following structures are to short to cover range given in fasta file: %s' % (','.join(map(str, bad_ent))))

        return Alignment.build_from_list_of_mers(structures, mrs) 


class Alignment(object):

    """Class that stores informations about aligned mers and provides basic operations on them."""

    # __metaclass__ = ABCMeta  # ??? nie mozna tego polaczyc z contexmeta

    @staticmethod
    def build_from_pair_alignments(list_of_pair_alignments):
        """Static method of building alignment from list of alignemnts.

        Argument:
        list_of_pair_alignments -- sequence of PairAlignment instances.

        Returns first element of given sequence if seqence length is 1. Otherwise returns MultipleAlignment instance.
        """
        if not len(list_of_pair_alignments) < 2:
            return MultipleAlignment(alignments=list_of_pair_alignments)    # pylint: disable=unexpected-keyword-arg, no-value-for-parameter
            # this arg is provided by ContextStateMeta class
            # since given parameteres steer kind of __init__ to be used - there is no point in giving 'structures' nor 'list_of_aligned_mers' parameters values
        else:
            return list_of_pair_alignments[0]

    @staticmethod
    def build_from_list_of_mers(list_of_structures, list_of_mers):
        """Static method of building alignment from given list of aligned mers.

        Arguments:
        list_of_structures -- ordered sequence of pydesc.structure.AbstractStructure instances.
        list_of_mers -- sequence of tuples containing aligned mers coming from appropriate structures from list_of_structures.

        Returns PairAlignment instance if there are two structures and two mers in each tuple in aligned_mers sequence.
        Otherwise returns MultipleAlignment instance.
        """
        if len(list_of_structures) == 2:
            return PairAlignment(list_of_structures, list_of_mers)
        return MultipleAlignment(structures=list_of_structures, list_of_aligned_mers=list_of_mers)

    def get_structures_names(self):
        """Returns a list of names of aligned structures in order given during initialization."""
        return [str(structure_obj) for structure_obj in self.structures]      # pylint: disable=no-member
        # all subclasses provides structures attr

    def get_specified_structures(self):
        return [sel.create_structure(stc) for sel, stc in zip(self.get_selections(), self.structures)]

    def get_selections(self):
        """Returns list of selections of mers used in current alignment in order refering to structures order."""
        selections = [[] for i in self.structures]      # pylint: disable=no-member
        # all subclasses provides structures attr
        for i, str_segments in enumerate(map(get_segments, [[mer for mer in mers if mer != '-'] for mers in map(set, zip(*self.aligned_mers))])):   # pylint: disable=no-member
            set_ = []
            for segment in str_segments:
                if segment[0] == segment[1]:
                    set_.append(segment[0])
                else:
                    selections[i].append(pydesc.selection.Range(*map(operator.methodcaller('get_pdb_id'), segment)))
            if len(set_) != 0:
                selections[i].append(pydesc.selection.Set(map(operator.methodcaller('get_pdb_id'), set_)))
        return [pydesc.selection.SelectionsUnion(sels) if len(sels) != 1 else sels[0] for sels in selections]

    def get_aligned_with(self, mer):
        """Returns narrowed alignment containing only mers alignmed with given mer.

        Argument:
        mer -- instance of pydesc.monomer.MonomerChainable subclass.
        """
        res = zip(*[mers for mers in zip(self.structures, *[tup for tup in self if mer in tup]) if any(mer != "-" for mer in mers[1:])])  # pylint:disable=no-member
        # res[0] is a list of filtered structures, res[1:] is a list of lists of aligned mers; structures that have no aligned mers are excluded
        if res == []:
            res = [[], []]
        return self.build_from_list_of_mers(res[0], res[1:])

    def sum_domain_with(self, alignment_obj):
        """Joins current alignment with given alignemnt and returns MultipleAlignment instance.

        Argument:
        alignment_obj -- instance of Alignment class.

        Raises AttributeError if current and given alignments do not align at least one common structure.
        Structures are considered the same when they are the same python object.
        User can provied that by loading alignments with the same AlignmentLoader with reload_ argument set on False.
        """
        #~ if not any(str(structure_obj) in alignment_obj.get_structures_names() for structure_obj in self.structures):    # pylint: disable=no-member
            # checking, if any structure occures in all alignments
            #~ raise AttributeError("Set of structures in given pair_alignments are disjoint, cannot sum pair_alignments' domains")
        pair_alignments = self.alignments + alignment_obj.alignments    # pylint: disable=no-member
        return Alignment.build_from_pair_alignments(pair_alignments)

    def save_csv(self, file_name):
        """Saves alignment in csv format file.

        file_name -- string; path to file to be created/overridden.
        """

        def mer_to_csv(mer):
            """Turns given mer into csv string that describes mer PDB id."""
            str_ = lambda x: '' if x is None else str(x)
            if mer == "-":
                return mer
            pdb_id = mer.get_pdb_id()
            return pdb_id[0] + ":" + "".join(map(str_, pdb_id[1:])) + ":" + mer.seq

        lines = ["\t".join(self.get_structures_names())]
        for aligned_mers in self:
            lines.append("\t".join([mer_to_csv(mer) for mer in aligned_mers]))
        with open(file_name + ".csv", "w") as file_:
            for line in lines:
                file_.write(line + "\n")

    def save_pal(self, file_name):
        """Saves alignment in pal format file.

        file_name -- string; path to file to be created/overridden.
        """
        parts = []
        names = []
        for pair_alignment in self.alignments:  # pylint: disable=no-member
            parts.append(pair_alignment._create_pal_string())   # pylint: disable=protected-access
            # the point of keeping this method is to access it here, but it is protected for users to know it is not important for them
            for line in parts[-1].splitlines():
                if line.startswith(">"):
                    names.extend(line[1:].split())
        names = set(names)
        with open(file_name + ".pal", "w") as file_:
            file_.write(str(len(names)) + "\n")
            for name in names:
                file_.write(name + "\n")
            for part in parts:
                file_.write(part)

    def save_fasta(self, file_name):
        """Saves alignment in fasta format file.

        file_name -- string; path to file to be created/overridden.
        """
        # pylint: disable=no-member
        mers = dict((structure_obj, []) for structure_obj in self.structures)
        for alignment_tuple in self.aligned_mers:
            for index, mer in enumerate(alignment_tuple):
                mers[self.structures[index]].append(mer)
        regions = dict((structure_obj, []) for structure_obj in self.structures)
        for structure_obj, list_of_mers in mers.items():
            list_of_mers = [mer for mer in list_of_mers if not isinstance(mer, str)]
            try:
                start = list_of_mers[0]
            except IndexError:
                continue
            for pair in zip(list_of_mers[:-1], list_of_mers[1:]):
                if pair[0].next_monomer == pair[1]:
                    continue
                end = pair[0]
                regions[structure_obj].append((start, end))
                start = pair[1]
            regions[structure_obj].append((start, list_of_mers[-1]))
        with open(file_name, "w") as file_:
            dct = {}
            for structure_obj in self.structures:
                # pylint: enable=no-member
                # structures attr is set in all subclasses of Alignment
                switch_to_fasta = lambda (start, end): start.my_chain + ":" + start.pid[1:] + "-" + end.pid[1:]
                ranges = ", ".join(map(switch_to_fasta, regions[structure_obj]))
                dct.setdefault(structure_obj, []).append(
                    ">" + str(structure_obj) + "[" + ranges + "]" + "\n")
                dct[structure_obj].append(
                    "".join([mer.seq if type(mer) != str else "." for mer in mers[structure_obj]]) + "\n")
            for k in sorted(dct):
                file_.write(dct[k][0])
                file_.write(dct[k][1])

    def save_xml(self, file_name):
        """Saves alignment in xml format file.

        file_name -- string; path to file to be created/overridden.
        """

        def convert(mer):
            """Turns given mer into xml-formated string describing mer PDB id."""
            if mer == "-":
                return "---------"
            pdb_id = mer.structure.converter.get_pdb_id(mer.ind)
            return str(pdb_id[1]) + pdb_id[2] + ":" + mer.seq + ":" + pdb_id[0]

        # pylint: disable=no-member
        docu = xml.dom.minidom.getDOMImplementation().createDocument(None, "multiple-alignment", None)
        ma_xml = docu.childNodes[0]
        ma_xml.attributes['n'] = str(len(self.structures))
        ma_xml.appendChild(docu.createElement("description"))
        members = docu.createElement("members")
        for structure_obj in self.structures:
            member = docu.createElement("member")
            member.appendChild(docu.createTextNode(str(structure_obj)))
            members.appendChild(member)
        ma_xml.appendChild(members)
        alter = docu.createElement("alternative")
        alter.attributes['id'] = "1"
        alter.attributes["n"] = str(len(self.aligned_mers))
        ma_xml.appendChild(alter)
        meqs = docu.createElement("mequivalences")
        meqs.attributes["n"] = str(len(self.aligned_mers))
        ma_xml.getElementsByTagName("alternative")[0].appendChild(meqs)
        for alignment_tuple in self.aligned_mers:
            row = docu.createElement("row")
            for mer in alignment_tuple:
                meq = docu.createElement("meq")
                meq.appendChild(docu.createTextNode(convert(mer)))
                row.appendChild(meq)
            meqs.appendChild(row)
        with open(file_name, "w") as file_:
            docu.writexml(file_, addindent="  ", newl="\n")
        # pylint: enable=no-member
        # there is no point in arguing with minidom doc about its classes attrs

    def expand(self):
        """Add tuples containing all mers of all structures (aligned with no other mers from other structures)."""
        all_mers = [set(stc) for stc in self.structures]
        aligned_mers_by_stc = zip(*self.aligned_mers)
        mers_to_add = [full - set(alg) for full, alg in zip(all_mers, aligned_mers_by_stc)]
        n_stcs = len(self.structures)
        for stc_ind, mers in enumerate(mers_to_add):
            for mer in mers:
                row = list("-" * (n_stcs - 1))
                row.insert(stc_ind, mer)
                self.aligned_mers.append(tuple(row))

    def sort(self):
        """Sort tuples in self.aligned_mers by indexes of their mers in structures."""

        class TSLst(object):

            def __init__(self, c):
                self.c = c

            def lt(self, it1, it2, no=False):
                for i, j in zip(it1, it2):
                    if "-" in (i, j):
                        continue
                    return i.ind < j.ind
                return no

            def push_fwd(self, it):
                for n, i in enumerate(self.c):
                    if self.lt(it, i):
                        self.c.insert(n, it)
                        break

        tsl = TSLst(sorted([i for i in self.aligned_mers if i[0] != '-'], key=lambda x: x[0].ind))

        dct = {}
        for tup in self.aligned_mers:
            if tup in tsl.c:
                continue
            dct.setdefault(tup.count("-"), []).append(tup)

        def show(tup):
            print " ".join(map(str, [i if type(i) is str else i.ind for i in tup]))

        for key in sorted(dct):
            for tup in dct[key]:
                tsl.push_fwd(tup)

        self.aligned_mers = tsl.c


class PairAlignment(Alignment):

    """Class that stores information about alignment of two structures."""

    __metaclass__ = ABCMeta

    def __init__(self, list_of_structures, list_of_aligned_mers):
        """PairAlignment constructor.

        Argument:
        list_of_structures -- sequence of aligned structures.
        list_of_aligned_mers -- sequence of tuples containing aligned pydesc.Monomer instances representing mers from appropriate structures in list_of_structures.
        """
        self.structures = list_of_structures
        self.aligned_mers = sorted(list_of_aligned_mers, key=lambda aligned_mers: aligned_mers[0].ind)
        self.alignments = [self]
        self._hashdict = None

    def _make_hash(self):
        """Sets _hashdict attribute - dict translating mers from one structure for aligned mers form second structure."""
        self._hashdict = dict((structure_obj, {}) for structure_obj in self.structures)
        str1, str2 = self.structures
        for mer1, mer2 in self.aligned_mers:
            self._hashdict[str1][mer1] = mer2
            self._hashdict[str2][mer2] = mer1

    def _hash(self, mer):
        """Returns mer aligned with given mer."""
        if self._hashdict is None:
            self._make_hash()
        return self._hashdict[mer.structure][mer]

    def __iter__(self):
        """Returns iterator that iterates over list of tuples containing aligned mers."""
        return iter(self.aligned_mers)

    def __len__(self):
        """Returns length of list of aligned mers."""
        return len(self.aligned_mers)

    def __repr__(self):
        filling = (str(self.structures[0]), str(self.structures[1]), len(self.aligned_mers))
        return "<PairAlignment of %s and %s: %i aligned mers>" % filling

    def _get_aligned_selections(self):
        """Returns pair of Set selections or pair of SelectionsUnion selections of aligned structures.

        If unions ale returned - they consists of other unions or of Ranges selections (Set if to short) of mers in order given in current alignemt.
        Ranges of mers aligned with two ranges in aligned structure are splited into two subsequent ranges.

        Calling method iter_recursively on returned unions strongly recommended.
        Method used to pretty printing of pair alignemnts, thus - not supported.
        """
        if len(self.aligned_mers) == 1:
            return map(pydesc.selection.Set, map(operator.methodcaller('get_pdb_id'), self.aligned_mers[0]))
        seg0 = [[self.aligned_mers[0][0]]]
        seg1 = [[self.aligned_mers[0][1]]]
        for (prev0, prev1), (curr0, curr1) in zip(self.aligned_mers[:-1], self.aligned_mers[1:]):
            if prev0.next_monomer == curr0 and prev1.next_monomer == curr1:
                continue
            seg0[-1].append(prev0)
            seg0.append([curr0])
            seg1[-1].append(prev1)
            seg1.append([curr1])
        seg0[-1].append(self.aligned_mers[-1][0])
        seg1[-1].append(self.aligned_mers[-1][1])

        def mkrange(st_en):
            """Returns range selection"""
            try:
                return pydesc.selection.Range(*map(operator.methodcaller("get_pdb_id"), st_en))
            except ValueError:
                # risen if starting and ending mer is the same
                return pydesc.selection.Set([st_en[0].get_pdb_id()])

        uni0 = pydesc.selection.SelectionsUnion(map(mkrange, seg0))
        uni1 = pydesc.selection.SelectionsUnion(map(mkrange, seg1))
        return uni0, uni1

    def __str__(self):  # pylint:disable=too-many-locals
        # we need all of them and ther is no point in spliting mthod in two
        """Prints alignment in pretty format, segment by segment."""
        string = []
        selections = self._get_aligned_selections()

        if all(isinstance(s, pydesc.selection.CombinedSelection) for s in selections):
            list_of_aligned_sels = zip(*map(operator.methodcaller('iter_recursively'), selections))
        else:
            list_of_aligned_sels = [tuple(selections)]
        max_struc = max(map(len, map(str, self.structures)))

        def get_string_ranges(selection):
            """Returns string containing PDB ids of given segment starting and ending mers."""
            try:
                return "-".join(map(str, [selection.start, selection.end]))
            except AttributeError:  # non-Range selections have no start and end attr
                return ",".join(map(str, selection.ids))

        ranges_str = [map(get_string_ranges, aligned_sel_pair) for aligned_sel_pair in list_of_aligned_sels]
        max_range = max(map(len, reduce(operator.add, ranges_str)))
        for num, (aligned_selections, ranges) in enumerate(zip(list_of_aligned_sels, ranges_str)):
            aligned_range_string = []
            for struc, sel, range_str in zip(self.structures, aligned_selections, ranges):
                name = str(struc)
                name_margin = "".join([" " for i in range(max_struc - len(name) + 4)])
                range_margin = "".join([" " for i in range(max_range - len(range_str) + 4)])
                seq = sel.create_structure(struc).get_sequence()
                seq = seq if seq != '' else '(%s)' % sel.create_structure(struc)[0].name    #???
                aligned_range_string.append(name + name_margin + range_str + range_margin + seq + "\n")
            intermission1 = "".join([" " for i in range(max_struc + max_range + 8)])
            intermission2 = "".join(["-" if i % 5 != 0 else "|" for i in range(len(aligned_range_string[-1]) - len(intermission1) - 1)])
            intermission = intermission1 + intermission2 + "\n"
            string.append("%i.\n" % (num + 1,) + intermission.join(aligned_range_string))
        return "\n".join(string)

    def _create_pal_string(self):
        """Creates string containing alignment information in pal file format."""

        def mkk(mer1, mer2):
            return (str(mer1.structure) + mer1.my_chain, str(mer2.structure) + mer2.my_chain)

        def gt_pids(mer1, mer2):
            return (mer1.pid, mer2.pid)

        #~ def add_range(smer1, smer2):
            #~ ranges[mkk(smer1, smer2)] = [(smer1.pid, smer2.pid)]

        pair0 = self.aligned_mers[0]
        ranges = {mkk(*pair0): [gt_pids(*pair0)]}
        #~ add_range(*pair0)

        for current_mers, next_mers in zip(self.aligned_mers[:-1], self.aligned_mers[1:]):
            # two subsequent pairs of aligned mers are taken: n and n + 1
            ent = mkk(*current_mers)
            if current_mers[0].is_next(next_mers[0]) and current_mers[1].is_next(next_mers[1]):
                continue
            ranges[ent].append(gt_pids(*current_mers))
            ranges.setdefault(mkk(*next_mers), []).append(gt_pids(*next_mers))
        pairm1 = self.aligned_mers[-1]
        ranges[mkk(*pairm1)].append(gt_pids(*pairm1))
        sections = []
        for aligned_chains, rng in ranges.items():
            strings_to_join = [">" + " ".join(aligned_chains) + "\n"]
            strings_to_join += ["@%i\n" % (len(rng) / 2)]
            for (m1s, m2s), (m1e, m2e) in zip(rng[::2], rng[1::2]):   # m for mer, 1 or 2 for column, s or e for start or end
                inpt = "%s -- %s <--> %s -- %s (1,0)\n" % (m1s, m1e, m2s, m2e) # ??? co za 1,0?
                strings_to_join += [inpt]
            sections.append("".join(strings_to_join))
        return "".join(sections)

    def _prepare_comparison(self, pairalignment):
        """Returns sorted lists of aligned mers from current and given alignments.

        Argument:
        pairalignment -- instance of PairAlignment.

        Raises ValueError if structures aligned in both alignments are different.
        Returned lists are sorted by default key and contain sorted tuples of aligned mers.
        """
        if sorted(self.structures) != sorted(pairalignment.structures):
            raise ValueError('Different sets of structures aligned')
        foregin_sorted_mers = sorted(map(sorted, pairalignment.aligned_mers))
        native_sorted_mers = sorted(map(sorted, self.aligned_mers))
        return sorted(foregin_sorted_mers), sorted(native_sorted_mers)

    def __eq__(self, pairalignment):
        """Returns True if compared PairAlignments align same sets of mers."""
        try:
            foregin_sorted_mers, native_sorted_mers = self._prepare_comparison(pairalignment)
        except ValueError:
            return False
        except AttributeError:
            return False
        if foregin_sorted_mers == native_sorted_mers:
            return True
        return False

    def __ge__(self, pairalignment):
        """Returns True if right compared alignment contains all mers aligned in left alignment."""
        try:
            foregin_sorted_mers, native_sorted_mers = self._prepare_comparison(pairalignment)
        except ValueError:
            return False
        if all(aligned_mers in native_sorted_mers for aligned_mers in foregin_sorted_mers):
            return True
        return False

    def __gt__(self, pairalignment):
        """Returns True if right compared alignment contains all mers aligned in left alignment, but not only them."""
        if self >= pairalignment and not self == pairalignment:
            return True
        return False

    def transit(self, alignment_obj):
        """Returns PairAlignment instance that stores information about putative alignment between two structures that are aligned with structure common for both - current and given - alignments.

        Argument:
        alignment_obj -- instance of PairAlignment.

        Given alignment has to have exactly one structure common with current alignment object.
        """
        common_structure = [structure_obj for structure_obj in self.structures if structure_obj in alignment_obj.structures][0]
        transition_dict = {}
        # next step is to add mers from common structure present in both alignments as keys
        # respective mers from both structures aligned with common structure
        # will occure as values
        for alignment_obj_ in (self, alignment_obj):
            common_index = alignment_obj_.structures.index(common_structure)
            # getting aligned mer tuple index of mers to be put in dictionary
            # as key
            for aligned_mers_ in alignment_obj_:
                aligned_mers_ = list(aligned_mers_)
                key = aligned_mers_.pop(common_index)
                try:
                    transition_dict[key] += aligned_mers_
                except KeyError:
                    transition_dict[key] = aligned_mers_
        aligned_mers = [mers for mers in transition_dict.values() if None not in mers]
        list_of_structures = [mer.structure for mer in aligned_mers[0]]
        return PairAlignment(list_of_structures, aligned_mers)

    def is_consistent(self, pair_alignment):
        """Checks if current and given PairAlignments stores noncontradictory relation between sets of aminoacids.

        Argument:
        pair_alignment -- instance of PairAlignment to be compared.
        """
        strcs1 = map(operator.attrgetter('derived_from'), self.structures)
        strcs2 = map(operator.attrgetter('derived_from'), pair_alignment.structures)
        selfdict = dict(self)
        selfdict_rev = dict(zip(*reversed(zip(*self))))
        if strcs1 == strcs2:
            otherdict = dict(pair_alignment)
            otherdict_rev = dict(zip(*reversed(zip(*pair_alignment))))
        elif strcs1 == list(reversed(strcs2)):
            otherdict = dict(zip(*reversed(zip(*pair_alignment))))
            otherdict_rev = dict(pair_alignment)
        else:
            # when alignment aligns different sets of structures - they are consistent
            return True

        def comp(d1, d2, key):
            """Checks if two values are different in dwo dicts, but returns False if keys are not in dict."""
            try:
                return d1[key] != d2[key]
            except KeyError:
                return False

        def check_dict(alg1, d1, d1r, d2, d2r):
            """Returns True if all keys in two dicts returns the same value (if present), otherwise - returns False."""
            for m1, m2 in alg1:
                if comp(d1, d2, m1) or comp(d1r, d2r, m2):
                    return False
            return True

        return check_dict(self, selfdict, selfdict_rev, otherdict, otherdict_rev) and check_dict(pair_alignment, otherdict, otherdict_rev, selfdict, selfdict_rev)

    def calculate_tension(self, cmaps=None):
        """Compares contact maps of aligned structures.

        Argument:
        cmaps -- optional, list of contact maps to be used insted of ones present in structures.

        Method creats all possible pydesc.structure.Contact objects basing on contacts present in
        at least one of two compared contact maps. Average RMSD for such contacts is returned.
        """
        conts = set([])
        if cmaps:
            cm1, cm2 = cmaps
        else:
            try:
                cm1, cm2 = map(operator.attrgetter("contact_map"), self.structures)
            except AttributeError:
                map(pydesc.structure.AbstractStructure.set_contact_map, self.structures)
                cm1, cm2 = map(operator.attrgetter("contact_map"), self.structures)
        cstc1, cstc2 = self.get_specified_structures()
        alg = dict(self.aligned_mers)
        ralg = dict(map(reversed, self.aligned_mers))
        for m1, m2 in self.aligned_mers:
            for m1c in cm1.contacts[m1.ind]:
                try:
                    conts.add((pydesc.structure.Contact(pydesc.structure.ElementChainable(m1),\
                                                        pydesc.structure.ElementChainable(cstc1[m1c])
                                                        ),
                               pydesc.structure.Contact(pydesc.structure.ElementChainable(m2),\
                                                        pydesc.structure.ElementChainable(alg[cstc1[m1c]])
                                                        )
                              )
                             )
                except (ValueError, KeyError):
                    pass
            for m2c in cm2.contacts[m2.ind]:
                try:
                    conts.add((pydesc.structure.Contact(m1.ind, ralg[cstc2[m2c]].ind, cstc1), pydesc.structure.Contact(m2.ind, m2c, cstc2)))
                except (ValueError, KeyError):
                    pass
        def overfit(o1, o2):
            oft = Overfit()
            oft.add_structure(o1, o2)
            return oft.overfit()[0]
        rmsds = [overfit(c1, c2) for c1, c2 in conts]
        try:
            return sum(rmsds) / len(rmsds)
        except ZeroDivisionError:
            ZeroDivisionError("There is no contacts aligned.")

class AlignedTuples(Alignment):

    """Class representing multiple alignment as list of aligned mers from all aligned structures."""

    _state_attr = ['structures', 'aligned_mers']

    def __init__(self, structures, list_of_aligned_mers):
        """AlignedTuples constructor.

        Arguments:
        structures -- sequence of pydesc.structure.AbstractStructures instances.
        list_of_aligned_mers -- sequence of tuples containing aligned mers from appropriate structures given in structures list.
        """
        self.structures = structures
        self.aligned_mers = list_of_aligned_mers

    def _fill_from_MergedPairAlignments(self, obj):   # pylint: disable=invalid-name
        # method name is as described in contextmeta.py documentation
        """Method required and called by ContextStateMeta - metaclass of MultipleAlignment.

        Argument:
        obj -- instance of MergedPairAlignments class.

        See contextmeta documentation for more information.
        """
        self.structures = set(reduce(operator.add, [pair_alignment.structures for pair_alignment in obj.alignments]))
        self.structures = list(self.structures)
        self.aligned_mers = set()
        for pair in obj.iterate_pairs():
            t_tup = ['-' for dummy in self.structures]
            for mer in list(obj._hash(pair[0])) + [pair[0]]:    # pylint:disable=protected-access
                # that is the method that should use protected attr
                t_tup[self.structures.index(mer.structure)] = mer
            self.aligned_mers.add(tuple(t_tup))
        self.aligned_mers = list(self.aligned_mers)

    def _fill_from_GraphsAndColumns(self, obj):  # pylint:disable=invalid-name
        # contextmeta demands that name
        """Method required and called by ContextStateMeta - metaclass of MultipleAlignment.

        Argument:
        obj -- instance of GraphsAndColumns class.

        See contextmeta documentation for more information.
        """
        self.structures = obj.structures
        self.aligned_mers = []

        def mk_clq(p_clq):
            """Makes clicques.

            Argument:
            p_clq -- putative clicque; list of mers to be tested for being clicque.

            If given list is a clicque, turnes it into tuple and appends clicque list,
            otherwise call mk_clq for all slices lacking mers that are not connected with
            all verticles (other mers).
            """
            acc_n = 0
            for mer in p_clq:
                if all(m in graph[mer] or m == mer for m in p_clq):  # pylint:disable=undefined-loop-variable
                    # this variable is defined for sure
                    acc_n += 1
                    continue
                mk_clq([i for i in p_clq if i != mer])
                mk_clq([i for i in graph[mer] if i in p_clq] + [mer])  # pylint:disable=undefined-loop-variable
            if len(p_clq) != acc_n:
                return
            clq = ["-" for i in self.structures]
            for mer in p_clq:
                clq[self.structures.index(mer.structure)] = mer
            clicques.add(tuple(clq))

        for ory_col, graph in zip(obj.columns, obj.graphs):
            col = [i for i in ory_col if i != '-']
            if all(len(val) + 1 == len(col) for val in graph.values()):
                # when every vertices has no. of vertices - 1 neighbours - column equals tuple
                self.aligned_mers.append(tuple(ory_col))
                continue

            clicques = set()
            mk_clq(col)
            self.aligned_mers.extend(list(clicques))

    def __iter__(self):
        """Returns iterator that itartaes over list of aligned mers."""
        return iter(self.aligned_mers)

    def __repr__(self):
        return "<AlignedTuples of %i structures (%i merged mers)>" % (len(self.structures), len(self.aligned_mers))

    def close(self):
        """Returns closed copy of current alignment.

        If any mer from one structure is aligned with mers from more than one structure - this method merges all alignments in one alignment, if it is possible.
        Closed AlignedTupless generates all possible PairAlignments when changed to MultiplePairAlignments.
        """
        all_aligned_mers = dict((mer, None) for aligned_mers in self for mer in aligned_mers if mer != "-")
        hash_ = []
        for aligned_mers in self:
            aligned_mers = tuple(sorted([mer for mer in aligned_mers if mer != "-"]))
            mers_to_change_index = []
            for mer in aligned_mers:
                c_index = all_aligned_mers[mer]
                # c_index stands for current index - value of current mer in
                # all_aligned_mers dict
                if c_index is not None:
                    # means that mer has been seen in any tuple of aligned mers
                    extra_mers = hash_[c_index]
                    mers_to_change_index.extend(extra_mers)
                    aligned_mers = tuple(set(aligned_mers + extra_mers))
                mers_to_change_index.append(mer)
                # len(hash_) is index of alignment tuple recently added to hash_
                # list
            for mer in aligned_mers:
                all_aligned_mers[mer] = len(hash_)
            hash_.append(aligned_mers)
        structures_indexes = dict((structure_obj, self.structures.index(structure_obj)) for structure_obj in self.structures)
        gl_aligned_mers = []
        for alignment_tuple in reversed(hash_):
            if all_aligned_mers[alignment_tuple[0]] is not None:
                aligned_mers = ["-" for dummy in self.structures]
                for mer in alignment_tuple:
                    all_aligned_mers[mer] = None
                    aligned_mers[structures_indexes[mer.structure]] = mer
                gl_aligned_mers.append(aligned_mers)
        return Alignment.build_from_list_of_mers([i for i in self.structures], gl_aligned_mers)

    def get_index(self, mer):
        """Returns index of alignment tuple in which given mere is present.

        Argument:
        mer -- pydesc.monomer.Monomer subclass instance.

        Raises KeyError if mer is not present.
        """
        for n, mrs in enumerate(self.aligned_mers):
            if mer in mrs: return n
        raise KeyError('Given instance of %s is not present in alignment %s' % (str(mer.pid), str(self)))


class MergedPairAlignments(Alignment):

    """Class representing multiple alignment as list of (all possible) pair alignment between aligned structures."""

    _state_attr = ['alignments']

    def __init__(self, alignments):
        """MergedPairAlignment constructor.

        Argument:
        alignments -- sequence of PairAlignment instances.

        Raises ValueError if given sequence is shorter than 2.
        """
        if len(alignments) < 2:
            raise ValueError("Cannot create multiple alignment: less than two alignments given")
        self.alignments = list(alignments)

    def __repr__(self):
        return "<Merged %i alignments>" % len(self.alignments)

    def __str__(self):
        """Prints all stored PairAlignments one by one."""
        return "alignment:\n".join(map(str, self.alignments))

    def _hash(self, mer):
        """Returns mer aligned with given mer."""
        algd = []
        for alg in self.alignments:
            try:
                algd.append(alg._hash(mer))  # pylint:disable=protected-access
            except KeyError:
                pass
        return tuple(algd)

    def _fill_from_AlignedTuples(self, obj):   # pylint: disable=invalid-name
        # method name is as described in contextmeta.py documentation
        """Method required and called by ContextStateMeta - metaclass of MultipleAlignment.

        Argument:
        obj -- instance of AlignedTuples class.

        See contextmeta documentation for more information.
        """
        self.alignments = []
        for i_1, structure_1 in enumerate(obj.structures):
            for structure_2 in obj.structures[i_1 + 1:]:
                i_2 = obj.structures.index(structure_2)
                aligned_pairs = []
                for aligned_mers in obj.aligned_mers:
                    if aligned_mers[i_1] == "-" or aligned_mers[i_2] == "-":
                        continue
                    aligned_pairs.append((aligned_mers[i_1], aligned_mers[i_2]))
                if len(aligned_pairs) != 0:
                    self.alignments.append(PairAlignment([structure_1, structure_2], aligned_pairs))

    def _fill_from_GraphsAndColumns(self, obj):  # pylint:disable=invalid-name
        """Method required and called by ContextStateMeta - metaclass of MultipleAlignment.

        Argument:
        obj -- instance of GraphsAndColumns class.

        See contextmeta documentation for more information.
        """
        algs = dict(((str1, str2), []) for i, str1 in enumerate(obj.structures) for str2 in obj.structures[i + 1:])
        for graph in obj.graphs:
            for mer1, aligned in graph.items():
                for mer2 in aligned:
                    try:
                        algs[(mer1.structure, mer2.structure)].append((mer1, mer2))
                    except KeyError:
                        pass
                        # to avoid duplicated tuples in pair alignments
        self.alignments = [PairAlignment(key, aligned_mers) for key, aligned_mers in algs.items()]

    def intersect_with(self, alignment_obj):
        """Returns narrowed alignment that aligns all structures common for both, current and given, alignments.

        Argument:
        alignment_obj -- instace of Alignment subclasses.

        Returns alignment represented by one of Alignment subclasses
        """
        current_alignment_structures = list()
        for pair_alignment in self.alignments:
            current_alignment_structures.extend([structure_obj for structure_obj in pair_alignment.structures if structure_obj not in current_alignment_structures])
        alignments = [pair_alignment for pair_alignment in alignment_obj.alignments if all(structure_obj in current_alignment_structures for structure_obj in pair_alignment.structures)]
        return Alignment.build_from_pair_alignments(alignments)

    def iterate_pairs(self):
        """Returns generator that iterates over aligned pairs of mers in subsequent PairAlignments."""
        for alg in self.alignments:
            for pair in alg:
                yield pair

    def project(self, list_of_structures):
        """Returns alignment with domain narrowed to given structures.

        Argument:
        list_of_structures -- sequence of pydesc.structure.AbstractStructure instances.

        Returns alignment represented by one of Alignment subclasses
        """
        validated_alignments = [alignment_obj for alignment_obj in self.alignments if all(structure_obj in list_of_structures for structure_obj in alignment_obj.structures)]
        # taking only those alignments, whose structures occure in given list of structures
        return Alignment.build_from_pair_alignments(validated_alignments)

    def subtract_domains(self, alignment_obj):
        """Returns alignment with domain narrowed to structures that are NOT present in given alignment.

        Argument:
        alignment_obj -- instance of any Alignment subclass.

        Returns alignment represented by one of Alignment subclasses
        """
        def discard_common(pair_alignment):
            """Returns False if given PairAlignment is common for given and current alignments."""
            if any(structure_obj in alignment_obj.structures for structure_obj in pair_alignment.structures):
                return False
            return True

        filtered = filter(discard_common, self.alignments)
        if len(filtered) > 2:
            return MultipleAlignment(alignments=filtered)   # pylint: disable=no-value-for-parameter, unexpected-keyword-arg
        return filtered[0]

    def get_alignment_for(self, stc1, stc2):
        """Returns PairAlignment for given pair of structures.

        stc1, stc2 -- structure objects or their names.
        """
        nms = set([])
        for i in (stc1, stc2):
            try:
                nms.add(i.name)
            except AttributeError:
                nms.add(i)
        return max([i for i in self.alignments if set(map(operator.attrgetter('name'), i.structures)) == nms])


class GraphsAndColumns(object):

    """Class representing """

    _state_attr = ['columns', 'graphs', 'structures']

    def _fill_from_AlignedTuples(self, obj):    # pylint:disable=invalid-name
        """Method required and called by ContextStateMeta - metaclass of MultipleAlignment.

        Argument:
        obj -- instance of AlignedTuples class.

        See contextmeta documentation for more information.
        """
        closed_obj = obj.close()
        graphs = [{} for i in closed_obj.aligned_mers]
        hash_ = dict((mer, i) for i, col in enumerate(closed_obj.aligned_mers) for mer in col)
        for aligned_mers in obj.aligned_mers:
            aligned_mers = [i for i in aligned_mers if i != '-']
            for mer_obj in aligned_mers:
                col_hash = hash_[mer_obj]
                try:
                    graphs[col_hash][mer_obj].extend([i for i in aligned_mers if i != mer_obj])
                except KeyError:
                    graphs[col_hash][mer_obj] = [i for i in aligned_mers if i != mer_obj]
                # graphs is a list of dicts containing mers as keys and list of merf from other strucutures
                # which they are aligned with as keys.
                # hash_ is used to get index of proper dict on graphs list, so filling of dicts could be performed.
        self.structures = obj.structures    # pylint:disable=attribute-defined-outside-init
        self.columns = closed_obj.aligned_mers  # pylint:disable=attribute-defined-outside-init
        self.graphs = graphs    # pylint:disable=attribute-defined-outside-init
        # method replaces init somehow

    def _fill_from_MergedPairAlignments(self, obj):  # pylint:disable=invalid-name
        """Method required and called by ContextStateMeta - metaclass of MultipleAlignment.

        Argument:
        obj -- instance of MergedPairAlignments class.

        See contextmeta documentation for more information.
        """
        self._fill_from_AlignedTuples(obj)
        # state is changed due to usage of AlignedTuples attributes in _fill_from method


class MultipleAlignment(AlignedTuples, MergedPairAlignments, GraphsAndColumns):

    """Class representing alignment between more than two structures. Metaclass is ContextStateMeta."""

    __metaclass__ = contextmeta.ContextStateMeta
