# Copyright 2017 Pawel Daniluk, Grzegorz Firlik, Tymoteusz Oleniecki
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

from pydesc.warnexcept import *
from urllib import request
from urllib.error import HTTPError
import os.path
import gzip
import tarfile
from Bio.PDB import PDBParser

try:
    from Bio.PDB import MMCIFParser
except ImportError:

    def MMCIFParser(*a, **ka):
        return None

    warn(Info("No MMCIFParser in Bio.PDB"))
from pydesc.config import ConfigManager
from io import StringIO

ConfigManager.new_branch("dbhandler")
ConfigManager.dbhandler.set_default("cachedir", "./biodb/")
ConfigManager.dbhandler.set_default("pdb_handler", "./biodb/")


def add_db_dir(handler):
    def wrap_handler(self, *args, **kwargs):
        return self.get_cache(args[0]) + handler(self, *args, **kwargs)

    return wrap_handler


def validate_id(handler):
    def wrap_handler(self, *args, **kwargs):
        if not self.is_id_valid(args[0]):
            raise InvalidID(1)
        return handler(self, *args, **kwargs)

    return wrap_handler


class InvalidID(Exception):
    pass


class DBHandler:
    def __init__(self, mode):
        self.mode = mode

    def get_cache(self, val):
        return ConfigManager.dbhandler.cachedir + self.db_name + "/" + val[1:3] + "/"

    def is_id_valid(self, val):
        return True

    def is_file_valid(self, fh):
        return True

    @validate_id
    @add_db_dir
    def get_file_location(self, val):
        return self.FileTemplate % val

    def assert_val(self, val):
        pth = self.get_file_location(val)
        if not os.path.exists(pth):
            raise IOError("No such file: %s" % pth)

    @validate_id
    def get_file_url(self, val):
        return self.url_template % val

    @validate_id
    def download_file(self, val):
        try:
            u = request.urlopen(self.get_file_url(val))
        except HTTPError as e:
            if e.getcode() == 404:
                raise InvalidID(2)
            raise
        reading = StringIO(u.read().decode("utf-8"))
        self.save_stream(reading, val)

    @validate_id
    def get_from_local_db(self, val):
        raise NotImplemented

    def save_stream(self, stream, val):
        fname = self.get_file_location(val)

        dirname = os.path.dirname(fname)

        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        buf = stream.read()

        if not self.is_file_valid(buf):
            raise InvalidID(3)

        local_file = open(self.get_file_location(val), "w")
        local_file.write(buf)
        local_file.close()

    @validate_id
    def get_file(self, val, mode=None):
        """Returns a handler to the file containing structure of choice.

        Arguments:
        val -- structure id.
        mode -- list of subsequent access modes to be used by handler, by default set to None. If so initialisation argument is used as mode.
        In mode 1 handler tries to download a new file (and overwrite the on ein biodb directory, if needed).
        Mode 2, if available, works the same way, but gets file from local copy of db. Requires additional settings in ConfigManager.
        In mode 3 handler access file from biodb without earlier attempt at getting file from other sources.
        """
        mode = self.mode if mode is None else mode
        dct = {
            3: (self.assert_val, Info(f"Accessing cache to load {val}...")),
            2: (self.get_from_local_db, Info(f"Accessing local db to load {val}...")),
            1: (self.download_file, Info("Downloading {val} to cache...")),
        }
        for i in mode:
            try:
                mth, info = dct[i]
                warn(info, 4)
                mth(val)
                warn(Info("Done."), 4)
            except NotImplementedError:
                raise
            except Exception as e:
                warn(Info(f"Failed (due to {type(e).__name__}: {e})."))
                continue
            else:
                return [self._get_file(val)]
        raise IOError("No file to load for %s" % val)

    def _get_file(self, val):
        """Private method creating file handler.

        Argument:
        val -- string, structure code.

        Meant to be overwritten in db handlers that should return more than one handler.
        """
        return open(self.get_file_location(val), "r")


class SCOPHandler(DBHandler):
    def __init__(self, mode=0):
        self.db_name = "scop"
        self.url_template = (
            "http://scop.berkeley.edu/astral/pdbstyle/?ver=2.03&id=%s&output=pdb"
        )
        DBHandler.__init__(self, mode)

    def is_id_valid(self, val):
        if len(val) != 7:
            return False

        if (
            val[0] != "d"
            or not val[1].isdigit()
            or not val[2:4].isalnum()
            or not (val[5].isalnum() or val[5] in ["_", "."])
            or not (val[6].isalnum() or val[6] == "_")
        ):
            return False

        return True

    def is_file_valid(self, buf):
        if buf.startswith("ERROR"):
            return False

        return True

    @validate_id
    @add_db_dir
    def get_file_location(self, val):
        return "%s/%s.ent" % (val[2:4], val)


class PDBHandler(DBHandler):
    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "pdb"
        self.url_template = "http://www.rcsb.org/pdb/files/%s.pdb"
        DBHandler.__init__(self, mode)

    def is_id_valid(self, val):
        if len(val) != 4:
            return False

        if not val[0].isdigit() or not val[1:3].isalnum():
            return False

        return True

    def get_file(self, code, *args, **kwargs):
        return DBHandler.get_file(self, code.lower(), *args, **kwargs)

    @validate_id
    @add_db_dir
    def get_file_location(self, val):
        return "%s.pdb" % (val,)

    @validate_id
    def get_from_local_db(self, val):
        with gzip.open(
            ConfigManager.dbhandler.pdb_handler
            + "data/structures/divided/pdb/%s/pdb%s.ent.gz" % (val[1:3], val)
        ) as f:
            self.save_stream(f, val)


class PDBBundleHandler(DBHandler):
    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "pdb-bundle"
        self.url_template = "ftp://ftp.wwpdb.org/pub/pdb/compatible/pdb_bundle/"
        DBHandler.__init__(self, mode)

    def get_bundle(self, val):
        cache = self.get_cache(val.lower())
        with open(cache + val.lower() + "-chain-id-mapping.txt") as f:
            bundle = [
                i.replace(":", "") for i in map(str.strip, f.readlines()) if ".pdb" in i
            ]
        return bundle

    def get_file_location(self, val):
        cache = self.get_cache(val)
        bundle = self.get_bundle(val)
        return [cache + i for i in bundle]

    @validate_id
    def download_file(self, val):
        raise NotImplementedError

    def get_from_local_db(self, val):
        cache = self.get_cache(val)
        dir = ConfigManager.dbhandler.pdb_handler + "compatible/pdb_bundle/"
        tar = tarfile.open(
            dir
            + "%s/%s/%s-pdb-bundle.tar.gz"
            % tuple(map(str.lower, (val[1:3], val, val))),
            "r:gz",
        )
        try:
            os.makedirs(cache)
        except OSError:
            pass
        tar.extractall(cache)

    def assert_val(self, val):
        cache = self.get_cache(val)
        bundle = self.get_bundle(val)
        for file in bundle:
            if file not in os.listdir(cache):
                raise IOError("Lacking files for %s bundle." % val)

    def _get_file(self, val):
        """Returns list of handlers to bundle of pdb files."""
        bundle = self.get_bundle(val)
        return [open(i) for i in self.get_file_location(val)]

    def get_mapping(self, val):
        """Returns dict of mmCIF chain names as keys and tuples containing pdb file name and pdb
         chain name as values.
         """
        cache = self.get_cache(val.lower())
        dct = {}
        with open(cache + val.lower() + "-chain-id-mapping.txt") as f:
            for blk in f.read().split("\n\n")[1:]:
                ls = blk.split("\n")
                pdbn = ls[0].strip()[:-1]
                dct[pdbn] = dict(list(map(str.split, list(filter(bool, ls[1:])))))
        return dct


class MMCIFHandler(PDBHandler):
    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "mmCIF"
        self.url_template = "http://www.rcsb.org/pdb/files/%s.cif"
        DBHandler.__init__(self, mode)

    @validate_id
    @add_db_dir
    def get_file_location(self, val):
        return "%s/%s.cif" % (val[1:3], val)


class CATHHandler(DBHandler):
    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "CATH"
        # unpublished path, subject to change w/o notice
        self.url_template = "http://data.cathdb.info/v4_0_0/pdb/%s"
        DBHandler.__init__(self, mode)

    def is_id_valid(self, val):
        if len(val) != 4:
            return False

        if not val[0].isdigit() or not val[1:3].isalnum():
            return False

        return True

    @validate_id
    @add_db_dir
    def get_file_location(self, val):
        return "%s/%s.pdb" % (val[2:4], val)


class BioUnitHandler(DBHandler):
    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "PDBBioUnit"
        self.url_template = "http://www.rcsb.org/pdb/files/%s.pdb%d.gz"
        DBHandler.__init__(self, mode)

    def is_id_valid(self, val):
        if len(val) != 4:
            return False

        if not val[0].isdigit() or not val[1:3].isalnum():
            return False

        return True

    @validate_id
    @add_db_dir
    def get_file_location(self, val, unit):
        return "%s/%s.pdb%d.gz" % (val[2:4], val, unit)

    @validate_id
    def get_file_url(self, val, unit):
        return self.url_template % (val, unit)

    @validate_id
    def download_file(self, val, unit):
        try:
            u = request.urlopen(self.get_file_url(val, unit))
        except HTTPError as e:
            if e.getcode() == 404:
                raise InvalidID(4)
            raise

        f_name = self.get_file_location(val, unit)
        dir_name = os.path.dirname(f_name)

        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)

        buf = u.read()

        if not self.is_file_valid(buf):
            raise InvalidID(5)

        local_file = open(self.get_file_location(val, unit), "w")
        local_file.write(buf)
        local_file.close()

    @validate_id
    def get_file(self, val, unit, mode=(3, 2, 1)):
        import gzip

        if unit is not None:
            try:
                if mode == 1:
                    raise
                print((self.get_file_location(val, unit)))
                fh = gzip.open(self.get_file_location(val, unit), "r")
                print((val + "_" + str(unit) + ": accessing local copy..."))
            except:
                if mode == 2:
                    raise Exception(
                        "No local "
                        + str(val)
                        + " file. Check path or change access mode."
                    )
                self.download_file(val, unit)
                print(("Downloading " + val + "_" + str(unit) + "..."))
                fh = gzip.open(self.get_file_location(val, unit), "r")
            return [fh]

        else:
            unit = 1
            fh_list = []
            while True:
                try:
                    print(mode)
                    fh_list.append(self.get_file(val, unit, mode))
                    unit += 1
                except:
                    if not fh_list:
                        raise Exception(
                            "No local "
                            + str(val)
                            + " biounit files. Check path or change access mode."
                        )
                    return fh_list


class MetaHandler:
    def __init__(self, mode=(3, 2, 1)):
        self.mode = mode

    def get_file(self, val, mode=(3, 2, 1)):

        db, dummy, filename = val.partition("://")
        db = db.lower()

        # repetitions for testing purposes only
        db_tuple = ("pdb+mmCIF", "cath", "scop")

        mode = self.mode if mode is None else mode

        db_dict = {
            "cath": (CATHHandler(mode),),
            "pdb": (PDBHandler(mode),),
            "mmCIF": (MMCIFHandler(mode),),
            "pdb+mmCIF": (PDBHandler(mode), MMCIFHandler(mode)),
            "scop": (SCOPHandler(mode),),
            "unit": BioUnitHandler(mode),
        }

        if db in db_dict:
            if db == "unit":
                filename, dummy, unit = filename.partition("/")
                if not unit:
                    # when unitDB is to be searched, but no particular unit was given
                    return BioUnitHandler().get_file(filename, None, mode)
                # downloading particular units from unitDB
                return BioUnitHandler().get_file(filename, int(unit), mode)

            # when source is given and it is NOT unit
            for hdlr in db_dict[db]:
                try:
                    return hdlr.get_file(filename, mode)
                except InvalidID:
                    continue
            raise InvalidID

        # if no source -- lookup goes in order given by db_tuple
        else:
            for i, k in enumerate(db_tuple):
                for hdlr in db_dict[k]:
                    try:
                        return hdlr.get_file(db, mode)
                    except Exception as e:
                        continue
            raise InvalidID("%s" % db)


class MetaParser:
    def __init__(self, *args, **kwargs):
        self.parsers = [PDBParser(*args, **kwargs), MMCIFParser(*args, **kwargs)]

    def get_structure(self, stc, file, *args, **kwargs):
        for prsr in self.parsers:
            try:
                return prsr.get_structure(stc, file, *args, **kwargs)
            except ValueError:  # the only exception known to be raised when proper mmCIF is passed to PDBParser
                file.seek(0)
                continue
        raise ValueError(
            "None of parsers could get %s structure. Tried: %s."
            % (stc, ", ".join([type(i).__name__ for i in self.parsers]))
        )
