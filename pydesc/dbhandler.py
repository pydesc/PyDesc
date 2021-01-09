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

import gzip
import os.path
import tarfile
from contextlib import contextmanager
from io import StringIO
from urllib import request
from urllib.error import HTTPError

from Bio.PDB import PDBParser

from pydesc.config import ConfigManager
from pydesc.warnexcept import Info
from pydesc.warnexcept import warn

try:
    from Bio.PDB import MMCIFParser
except ImportError:

    def MMCIFParser(*a, **ka):
        return None

    warn(Info("No MMCIFParser in Bio.PDB"))

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


class OperationModeError(Exception):
    pass


class ContextManagerMixIn:
    @contextmanager
    def open(self, val, *args, **kwargs):
        files = self.get_file(val, *args, **kwargs)
        yield files
        for file_ in files:
            file_.close()


class DBHandler(ContextManagerMixIn):
    def __init__(self, mode):
        self.mode = mode

    def get_cache(self, val):
        return ConfigManager.dbhandler.cachedir + self.db_name + "/" + val[1:3] + "/"

    def is_id_valid(self, val):
        return True

    def is_file_valid(self, fh):
        return True

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
    def get_file(self, val):
        """Returns a handler to the file containing structure of choice.

        Arguments:
        val -- structure id.
        """
        dct = {
            3: (self.assert_val, Info(f"Accessing cache to load {val}...")),
            2: (self.get_from_local_db, Info(f"Accessing local db to load {val}...")),
            1: (self.download_file, Info(f"Downloading {val} to cache...")),
        }
        for i in self.mode:
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
    def __init__(self, mode=(1, 3)):
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
            content = request.urlopen(self.get_file_url(val, unit))
        except HTTPError as e:
            if e.getcode() == 404:
                raise InvalidID(4)
            raise

        f_name = self.get_file_location(val, unit)
        dir_name = os.path.dirname(f_name)

        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)

        buf = content.read()

        if not self.is_file_valid(buf):
            raise InvalidID(5)

        local_file = open(self.get_file_location(val, unit), "wb")
        local_file.write(buf)
        local_file.close()

    @validate_id
    def get_file(self, val, unit=None):
        if unit is not None:
            return self._get_file_unit(val, unit)
        unit = 1
        fh_list = []
        while True:
            unit_handlers = self._get_file_unit(val, unit)
            if unit_handlers is None:
                if not fh_list:
                    msg = (
                        f"No local {val} biounit files. "
                        f"Check path or change access mode."
                    )
                    raise ValueError(msg)
                return fh_list
            fh_list.extend(unit_handlers)
            unit += 1

    def _get_file_unit(self, val, unit):
        result = []
        if 3 in self.mode:
            warn(Info(f"Accessing local copy of {val}/{unit}..."))
            try:
                fh = gzip.open(self.get_file_location(val, unit), "r")
            except FileNotFoundError:
                warn(Info("No local copy found."))
            else:
                result.append(fh)
                warn(Info("Done."))
                return result
        if 1 in self.mode:
            warn(Info(f"Downloading {val}/{unit}..."))
            try:
                self.download_file(val, unit)
            except Exception as e:
                warn(Info(f"Failed due to {type(e)}:{str(e)}"))
            else:
                fh = gzip.open(self.get_file_location(val, unit), "r")
                result.append(fh)
                warn(Info("Done."))
                return result
        if 2 in self.mode:
            raise OperationModeError("BioUnit handler cannot operate in mode 2.")


class MetaHandler(ContextManagerMixIn):
    def __init__(self, mode=(3, 2, 1)):
        self.mode = mode

    def get_file(self, val):
        db, dummy, filename = val.partition("://")
        db = db.lower()

        db_tuple = ("pdb+mmCIF", "cath", "scop")
        mode = self.mode
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
                    return BioUnitHandler().get_file(filename, None)
                return BioUnitHandler().get_file(filename, int(unit))

            # when source is given and it is NOT unit
            for handler in db_dict[db]:
                try:
                    return handler.get_file(filename)
                except InvalidID:
                    continue
            raise InvalidID

        # if no source -- lookup goes in order given by db_tuple
        filename = db
        for i, db in enumerate(db_tuple):
            for handler in db_dict[db]:
                try:
                    return handler.get_file(filename)
                except Exception as e:
                    msg = f"Failed to load structure {val} from bd {db}"
                    info = Info(msg)
                    warn(info)
                    continue
        raise InvalidID("%s" % db)


class MetaParser:
    def __init__(self, *args, **kwargs):
        self.parsers = [PDBParser(*args, **kwargs), MMCIFParser(*args, **kwargs)]

    def get_structure(self, stc, file, *args, **kwargs):
        for prsr in self.parsers:
            try:
                return prsr.get_structure(stc, file, *args, **kwargs)
            except ValueError:  # the only exception known to be raised when
                # proper mmCIF is passed to PDBParser
                file.seek(0)
                continue
        raise ValueError(
            "None of parsers could get %s structure. Tried: %s."
            % (stc, ", ".join([type(i).__name__ for i in self.parsers]))
        )
