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

"""Database handler responsible for downloading and reading cache of local structure
files."""

import gzip
import os.path
import tarfile
from abc import ABC
from abc import abstractmethod
from contextlib import contextmanager
from functools import wraps
from io import StringIO
from pathlib import Path
from urllib import request
from urllib.error import HTTPError

from pydesc.config import ConfigManager
from pydesc.warnexcept import Info
from pydesc.warnexcept import InvalidID
from pydesc.warnexcept import OperationModeError
from pydesc.warnexcept import warn

ConfigManager.new_branch("dbhandler")
ConfigManager.dbhandler.set_default("cachedir", "./biodb/")
ConfigManager.dbhandler.set_default("pdb_handler", "./biodb/")


def add_db_dir(method):
    """Glues path to db's parental directory before path returned by original method."""

    @wraps(method)
    def wrapped_method(self, val, *args, **kwargs):
        return self.get_db_parent_dir() / method(self, val, *args, **kwargs)

    return wrapped_method


def validate_id(method):
    """Calls entry id validation before calling wrapped method."""

    @wraps(method)
    def wrapped_method(self, entry_id, *args, **kwargs):
        if not self.is_id_valid(entry_id):
            handler_type = type(self).__name__
            msg = f"Invalid entry id ({entry_id} passed to db handler ({handler_type})."
            raise InvalidID(msg)
        return method(self, entry_id, *args, **kwargs)

    return wrapped_method


class ContextManagerMixIn:
    """MixIn for database handlers making them usable as context managers.

    Adds method *open*, which takes val and any additional args and passes them to
    *get_file* method.
    It assumes result is a list of file handlers.
    They are guaranteed to be closed when closing context manager.
    """

    @contextmanager
    def open(self, val, *args, **kwargs):
        """Return context manager calling *get_file* method with passed arguments,
        that closes all file handlers returned by this method while leaving context
        manager."""
        files = self.get_file(val, *args, **kwargs)
        yield files
        for file_ in files:
            file_.close()


class DBHandler(ContextManagerMixIn, ABC):
    """Abstract database handler.

    When extending, make sure to implement:
    * get_entry_path(entry_id) method;
    * db_name class attribute;
    * url_template class attribute or overwrite get_entry_url(entry_id) method.

    * is_id_valid(entry_id) method by default returns True, so it is also good idea to
    overwrite this method.

    Overwriting or extending other methods is optional.
    Default implementation is often suitable for database giving single file per entry.

    """

    def __init__(self, mode):
        self.mode = mode

    def get_db_parent_dir(self):
        """Return path to db directory in local cache."""
        cache_dir = Path(ConfigManager.dbhandler.cachedir)
        path = cache_dir / self.db_name
        return path

    def is_id_valid(self, entry_id):
        """Return True if entry id is valid."""
        return True

    def assert_entry_downloaded(self, entry_id):
        """Raise IOError if given entry has no file in local cache."""
        path = self.get_entry_path(entry_id)
        if not path.is_file():
            raise IOError("No such file: %s" % path)

    @abstractmethod
    def get_entry_path(self, entry_id):
        """Get path to a file with structure of given id.
        This method just creates a path, getting a result does not mean a file
        exists."""
        pass

    @validate_id
    def get_entry_url(self, entry_id):
        """Fill url template with given structure id."""
        return self.url_template % entry_id

    @validate_id
    def download_file(self, entry_id):
        """Downloads content and write it to a local file."""
        try:
            u = request.urlopen(self.get_entry_url(entry_id))
        except HTTPError as e:
            if e.getcode() == 404:
                raise InvalidID(2)
            raise
        reading = StringIO(u.read().decode("utf-8"))
        self.save_stream(reading, entry_id)

    @validate_id
    def get_from_local_db(self, entry_id):
        """Method returning file for given structure id from local database.

        Raises:
            NotImplemented if handler cannot serve local db.
        """
        raise NotImplemented

    def save_stream(self, stream, entry_id):
        """Write given stream to a file in local cache representing structure of
        given id."""
        file_name = self.get_entry_path(entry_id)
        dir_name = os.path.dirname(file_name)
        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)
        buffer = stream.read()
        with open(self.get_entry_path(entry_id), "w") as local_file:
            local_file.write(buffer)

    @validate_id
    def get_file(self, entry_id):
        """Returns a handler to the file containing structure of choice.

        Args:
            entry_id(str): structure id.

        Returns:
            list: of file handlers.

        """
        dct = {
            3: (
                self.assert_entry_downloaded,
                Info(f"Accessing cache to load {entry_id}..."),
            ),
            2: (
                self.get_from_local_db,
                Info(f"Accessing local db to load {entry_id}..."),
            ),
            1: (self.download_file, Info(f"Downloading {entry_id} to cache...")),
        }
        for i in self.mode:
            try:
                mth, info = dct[i]
                warn(info, 4)
                mth(entry_id)
                warn(Info("Done."), 4)
            except NotImplementedError:
                raise
            except Exception as e:
                warn(Info(f"Failed (due to {type(e).__name__}: {e})."))
                continue
            else:
                return [self._get_file(entry_id)]
        raise IOError("No file to load for %s" % entry_id)

    def _get_file(self, entry_id):
        """Method returning a file.

        Meant to be overwritten in db handlers that should return more than one handler.

        Args:
            entry_id(str): structure id.

        Returns:
            : an open file handler.
        """
        return open(self.get_entry_path(entry_id), "r")


class SCOPHandler(DBHandler):
    """Handler providing PDB files from SCOP database."""

    def __init__(self, mode=0):
        self.db_name = "scop"
        self.url_template = (
            "http://scop.berkeley.edu/astral/pdbstyle/?ver=2.03&id=%s&output=pdb"
        )
        super().__init__(mode)

    def is_id_valid(self, entry_id):
        if len(entry_id) != 7:
            return False

        if (
            entry_id[0] != "d"
            or not entry_id[1].isdigit()
            or not entry_id[2:4].isalnum()
            or not (entry_id[5].isalnum() or entry_id[5] in ["_", "."])
            or not (entry_id[6].isalnum() or entry_id[6] == "_")
        ):
            return False
        return super().is_id_valid(entry_id)

    @validate_id
    @add_db_dir
    def get_entry_path(self, entry_id):
        super().get_entry_path(entry_id)
        return "%s/%s.ent" % (entry_id[2:4], entry_id)


class PDBHandler(DBHandler):
    """Handler providing PDB files from RCSB Protein Data Bank."""

    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "pdb"
        self.url_template = "http://www.rcsb.org/pdb/files/%s.pdb"
        super().__init__(mode)

    def is_id_valid(self, entry_id):
        if len(entry_id) != 4:
            return False
        if not entry_id[0].isdigit() or not entry_id[1:3].isalnum():
            return False
        return super().is_id_valid(entry_id)

    def get_file(self, entry_id, *args, **kwargs):
        result = super().get_file(entry_id.lower(), *args, **kwargs)
        return result

    @validate_id
    @add_db_dir
    def get_entry_path(self, entry_id):
        super().get_entry_path(entry_id)
        return "%s/%s.pdb" % (entry_id[1:3], entry_id)

    @validate_id
    def get_from_local_db(self, entry_id):
        with gzip.open(
            ConfigManager.dbhandler.pdb_handler
            + "data/structures/divided/pdb/%s/pdb%s.ent.gz" % (entry_id[1:3], entry_id)
        ) as f:
            self.save_stream(f, entry_id)


class MMCIFHandler(PDBHandler):
    """Handler providing mmCIF files from RCSB Protein Data Bank."""

    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "mmCIF"
        self.url_template = "http://www.rcsb.org/pdb/files/%s.cif"
        super().__init__(mode)

    @validate_id
    @add_db_dir
    def get_entry_path(self, entry_id):
        super().get_entry_path(entry_id)
        return "%s/%s.cif" % (entry_id[1:3], entry_id)


class PDBBundleHandler(DBHandler):
    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "pdb-bundle"
        self.url_template = "ftp://ftp.wwpdb.org/pub/pdb/compatible/pdb_bundle/"
        DBHandler.__init__(self, mode)

    def get_bundle(self, val):
        cache = self.get_db_parent_dir(val.lower())
        with open(cache + val.lower() + "-chain-id-mapping.txt") as f:
            bundle = [
                i.replace(":", "") for i in map(str.strip, f.readlines()) if ".pdb" in i
            ]
        return bundle

    def get_entry_path(self, entry_id):
        cache = self.get_db_parent_dir(entry_id)
        bundle = self.get_bundle(entry_id)
        return [cache + i for i in bundle]

    @validate_id
    def download_file(self, entry_id):
        raise NotImplementedError

    def get_from_local_db(self, entry_id):
        cache = self.get_db_parent_dir(entry_id)
        dir = ConfigManager.dbhandler.pdb_handler + "compatible/pdb_bundle/"
        tar = tarfile.open(
            dir
            + "%s/%s/%s-pdb-bundle.tar.gz"
            % tuple(map(str.lower, (entry_id[1:3], entry_id, entry_id))),
            "r:gz",
        )
        try:
            os.makedirs(cache)
        except OSError:
            pass
        tar.extractall(cache)

    def assert_entry_downloaded(self, entry_id):
        cache = self.get_db_parent_dir(entry_id)
        bundle = self.get_bundle(entry_id)
        for file in bundle:
            if file not in os.listdir(cache):
                raise IOError("Lacking files for %s bundle." % entry_id)

    def _get_file(self, entry_id):
        """Returns list of handlers to bundle of pdb files."""
        bundle = self.get_bundle(entry_id)
        return [open(i) for i in self.get_entry_path(entry_id)]

    def get_mapping(self, val):
        """Returns dict of mmCIF chain names as keys and tuples containing pdb file name and pdb
         chain name as values.
         """
        cache = self.get_db_parent_dir(val.lower())
        dct = {}
        with open(cache + val.lower() + "-chain-id-mapping.txt") as f:
            for blk in f.read().split("\n\n")[1:]:
                ls = blk.split("\n")
                pdbn = ls[0].strip()[:-1]
                dct[pdbn] = dict(list(map(str.split, list(filter(bool, ls[1:])))))
        return dct


class CATHHandler(DBHandler):
    """Handler providing PDB files from CATH database."""

    def __init__(self, mode=(1, 2, 3)):
        self.db_name = "CATH"
        # unpublished path, subject to change w/o notice
        self.url_template = "http://data.cathdb.info/v4_0_0/pdb/%s"
        super().__init__(mode)

    def is_id_valid(self, entry_id):
        if len(entry_id) != 4:
            return False
        if not entry_id[0].isdigit() or not entry_id[1:3].isalnum():
            return False
        return super().is_id_valid(entry_id)

    @validate_id
    @add_db_dir
    def get_entry_path(self, entry_id):
        super().get_entry_path(entry_id)
        return "%s/%s.pdb" % (entry_id[2:4], entry_id)


class BioUnitHandler(DBHandler):
    """Handler providing PDB files from RCSB Protein Data Bank BioUnit."""

    def __init__(self, mode=(1, 3)):
        self.db_name = "PDBBioUnit"
        self.url_template = "http://www.rcsb.org/pdb/files/%s.pdb%d.gz"
        DBHandler.__init__(self, mode)

    def is_id_valid(self, entry_id):
        if len(entry_id) != 4:
            return False
        if not entry_id[0].isdigit() or not entry_id[1:3].isalnum():
            return False
        return super().is_id_valid(entry_id)

    @validate_id
    @add_db_dir
    def get_entry_path(self, entry_id, unit):
        super().get_entry_path(entry_id)
        return "%s/%s.pdb%d.gz" % (entry_id[2:4], entry_id, unit)

    @validate_id
    def get_entry_url(self, entry_id, unit):
        """Get URL to entry with unit number."""
        return self.url_template % (entry_id, unit)

    @validate_id
    def download_file(self, entry_id, unit):
        """Download a file for given entry and unit number."""
        try:
            content = request.urlopen(self.get_entry_url(entry_id, unit))
        except HTTPError as e:
            if e.getcode() == 404:
                raise InvalidID(4)
            raise

        f_name = self.get_entry_path(entry_id, unit)
        dir_name = os.path.dirname(f_name)

        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)

        buf = content.read()

        local_file = open(self.get_entry_path(entry_id, unit), "wb")
        local_file.write(buf)
        local_file.close()

    @validate_id
    def get_file(self, entry_id, unit=None):
        """Get a file for given entry id and unit number, or set unit to None (
        default) to get all units."""
        if unit is not None:
            return self._get_file_unit(entry_id, unit)
        unit = 1
        fh_list = []
        while True:
            unit_handlers = self._get_file_unit(entry_id, unit)
            if unit_handlers is None:
                if not fh_list:
                    msg = (
                        f"No local {entry_id} biounit files. "
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
                fh = gzip.open(self.get_entry_path(val, unit), "r")
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
                fh = gzip.open(self.get_entry_path(val, unit), "r")
                result.append(fh)
                warn(Info("Done."))
                return result
        if 2 in self.mode:
            raise OperationModeError("BioUnit handler cannot operate in mode 2.")


class MetaHandler(ContextManagerMixIn):
    """Handler running over all possible handlers and trying out all of them.

    Interface of this class provides get_file method similar to other handlers.
    It accepts entry_id with slightly extended syntax though.
    """

    def __init__(self, mode=(3, 2, 1)):
        self.mode = mode

    def get_file(self, entry_id):
        """Try to get file handler for structure from first database which entry_id
        will match.

        To force using particular handler, use special entry format:
        <handler>://<entry-id>
        Available handlers:
        * pdb -- PDB with PDB files only
        * mmCIF -- PDB with mmCIF files only
        * pdb+mmCIF -- PDB with any file type
        * scop -- SCOP with PDB file only
        * unit -- PDBBioUnit with PDB file only
        E.g.:
        pdb://1no5
        mmCIF://1no5

        In case of 'unit' db, one can additionally add unit after '/', e.g.:
        unit://3pi2/3
        If unit was not passed, all unit were downloaded.

        Args:
            entry_id(str): entry id or entry with db and other special options.

        Returns:
            list: sequence of file handlers to local cache.

        """
        db, dummy, filename = entry_id.partition("://")
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
                    msg = f"Failed to load structure {entry_id} from bd {db}"
                    info = Info(msg)
                    warn(info)
                    continue
        raise InvalidID("%s" % db)
