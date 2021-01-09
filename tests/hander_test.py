from shutil import rmtree

import pytest

from pydesc.config import ConfigManager
from pydesc.dbhandler import BioUnitHandler
from pydesc.dbhandler import MMCIFHandler
from pydesc.dbhandler import PDBHandler
from pydesc.dbhandler import SCOPHandler


class TestPDBHandler:
    @pytest.fixture(scope="function", autouse=True)
    def clear_cache(self):
        cache_dir = ConfigManager.dbhandler.cachedir
        try:
            rmtree(cache_dir)
        except FileNotFoundError:
            pass
        yield
        rmtree(cache_dir)

    def test_download(self):
        handler = PDBHandler(mode=[1])

        files = handler.get_file("1no5")

        assert len(files) == 1
        reading = files[0].read()
        assert 1672 == reading.count("ATOM")

    def test_mmcif_download(self):
        handler = MMCIFHandler(mode=[1])
        files = handler.get_file("1no5")
        assert len(files) == 1

    def test_bio_download(self):
        handler = BioUnitHandler(mode=[1])
        files = handler.get_file("3pi2", 1)
        assert len(files) == 1


class TestSCOPHandler:
    def test_download(self):
        handler = SCOPHandler(mode=[1])
        files = handler.get_file("d1nkla_")
        assert len(files) == 1
