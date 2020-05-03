from shutil import rmtree

import pytest

from pydesc.config import ConfigManager
from pydesc.dbhandler import PDBHandler


class TestPDBHandler:

    @pytest.fixture(scope="function", autouse=True)
    def clear_cache(self):
        cache_dir = ConfigManager.dbhandler.cachedir
        rmtree(cache_dir)
        yield
        rmtree(cache_dir)

    def test_download(self):
        handler = PDBHandler(mode=[1])

        files = handler.get_file("1no5")

        assert len(files) == 1
        reading = files[0].read()
        assert 1672 == reading.count("ATOM")
