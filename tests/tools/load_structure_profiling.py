import cProfile
import datetime
import os.path

from pydesc.structure import StructureLoader
from tests.conftest import TEST_STRUCTURES_DIR


def run():
    sl = StructureLoader()
    pth = os.path.join(TEST_STRUCTURES_DIR, "prots_only_nmr", "1A24.pdb")
    s = sl.load_structures(path=pth)
    return s


if __name__ == "__main__":
    now = datetime.datetime.now().strftime("%d-%m_%H-%M-%S")
    cProfile.run("run()", "profile_%s" % now)
