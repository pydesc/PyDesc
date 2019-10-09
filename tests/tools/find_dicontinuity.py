from pydesc.structure import StructureLoader as SL

from tests.conftest import PDB_FILES_DICT, TEST_STRUCTURES_DIR
import os.path

sl = SL()

for fname in PDB_FILES_DICT["prots_only"]:
    pth = os.path.join(TEST_STRUCTURES_DIR, "prots_only", fname)
    stc_l = sl.load_structures(path=pth)
    for stc in stc_l:
        for chain in stc.chains:
            chainable = [mer for mer in chain if mer.is_chainable()]
            last_chainable = max(chainable, key=lambda mer: mer.ind)
            for mer, merp1 in zip(chainable, chainable[1:]):
                if mer.next_mer != merp1:
                    info = "%s %s -| |- %s (inds: %i, %i)"
                    tuple_ = (
                        str(chain.name),
                        str(mer.get_pdb_id()),
                        str(merp1.get_pdb_id()),
                        mer.ind,
                        merp1.ind,
                    )
                    print((info % tuple_))
