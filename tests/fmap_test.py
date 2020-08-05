import numpy as np
import os.path

from pydesc.api.cmaps import calculate_contact_map
from pydesc.api.cmaps import create_frequency_map_from_contact_maps
from pydesc.api.structure import get_structures_from_file


def test_fmap_from_cmap(structures_dir):
    path = os.path.join(structures_dir, "prots_only_nmr", "2LJP.pdb")
    structures = get_structures_from_file(path)

    cms = [calculate_contact_map(stc) for stc in structures]

    fmap = create_frequency_map_from_contact_maps(cms)
    pair = (0, 1)
    for pair, freq in fmap:
        assert 0.0 <= freq <= 1.0
    assert fmap.get_dok_matrix().dtype == np.float64
    freq = fmap.get_contact_frequency(*pair)
    occ = fmap.get_contact_occurs(*pair)
    assert freq == (occ / fmap.n_frames)
    assert fmap.n_frames == 20
    freqs = fmap.get_contacts_frequencies(pair[0])
    occs = fmap.get_contacts_occurs(pair[0])
    assert len(freqs) == len(occs)
    assert [i[0] for i in occs] == [j[0] for j in freqs]
