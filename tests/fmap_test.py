import numpy as np
import os.path

from pydesc.api.cmaps import calculate_contact_map
from pydesc.api.cmaps import create_frequency_map_from_contact_maps
from pydesc.api.structure import get_structures_from_file
from pydesc.selection import ChainSelection

from pydesc.structure import TrajectoryLoader
from pydesc.contacts.maps import FrequencyMapCalculator
from pydesc.api.criteria import get_default_protein_criterion, get_gc_distance_criterion


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


def test_fmap_calculator(trajectories_path, topologies_path):
    topo_path = topologies_path.format(topo="mdm2")
    traj_path = trajectories_path.format(type="dcd", traj="mdm2_5frames")

    topo = get_structures_from_file(topo_path)[0]
    traj = TrajectoryLoader().load_trajectory(traj_path, topo)

    crit = get_default_protein_criterion()
    calc = FrequencyMapCalculator(traj, crit)
    fmap = calc.calculate_frequency_map()

    frequencies = {i[-1] for i in fmap}

    assert 0 < min(frequencies) <= 1.0
    assert 0 < max(frequencies) <= 1.0

    contacts_mer_2 = fmap.get_contacts_frequencies(2)
    mers_in_contact_w_2 = [i[0] for i in contacts_mer_2]
    assert 3 in mers_in_contact_w_2
    assert 84 in mers_in_contact_w_2

    assert dict(contacts_mer_2)[3] == 1.0
    assert dict(contacts_mer_2)[84] == 0.1


def test_partial_fmap(trajectories_path, topologies_path):
    topo_path = topologies_path.format(topo="trf1-dna")
    traj_path = trajectories_path.format(type="xtc", traj="trf1-dna_5frames")

    topo = get_structures_from_file(topo_path)[0]
    traj = TrajectoryLoader().load_trajectory(traj_path, topo)

    crit = get_default_protein_criterion()
    calc_whole = FrequencyMapCalculator(traj, crit)
    chain_a = ChainSelection('A').create_structure(traj)     # dna chain
    calc_chain_a = FrequencyMapCalculator(chain_a, crit)
    chain_b = ChainSelection('B').create_structure(traj)     # protein chain
    calc_chain_b = FrequencyMapCalculator(chain_b, crit)

    fmap_chain_a = calc_chain_a.calculate_frequency_map()
    fmap_chain_b = calc_chain_b.calculate_frequency_map()
    fmap_whole = calc_whole.calculate_frequency_map()

    assert len(fmap_chain_a) == 0
    assert list(fmap_chain_b) == list(fmap_whole)


def test_inter_fmap(trajectories_path, topologies_path):
    topo_path = topologies_path.format(topo="trf1-dna")
    traj_path = trajectories_path.format(type="xtc", traj="trf1-dna_5frames")

    topo = get_structures_from_file(topo_path)[0]
    traj = TrajectoryLoader().load_trajectory(traj_path, topo)

    crit = get_gc_distance_criterion()
    selections = (
        ChainSelection("A"),
        ChainSelection("B"),
    )
    calc = FrequencyMapCalculator(traj, crit, selections)
    fmap = calc.calculate_frequency_map()

    assert len(fmap) == 10
    for _, value in fmap:
        assert 0.0 <= value <= 1.0
