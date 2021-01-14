from itertools import cycle
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.structure import get_structures_from_file
from pydesc.api.trajectory import freeze_frame
from pydesc.chemistry.factories import BioPythonAtomSetFactory
from pydesc.chemistry.full_atom import MonoatomicIon
from pydesc.chemistry.martini import MartiniResidue
from pydesc.comparison import RMSCLDCalculator
from pydesc.contacts import ContactMapCalculator
from pydesc.contacts.geometrical import PointsDistanceCriterion
from pydesc.selection import AtomSetExactType
from pydesc.selection import ChainSelection
from pydesc.structure import StructureLoader
from pydesc.structure import TrajectoryLoader


def make_frame_with_points(points1, points2):
    frame = MagicMock()
    mers_call1 = []
    for coord in points1:
        point = MagicMock(vector=coord)
        mer = MagicMock(gc=point)
        mer.get_point.return_value = point
        mers_call1.append(mer)
    mers_call2 = []
    for coord in points2:
        point = MagicMock(vector=coord)
        mer = MagicMock(gc=point)
        mer.get_point.return_value = point
        mers_call2.append(mer)
    frame.__getitem__.side_effect = cycle(mers_call1 + mers_call2)
    return frame


class TestUnit:
    @pytest.fixture(scope="function")
    def frame1(self):
        coords = (
            [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0],],
            [[0.0, 3.0, 0.0], [5.0, 3.0, 0.0], [10.0, 3.0, 0.0],],
        )
        frame = make_frame_with_points(*coords)
        return frame

    @pytest.fixture(scope="function")
    def frame2(self):
        coords = (
            [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0],],
            [[0.0, 6.0, 0.0], [5.0, 6.0, 0.0], [10.0, 6.0, 0.0],],
        )
        frame = make_frame_with_points(*coords)
        return frame

    @pytest.fixture(scope="function")
    def frame3(self):
        coords = (
            [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0],],
            [[0.0, 3.2, 0.0], [5.0, 3.3, 0.0], [10.0, 3.1, 0.0],],
        )
        frame = make_frame_with_points(*coords)
        return frame

    @pytest.fixture(scope="class")
    def map1(self):
        map1 = MagicMock()
        contacts = [
            ((0, 3), 2),
            ((1, 4), 2),
            ((2, 5), 2),
        ]
        map1.combine.return_value.__iter__.return_value = contacts
        return map1

    @pytest.fixture(scope="class")
    def map2(self):
        map2 = MagicMock()
        contacts = [
            ((0, 3), 1),
            ((1, 4), 1),
            ((2, 5), 1),
        ]
        map2.combine.return_value.__iter__.return_value = contacts
        return map2

    @pytest.fixture(scope="class")
    def empty_map(self):
        map1 = MagicMock()
        map1.combine.return_value.__iter__.return_value = []
        return map1

    def test_compare_3frames(self, frame1, frame2, frame3, map1):
        calc1 = RMSCLDCalculator(frame1, frame2, map1, map1)
        calc2 = RMSCLDCalculator(frame1, frame3, map1, map1)
        score1 = calc1.calculate()
        score2 = calc2.calculate()
        assert score2 < 1
        assert score1 > 1
        assert score1 == pytest.approx(3.0)

    def test_no_contacts(self, frame1, empty_map):
        with pytest.raises(ValueError):
            RMSCLDCalculator(frame1, frame1, empty_map, empty_map)

    def test_no_certain(self, frame1, frame2, map2):
        calc = RMSCLDCalculator(frame1, frame2, map2, map2)
        score = calc.calculate()
        assert score < 3.0


class TestSystem:
    def test_martini(self, trajectories_path, structures_dir):
        stc_path = Path(structures_dir) / "martini" / "gpcr_d.pdb"
        factory = BioPythonAtomSetFactory([MartiniResidue, MonoatomicIon])
        with open(stc_path) as fh:
            (topology,) = StructureLoader(atom_set_factory=factory).load_structures(
                [fh]
            )
        traj_path = trajectories_path.format(type="xtc", traj="gpcr_5frames")
        trajectory = TrajectoryLoader().load_trajectory(traj_path, topology)
        frame1 = freeze_frame(trajectory)
        trajectory.set_frame(4)
        frame2 = freeze_frame(trajectory)
        crit = PointsDistanceCriterion("BB", 33.0, 2.0)
        sele1 = AtomSetExactType(MartiniResidue) * ChainSelection("A")
        sele2 = AtomSetExactType(MartiniResidue) * ChainSelection("B")
        f1_cmc = ContactMapCalculator(frame1, crit, selections=(sele1, sele2))
        f2_cmc = ContactMapCalculator(frame2, crit, selections=(sele1, sele2))
        cm1 = f1_cmc.calculate_contact_map()
        cm2 = f2_cmc.calculate_contact_map()
        calc = RMSCLDCalculator(frame1, frame2, cm1, cm2, point="last_sc")
        score = calc.calculate()
        assert 0 < score < 5.0

    def test_full_atom(self, trajectories_path, topologies_path):
        traj_path = trajectories_path.format(type="xtc", traj="mdm2_5frames")
        topo_path = topologies_path.format(topo="mdm2")
        (topology,) = get_structures_from_file(topo_path)
        trajectory = TrajectoryLoader().load_trajectory(traj_path, topology)
        frame1 = freeze_frame(trajectory)
        trajectory.set_frame(4)
        frame2 = freeze_frame(trajectory)

        crit = get_default_protein_criterion()
        f1_cmc = ContactMapCalculator(frame1, crit)
        f2_cmc = ContactMapCalculator(frame2, crit)
        cm1 = f1_cmc.calculate_contact_map()
        cm2 = f2_cmc.calculate_contact_map()

        n_contacts = len(tuple(cm1.combine(cm2)))
        assert n_contacts > 0

        calc = RMSCLDCalculator(frame1, frame2, cm1, cm2)
        score = calc.calculate()
        assert 10 > score > 0
