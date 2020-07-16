from pydesc.contacts.base import ContactsAlternative
from pydesc.contacts.base import ContactsConjunction
from pydesc.contacts.geometrical import DistancesDifferenceCriterion
from pydesc.contacts.geometrical import PointsDistanceCriterion
from pydesc.mers.full_atom import Residue
from pydesc.selection import MerSubclasses

RESIDUE_SELECTION = MerSubclasses(Residue)

CA_DISTANCE = PointsDistanceCriterion("ca", 6.0, 0.5)
CA_DISTANCE.set_selections(RESIDUE_SELECTION, RESIDUE_SELECTION)
CBX_DISTANCE = PointsDistanceCriterion("cbx", 6.5, 0.5)
CBX_DISTANCE.set_selections(RESIDUE_SELECTION, RESIDUE_SELECTION)
RC_DISTANCE = PointsDistanceCriterion("rc", 7.5, 0.5)

CA_CBX_DISTANCE_DIFFERENCE = DistancesDifferenceCriterion("ca", "cbx", 0.75, 0.05)

DEFAULT_PROTEIN = ContactsAlternative(
    CA_DISTANCE, ContactsConjunction(CBX_DISTANCE, CA_CBX_DISTANCE_DIFFERENCE)
)
