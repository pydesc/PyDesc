import numpy

from pydesc.api.criteria import get_rc_distance_criterion
from pydesc.contacts.maps import ContactMapCalculator
from pydesc.contacts.maps import FrequencyMap


def calculate_contact_map(structure, criterion=None):
    """Calculate contact map for given structure.

    By default if uses rc distance criterion, unless other criterion was passed.

    Args:
        structure: any pydesc structure or sub structure.
        criterion: optional; any contact criterion. By default rc distance criterion
        is used.

    Returns:
        ContactMap: contact map instance.

    """
    if criterion is None:
        criterion = get_rc_distance_criterion()
    calculator = ContactMapCalculator(structure, criterion)
    contact_map = calculator.calculate_contact_map()

    return contact_map


def create_frequency_map_from_contact_maps(contact_maps):
    """Create frequency map from sequence of contact maps of the same structure,
    e.g. for NMR frames.

    Structure of first contact map is passed to newly created instance of FrequencyMap.
    Contacts of value 1 are taken into account with weight 0.5.

    Args:
        contact_maps: sequence of pydesc contact maps.

    Returns:
        FrequencyMap: map of contacts and their frequency in (pseudo)trajectory (any
        set of different frames of the same structure).

    """
    trajectory = contact_maps[0].converter
    n_frames = len(contact_maps)
    matrix = contact_maps[0].get_dok_matrix().astype(numpy.float64)
    for contact_map in contact_maps[1:]:
        matrix += contact_map.get_dok_matrix()
    matrix /= 2
    frequency_map = FrequencyMap(matrix, trajectory, n_frames)
    return frequency_map
