from pydesc.contacts.maps import ContactMapCalculator
from pydesc.api.criteria import get_rc_distance_criterion


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
