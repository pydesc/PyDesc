# Copyright 2017 Tymoteusz Oleniecki
#
# This file is part of PyDesc.
#
# PyDesc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyDesc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.

"""
PyDesc Warnings and Exceptions.

Usage: to use Exceptions simply import this module and raise them.
To use configurable warnings import this module and use static method warn
in class Warn.
In configuration manager set filter values in branch warnings.class_filters.
<Warning Class> to one of values from list below:
error -- raises error instead of printing warning.
ignore -- skips warnings.
always -- always prints warnings.

Exception classes:
DiscontinuityError -- exceptions raised when chain discontinuity is found.
WrongAtomDistances -- exceptions raised when wrong distances between atoms
occurs in monomer.

Warnings classes:
CopyDownload -- information about creation of local copy of file from data
base.
Info -- any information.
DeprecationWarning -- warning against unsupported methods.
LocalCopyAccess -- information about access to local copy of file.
MonomerCreationFailed -- information about failure during monomer creation.
NoConfiguration -- information about lack of configuration in config manager.
UnknownParticleName -- information about occurrence of particle of unknown
name.
UserWarning -- default warning.

created: 04.02.2014 , Tymoteusz 'hert' Oleniecki
"""

import warnings

from pydesc.config import ConfigManager

# pylint: disable=no-member
ConfigManager.new_branch("warnings")
ConfigManager.warnings.set_default("quiet", False)
ConfigManager.warnings.new_branch("class_filters")
ConfigManager.warnings.class_filters.set_default(
    "IncompleteChainableParticle", "always"
)
ConfigManager.warnings.class_filters.set_default("Info", "always")
ConfigManager.warnings.class_filters.set_default("DeprecationWarning", "default")
ConfigManager.warnings.class_filters.set_default("NoConfiguration", "error")
ConfigManager.warnings.class_filters.set_default("UnknownParticleName", "always")
ConfigManager.warnings.class_filters.set_default("UserWarning", "always")
ConfigManager.warnings.class_filters.set_default("WrongMonomerType", "always")


# pylint: enable=no-member


def warn(warn_inst, stack_level=0):
    """Function that deals with warnings in a way set in configuration manager.

    Arguments:
    warn_inst -- instance of warning to be called.
    stack_level -- int, describes level of traceback to be printed.
    """
    if not ConfigManager.warnings.quiet:  # pylint: disable=no-member
        stack_level += 2
        warnings.warn(warn_inst, stacklevel=stack_level)


def set_filters():
    """Sets filters in filter context manager."""
    warnings.resetwarnings()

    # TODO: Scan for warnings and exceptions automatically. Add checks for
    #  warnings without entries in ConfigManager.
    for wrn_cls in (
        Info,
        DeprecationWarning,
        NoConfiguration,
        UnknownParticleName,
        UserWarning,
    ):
        warn_name = wrn_cls.__name__
        filter_ = getattr(ConfigManager.warnings.class_filters, warn_name)
        warnings.filterwarnings(filter_, category=wrn_cls)


class WarnManager(warnings.catch_warnings):
    """Context manager that handles warnings raised in different contexts.

    Used in pydesc.structure and pydesc.monomer modules to avoid raising
    warnings during creation of unused monomer instances.
    Subclass of warnings.catch_warnings class. Extends __init__ and
    __enter__ methods. Overrides __repr__, __call__ and __exit__ methods.
    """

    def __init__(self, obj=None):
        """WarnManager constructor.

        Argument:
        obj -- an object that is to be tested in different contexts. Used
        only for recognition purpose. None by default.
        """
        self.obj = obj
        self.exceptions = {None: [], "__internal__": []}
        self.last_context = None
        warnings.catch_warnings.__init__(self, record=True)

    def __repr__(self):
        return "<WarnManager for %s object>" % str(self.obj)

    def __call__(self, context):
        """Gets new context and returns self.

        Argument:
        context -- any object (e.g. string) that specifies context.
        """
        self.last_context = context
        self.exceptions[context] = []
        return self

    def __enter__(self):
        """Extended super class method called at the beginning of 'with'
        statement code block.

        This method provides possibility of catching warnings as
        warnings.catch_warnings do.
        """
        self._entered = False
        self.exceptions["__internal__"] = warnings.catch_warnings.__enter__(self)
        return self

    def __exit__(self, warning_type, warning_instance, traceback, *args, **kwargs):
        """Overridden super class method called in 'with' statement when
        warning or exception was raised in code block.

        Arguments:
        warning_type -- type of raised warning or exception.
        warning_value -- value of raised warning or exception.
        traceback -- raised exception traceback.

        Exceptions that are not Warning subclass are raised immediately,
        while warnings are stored for future usage.
        """
        self.exceptions[self.last_context].extend(
            [
                self.exceptions["__internal__"].pop(0).message
                for _ in self.exceptions["__internal__"]
            ]
        )
        warnings.catch_warnings.__exit__(self)
        if warning_type is None:
            return True
            # that means there was no error
        if not issubclass(warning_type, Warning):
            return False
            # error was not a warning and it should be raised
        self.exceptions[self.last_context].append(warning_instance)
        return True
        # error was just a warning and we want it to be added to the list

    def raise_all(self, context=None):
        """Method that raises all stored warnings.

        Argument:
        context -- any object (e.g. string) that specifies context of
        warnings to be raise.

        Method use pydesc.warnexcept.warn function to throw warnings,
        therefore pydesc configuration affects warnings filtering.
        """
        for warning in self.exceptions.get(context, []):
            warn(warning, 4)


class CannotCalculateContact(Exception):
    """Class of exceptions raised by contact criteria whenever given mers
    lack attributes or properties needed to calculate contact."""


class DiscontinuityError(Exception):
    """Class of exceptions raised due to discontinuity between chainable mers."""

    pass


class FrameNotAvailable(Exception):
    """Raised when trajectory frame to be set exceeds available range."""

    pass


class IncompleteParticle(Exception):
    """Class of exceptions raised when incomplete particle is given to
    create monomer instance."""


class UnknownPDBid(Exception):
    """Unknown PDB id error."""

    pass


class WrongAtomDistances(Exception):
    """Class of exceptions raised due to wrong distance between atoms
    occurring in structures."""

    pass


class WrongElement(Exception):
    """Error shown when given element is incorrect."""

    pass


class WrongAtomSetType(Exception):
    """Class of warnings. Warning against wrong type of mers given for
    contact calculation under specific criteria. Printed by default.
    """

    pass


class NotASlice(Exception):
    pass


class Info(Warning):
    """Class of warnings. Standard information. Printed by default."""


class NoConfiguration(Warning):
    """Class of warnings given while no appropriate configuration in
    configuration manager is found. Raised as an error by default.
    """


class UnknownParticleName(Warning):
    """Class of warnings. Information about particles that have unknown
    names. Printed by default.
    """


class InvalidID(Exception):
    """Database handler exception for wrong structure ID."""

    pass


class OperationModeError(Exception):
    """Database handler exception for not served mode of getting files."""

    pass
