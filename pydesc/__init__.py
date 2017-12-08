# Copyright 2017 Pawel Daniluk
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

class PropertiesRecord(object):

    """Class that stores properties values."""

    def __init__(self, initial_annotations):
        """Record constructor.

        Argument:
        initial_annotations -- dict; keys are attrs names (strings), while values are any pyhon objects.
        """
        for name, value in initial_annotations.items():
            setattr(self, name, value)

    def __repr__(self):
        return "".join(["<Rec:"] + [name + " - " + str(value) + "; " for name, value in self.__dict__.items()] + [">"])
