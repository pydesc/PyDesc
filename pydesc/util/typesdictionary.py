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


class TypesDictionary(dict):
    """
    Dictionary which returns value for a given key or a value for a key which is a name of a super class of a given key.
    """

    def __getitem__(self, key):
        """
        Overrides dictionary __getitem__() method.
        """
        if key in self:
            return dict.__getitem__(self, key)
        else:
            try:
                return dict.__getitem__(self, key.__bases__[0])
            except:
                raise KeyError(key)
