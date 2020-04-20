# Copyright 2017 Agnieszka Mykowiecka, Tymoteusz Oleniecki
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
Classes that let PyDesc deal with coordinates calculations.

created: 10.07.2013, Agnieszka Mykowiecka, Tymoteusz 'hert' Oleniecki
"""

import numpy
import scipy.linalg

# pylint: disable=no-member, invalid-name
norm = scipy.linalg.get_blas_funcs("nrm2")
axpy = scipy.linalg.get_blas_funcs("axpy")


# pylint: enable=no-member, invalid-name


class TRTMatrix:
    """Translation-Rotation-Translation matrices class.

    Stores prerotational translation vector, rotation matrix and
    post-rotational translation vector.
    Attributes:
    pre_vector -- list of three floats.
    rotation_matrix -- list of three lists containing three floats.
    post_vector -- post-rotational translation vector, list of three
    floats.
    """

    # pylint: disable=no-member
    def __init__(self):
        """Translation-Rotation-Translation matrix objects' constructor.

        Takes no arguments. Sets three attributes:
        pre_vector -- numpy.array of three floats; by
        default [.0, .0, .0].
        post_vector -- numpy array as above. Represents
        post-rotational vector.
        rotation_matrix -- numpy.array. diagonal identity matrix 3x3 by
        default.
        """
        self.pre_vector = numpy.array([0.0, 0.0, 0.0])
        self.rotation_matrix = numpy.array(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        )
        self.post_vector = numpy.array([0.0, 0.0, 0.0])

    def add_rotation(self, rotation_matrix):
        """Adds rotation to stored rotation matrix.

        Argument:
        rotation_matrix -- numpy.array of floats 3x3 shape.
        """
        rotation_matrix = numpy.array(rotation_matrix).reshape(3, 3)
        self.rotation_matrix = rotation_matrix * self.rotation_matrix

    def add_translation(self, vector):
        """Adds vector to post-rotational translation vector.

        Argument:
        vector -- list of three floats.
        """
        vector = numpy.array(vector)
        self.post_vector += vector

    def add_prerotational_translation(self, vector):
        """Adds vector to prerotational translation vector.

        Argument:
        vector -- list of three floats.
        """
        vector = numpy.array(vector)
        self.pre_vector += vector

    def transform(self, x=None, y=None, z=None, vec=None):
        """Returns transformed coordinates.

        Moves given point by pre_vector, rotates the
        result in accordance with rotation_matrix
        and in the end moves the point by post_vector.

        Arguments:
        x, y, z -- spacial coordinates (optional).
        vec -- numpy.array (optional).

        Needs vec or set of coord to proceed.
        """
        if vec is None:
            vec = numpy.array((x, y, z))
        pre_vector = self.pre_vector + vec
        rotated_vector = numpy.squeeze(self.rotation_matrix.dot(pre_vector).A)
        return rotated_vector + self.post_vector

    # pylint: enable=invalid-name
    # names x, y and z are the best possible here

    def reset_translation(self):
        """Resets translation vector to zero vector."""
        self.post_vector = numpy.array([0.0, 0.0, 0.0])

    def reset_rotation(self):
        """Resets rotation matrix to diagonal identity matrix."""
        self.rotation_matrix = numpy.array(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        )

    def reset_prerotational_translation(self):
        """Resets prerotational translation vector to zero vector."""
        self.pre_vector = numpy.array([0.0, 0.0, 0.0])

    def reset(self):
        """Resets rotation, translation and prerotational translation to
        default values."""
        self.reset_translation()
        self.reset_rotation()
        self.reset_prerotational_translation()

    # pylint: enable=no-member
    # matrix and array are present among numpy attrs

    def combine(self, other):
        ret = TRTMatrix()
        rotation = numpy.dot(other.rotation_matrix, self.rotation_matrix)
        ret.add_rotation(rotation)
        ret.add_translation(self.post_vector + other.post_vector)
        ret.add_prerotational_translation(self.pre_vector + other.pre_vector)
        return ret


class Coord:
    """Stores coordinates as position vector.

    Note, that addition and subtraction of coord instances can be performed.
    It is also possible to multiply each coordinate by number.
    """

    def __init__(self, x=0.0, y=0.0, z=0.0, numpy_vec=None):
        """Coord constructor.

        Arguments:
        x, y, z -- spatial coordinates.
        numpy_vec -- NumPy array containing coordinates (if provided XYZ are
        ignored)
        """
        if numpy_vec is not None:
            self.vector = numpy_vec
        else:
            self.vector = numpy.array(
                (float(x), float(y), float(z))
            )  # pylint: disable=no-member

    def __add__(self, coord_obj):
        """Returns result of addition of two vectors represented by Coord
        instances as new Coord instance.

        Argument:
        coord_obj -- instance of Coord class.
        """
        return Coord(numpy_vec=(self.vector + coord_obj.vector))

    def __sub__(self, coord_obj):
        """Returns result of subtraction of two vectors represented by Coord
        instances as new Coord instance.

        Argument:
        coord_obj -- instance of Coord class.
        """
        return Coord(numpy_vec=(self.vector - coord_obj.vector))

    def __mul__(self, factor):
        """Returns result of multiplication of vector a by given number as
        new Coord instance.

        Argument:
        factor -- a float or an integer.
        """
        return Coord(numpy_vec=(self.vector * factor))

    def __div__(self, denominator):
        """Returns result of division of vector by a given number as new
        Coord instance.

        Argument:
        denominator -- a float or an integer.
        """
        return Coord(numpy_vec=(self.vector / denominator))

    def __iter__(self):
        """Returns iterator that iterates over stored coordinates: x, y, z."""
        return iter(self.vector)

    def __repr__(self):
        return "<Coord: %f %f %f>" % tuple(self.vector)

    def __abs__(self):
        """Returns Coord's vector length."""
        return self.calculate_length()

    def copy(self):
        """Return copy of self."""
        coord_type = type(self)
        return coord_type(numpy_vec=self.vector)

    def calculate_length(self):
        """Calculates the length of the position vector.

        Returns the distance in a given unit of length.
        """
        return norm(self.vector)

    def extend(self, length):
        """Returns position vector extended by given length.

        Argument:
        length -- float in position vector unit.
        """
        return self * (1 + length / self.calculate_length())

    def get_coord(self, transformation_matrix=None):
        """Returns a list of stored spacial coordinates.

        Argument:
        transformation_matrix -- instance of TRTMatrix to be used in
        coordinates transformation.
        """
        if transformation_matrix is not None:
            return tuple(transformation_matrix.transform(*self.vector))
        return tuple(self.vector)

    def get_unit_vector(self):
        return self.vector / self.calculate_length()

    # pylint: disable=invalid-name
    @property
    def x(self):
        """Returns x coord."""
        return self.vector[0]

    @x.setter
    def x(self, value):
        """Sets x coord."""
        self.vector[0] = value

    @property
    def y(self):
        """Returns y coord."""
        return self.vector[1]

    @y.setter
    def y(self, value):
        """Sets y coord."""
        self.vector[1] = value

    @property
    def z(self):
        """Returns z coord."""
        return self.vector[2]

    @z.setter
    def z(self, value):
        """Sets z coord."""
        self.vector[2] = value


# pylint: enable=invalid-name
# names x, y and z are the best possible here


class Plane:
    """Stores coordinates as position vector.

    Methods:
    ort_projection -- returns object of class Coord with coordinates of the
    orthogonal projection of given point onto plane.
    dihedral_angle_cos -- returns cosines of the angle between two planes.
    bisection_plane -- returns the bisection plane for two given planes.
    """

    @staticmethod
    def build(v1, v2, v3):
        vec = numpy.cross((v1.vector - v2.vector), (v3.vector - v2.vector))

        d = -numpy.dot(vec, v1.vector)
        plane = Plane(numpy_vec=vec, d=d)

        return plane

    # pylint: disable=invalid-name
    def __init__(self, a=0, b=0, c=0, d=0, numpy_vec=None):
        """Plane constructor.

        Normal vector is normalized during initialization.

        Arguments:
        a, b, c, d -- plane coefficients.
        """

        if numpy_vec is not None:
            vec = numpy_vec
        else:
            vec = numpy.array((float(a), float(b), float(c)))

        nrm = norm(vec)

        self.norm_vec = vec / nrm
        self.d = float(d) / nrm

    def set_perpendicular_vector(self, x=0.0, y=0.0, z=0.0, vec=None):
        """
        """
        pass

    @property
    def a(self):
        return self.norm_vec[0]

    @property
    def b(self):
        return self.norm_vec[1]

    @property
    def c(self):
        return self.norm_vec[2]

    @property
    def perpendicular_x(self):
        return self.norm_vec[0]

    @property
    def perpendicular_y(self):
        return self.norm_vec[1]

    @property
    def perpendicular_z(self):
        return self.norm_vec[2]

    @property
    def perpendicular_vec(self):
        return self.norm_vec

    def ort_projection(self, point):
        """Returns coordinates of the orthogonal projection of given point onto
        plane.

        Argument:
        point -- object of class Coord.
        """

        r = numpy.dot(self.perpendicular_vec, point.vector) + self.d

        res = point.vector - self.norm_vec * r

        return Coord(numpy_vec=res)

    def dihedral_angle(self, plane2):
        """
        """
        nrm = 1
        cos = numpy.dot(self.perpendicular_vec, plane2.perpendicular_vec) / nrm
        c_prod = numpy.cross(self.perpendicular_vec, plane2.perpendicular_vec)
        c_p_norm = norm(c_prod)
        sin = c_p_norm / nrm
        angle = numpy.arctan2(sin, cos)

        return angle

    def dihedral_angle_cos(self, plane2):
        """Returns value of the cosinus of the angle between two planes.

        Argument:
        plane2 -- object of class Plane.
        """
        return numpy.dot(self.perpendicular_vec, plane2.perpendicular_vec)

    def bisection_plane(self, plane2):
        """Returns coefficients of bisection plane.

        Argument:
        plane2 -- object of class Plane.
        """
        v_prim, d_prim = self.norm_vec, self.d
        v1prim, d1prim = plane2.norm_vec, plane2.d

        cos3 = numpy.dot(v_prim, v1prim)

        if cos3 < 0:
            v = v_prim - v1prim
            d = d_prim - d1prim
        else:
            v = v_prim + v1prim
            d = d_prim + d1prim

        res = Plane(numpy_vec=v, d=d)

        return res


# pylint: enable=invalid-name
# arguments and local names are informative here
