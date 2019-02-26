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
Unit tests_legacy for geometry.py.

Usage:
    python geometry_test.py [-v]

    or

    python -m unittest [-v] geometry_test
"""

import tests_legacy.syntax_check as syntax_check
from tests_legacy.syntax_check import testing, notest
import unittest

import random
import string
import math
import operator
import itertools


import pydesc.geometry as geometry

syntax_check.module = geometry

TestSyntax = syntax_check.module_syntax()

randomcoord = lambda n=3: [random.uniform(-5, 5) for i in range(n)]
randomword = lambda length: ''.join(random.choice(string.lowercase) for i in range(length))
randomangle = lambda: random.uniform(-math.pi, math.pi)

# pylint: disable=C0111

notest(geometry.Coord.get_coord)
notest(geometry.Coord.__abs__)
notest(geometry.TRTMatrix.reset)
notest(geometry.TRTMatrix.reset_prerotational_translation)
notest(geometry.TRTMatrix.reset_rotation)
notest(geometry.TRTMatrix.reset_translation)


def rotmat(theta, axis):
    c = math.cos(theta)
    s = math.sin(theta)

    if axis == 0:
        res = [[1, 0, 0], [0, c, s], [0, -s, c]]
    elif axis == 1:
        res = [[c, 0, s], [0, 1, 0], [-s, 0, c]]
    else:
        res = [[c, s, 0], [-s, c, 0], [0, 0, 1]]

    return res


dot = lambda v1, v2: sum(map(operator.mul, v1, v2))


class CoordBasic(unittest.TestCase):

    @testing(geometry.Coord)
    @testing(geometry.Coord.__init__)
    @testing(geometry.Coord.x)
    @testing(geometry.Coord.y)
    @testing(geometry.Coord.z)
    def test_init(self):
        x, y, z = randomcoord()

        c1 = geometry.Coord(x, y, z)
        self.assertEquals((x, y, z), (c1.x, c1.y, c1.z))

        c2 = geometry.Coord(x=x, y=y, z=z)
        self.assertEquals((x, y, z), (c2.x, c2.y, c2.z))


class CoordTest(unittest.TestCase):

    def setUp(self):
        self.c1 = geometry.Coord(*randomcoord())
        self.c2 = geometry.Coord(*randomcoord())

    @testing(geometry.Coord.__add__)
    def test_add(self):
        c1 = self.c1
        c2 = self.c2
        for i in 'xyz':
            self.assertAlmostEqual((getattr(c1, i) + getattr(c2, i)), getattr(c1 + c2, i))

    @testing(geometry.Coord.__sub__)
    def test_sub(self):
        c1 = self.c1
        c2 = self.c2
        for i in 'xyz':
            self.assertAlmostEqual((getattr(c1, i) - getattr(c2, i)), getattr(c1 - c2, i))

    @testing(geometry.Coord.__add__)
    @testing(geometry.Coord.__sub__)
    def test_addsub(self):
        c1 = self.c1
        c2 = self.c2
        for i in 'xyz':
            self.assertAlmostEqual(getattr(c1, i), getattr(c1 - c2 + c2, i))

    @testing(geometry.Coord.__mul__)
    def test_mul(self):
        c1 = self.c1
        f = random.uniform(-5, 5)
        for i in 'xyz':
            self.assertAlmostEqual(getattr(c1, i) * f, getattr(c1 * f, i))

    @testing(geometry.Coord.__div__)
    def test_div(self):
        c1 = self.c1
        f = random.uniform(-5, 5)
        for i in 'xyz':
            self.assertAlmostEqual(getattr(c1, i) / f, getattr(c1 / f, i))

    @testing(geometry.Coord.__iter__)
    def test_iter(self):
        c1 = self.c1
        self.assertEqual([c1.x, c1.y, c1.z], list(c1))

    @testing(geometry.Coord.calculate_length)
    def test_length(self):
        c1 = self.c1
        mod = math.sqrt(dot(c1, c1))
        self.assertAlmostEqual(c1.calculate_length(), mod)

    @testing(geometry.Coord.extend)
    def test_extend(self):
        c1 = self.c1

        ff = random.uniform(-5, 5)
        mod = c1.calculate_length()

        for f in [ff, -ff, 0]:
            c = c1.extend(mod * f)
            for i in 'xyz':
                self.assertAlmostEqual(getattr(c1, i) * (f + 1), getattr(c, i))


class CoordEnhancements(unittest.TestCase):

    def setUp(self):
        self.c1 = geometry.Coord(*randomcoord())
        self.c2 = geometry.Coord(*randomcoord())

    @unittest.expectedFailure
    def test_abs(self):
        c = self.c1

        self.assertTrue(hasattr(c, '__abs__'))
        self.assertAlmostEqual(abs(c), c.compute_length())

    @unittest.expectedFailure
    def test_dot(self):
        self.assertTrue(hasattr(self.c1, 'dot'))

        cc = self.c1.dot(self.c2)
        self.assertAlmostEqual(cc, dot(self.c1, self.c2))

    @unittest.expectedFailure
    def test_deprec(self):
        c = self.c1
        syntax_check.test_deprec(self, c.get_coord, "get_coord")
        syntax_check.test_deprec(self, c.get_transformed_coord, "get_transformed_coord")
        syntax_check.test_deprec(self, c.calculate_length, "calculate_length")


class PlaneBasic(unittest.TestCase):

    @testing(geometry.Plane)
    @testing(geometry.Plane.__init__)
    def test_init(self):
        a, b, c, d = randomcoord(4)

        p = geometry.Plane(a, b, c, d)
        self.assertEquals((a, b, c, d), (p.a, p.b, p.c, p.d))


class PlaneTest(unittest.TestCase):

    def setUp(self):
        rc = randomcoord(4)
        self.p1 = geometry.Plane(*rc)
        self.np1 = geometry.Plane(*[-x for x in rc])
        rc = randomcoord(4)
        self.p2 = geometry.Plane(*rc)
        self.np2 = geometry.Plane(*[-x for x in rc])
        self.c1 = geometry.Coord(*randomcoord())
        self.c2 = geometry.Coord(*randomcoord())

    @testing(geometry.Plane.ort_projection)
    def test_ort_projection(self):
        cp = self.p1.ort_projection(self.c1)
        cnp = self.np1.ort_projection(self.c1)

        self.assertAlmostEqual(list(cp), list(cnp))

    @testing(geometry.Plane.ort_projection)
    def test_ort_projection1(self):
        self.p1.d = 0
        self.np1.d = 0
        cp = self.p1.ort_projection(self.c1)
        cpp = self.c1 - cp

        self.assertAlmostEqual(dot(cp, cpp), 0)
        cnp = self.np1.ort_projection(self.c1)
        cpp = self.c1 - cnp
        self.assertAlmostEqual(dot(cnp, cpp), 0)

        self.assertAlmostEqual(list(cp), list(cnp))

    @testing(geometry.Plane.bisection_plane)
    @testing(geometry.Plane.dihedral_angle_cos)
    def test_bisection_plane(self):
        bpl = []
        bpl.append(self.p1.bisection_plane(self.p2))
        bpl.append(self.p2.bisection_plane(self.p1))
        bpl.append(self.np1.bisection_plane(self.p2))
        bpl.append(self.np2.bisection_plane(self.p1))
        bpl.append(self.p1.bisection_plane(self.np2))
        bpl.append(self.p2.bisection_plane(self.np1))
        bpl.append(self.np1.bisection_plane(self.np2))
        bpl.append(self.np2.bisection_plane(self.np1))

        for bp in bpl[1:]:
            self.assertAlmostEqual(abs(bpl[0].dihedral_angle_cos(bp)), 1)


class TRTMatrixTest(unittest.TestCase):
    def setUp(self):
        self.l = [(randomangle(), random.choice(range(3))) for dummy in range(10)]

    @testing(geometry.TRTMatrix)
    @testing(geometry.TRTMatrix.__init__)
    def test_init(self):
        geometry.TRTMatrix()

    @testing(geometry.TRTMatrix.add_translation)
    @testing(geometry.TRTMatrix.add_prerotational_translation)
    @testing(geometry.TRTMatrix.add_rotation)
    def test_basicops(self):
        mat = geometry.TRTMatrix()
        mat.add_translation(randomcoord())
        mat.add_prerotational_translation(randomcoord())
        mat.add_rotation(rotmat(random.uniform(-math.pi, math.pi), random.choice(range(3))))

    @testing(geometry.TRTMatrix.add_rotation)
    def test_rot(self):
        one = geometry.TRTMatrix()
        mat = geometry.TRTMatrix()
        for (theta, axis) in self.l:
            mat.add_rotation(rotmat(theta, axis))
        for (theta, axis) in reversed(self.l):
            mat.add_rotation(rotmat(-theta, axis))
        map(self.assertAlmostEqual, itertools.chain(*mat.rotation_matrix.tolist()), itertools.chain(*one.rotation_matrix.tolist()))

    @testing(geometry.TRTMatrix.add_rotation)
    def test_transpose(self):
        one = geometry.TRTMatrix()
        mat = geometry.TRTMatrix()
        for (theta, axis) in self.l:
            mat.add_rotation(rotmat(theta, axis))
        trmat = geometry.TRTMatrix().rotation_matrix
        for i in range(3):
            for j in range(3):
                trmat[i, j] = mat.rotation_matrix[j, i]
        mat.add_rotation(trmat)
        map(self.assertAlmostEqual, itertools.chain(*mat.rotation_matrix.tolist()), itertools.chain(*one.rotation_matrix.tolist()))

    @testing(geometry.TRTMatrix.transform)
    def test_transform(self):
        for (theta, axis) in self.l:
            mat = geometry.TRTMatrix()
            mat.add_rotation(rotmat(theta, axis))
            vec = [0, 0, 0]
            vec[axis] = randomcoord(1)[0]
            map(self.assertAlmostEqual, mat.transform(0, 0, 0), (0, 0, 0))
            map(self.assertAlmostEqual, mat.transform(*vec), vec)

    @testing(geometry.TRTMatrix.add_rotation)
    @testing(geometry.TRTMatrix.transform)
    def test_linearity(self):
        mat = geometry.TRTMatrix()
        for (theta, axis) in self.l:
            mat.add_rotation(rotmat(theta, axis))
        scale = lambda vec, f: [x * f for x in vec]
        for dummy in range(100):
            vec = randomcoord()
            for f in randomcoord(100):
                map(self.assertAlmostEqual, mat.transform(*scale(vec, f)), scale(mat.transform(*vec), f))

    @testing(geometry.TRTMatrix.add_rotation)
    @testing(geometry.TRTMatrix.transform)
    def test_dotprod(self):
        mat = geometry.TRTMatrix()
        for (theta, axis) in self.l:
            mat.add_rotation(rotmat(theta, axis))
        for dummy in range(100):
            v1 = randomcoord()
            v2 = randomcoord()
            self.assertAlmostEqual(dot(mat.transform(*v1), mat.transform(*v2)), dot(v1, v2))

    @testing(geometry.TRTMatrix.add_rotation)
    @testing(geometry.TRTMatrix.add_translation)
    @testing(geometry.TRTMatrix.add_prerotational_translation)
    @testing(geometry.TRTMatrix.transform)
    def test_trans(self):
        mat1 = geometry.TRTMatrix()
        mat2 = geometry.TRTMatrix()
        mat = geometry.TRTMatrix()
        for (theta, axis) in self.l:
            mat.add_rotation(rotmat(theta, axis))
            mat1.add_rotation(rotmat(theta, axis))
            mat2.add_rotation(rotmat(theta, axis))
        for dummy in range(100):
            v1 = randomcoord()
            v2 = mat.transform(*v1)
            if random.choice([True, False]):
                mat1.add_prerotational_translation(v1)
                mat2.add_translation(v2)
            else:
                mat2.add_prerotational_translation(v1)
                mat1.add_translation(v2)
            for dummy in range(100):
                vec = randomcoord()
                map(self.assertAlmostEqual, mat1.transform(*vec), mat2.transform(*vec))

if __name__ == '__main__':
    unittest.main()
