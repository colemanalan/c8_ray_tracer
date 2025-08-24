import unittest
from c8_tracer_py import Vec2, Vec3


class TestVec2(unittest.TestCase):
    def test_construction_and_repr(self):
        v = Vec2(1.0, 2.0)
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 2.0)
        repr(v)

    def test_equality_and_negation(self):
        v1 = Vec2(3.0, -4.0)
        v2 = Vec2(3.0, -4.0)
        self.assertEqual(v1, v2)
        self.assertEqual(-v1, Vec2(-3.0, 4.0))

    def test_indexing_and_iteration(self):
        v = Vec2(5.0, 6.0)
        self.assertEqual(v[0], 5.0)
        self.assertEqual(v[1], 6.0)
        with self.assertRaises(IndexError):
            _ = v[2]
        v[0] = 7.0
        self.assertEqual(v.x, 7.0)
        self.assertEqual(len(v), 2)
        self.assertEqual(list(v), [7.0, 6.0])

    def test_arithmetic(self):
        v1 = Vec2(1.0, 2.0)
        v2 = Vec2(3.0, 4.0)
        self.assertEqual(v1 + v2, Vec2(4.0, 6.0))
        self.assertEqual(v2 - v1, Vec2(2.0, 2.0))
        self.assertEqual(v1 * 2.0, Vec2(2.0, 4.0))
        self.assertEqual(v2 / 2.0, Vec2(1.5, 2.0))

    def test_dot_length_norm(self):
        v = Vec2(3.0, 4.0)
        self.assertAlmostEqual(v.dot(v), 25.0)
        self.assertAlmostEqual(v.length(), 5.0)
        self.assertAlmostEqual(v.norm(), 5.0)
        self.assertEqual(v.normalized(), Vec2(0.6, 0.8))
        self.assertEqual(v.normalized(), Vec2(3.0, 4.0) / v.norm())


class TestVec3(unittest.TestCase):
    def test_construction_and_repr(self):
        v = Vec3(1.0, 2.0, 3.0)
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 2.0)
        self.assertEqual(v.z, 3.0)
        repr(v)

    def test_equality_and_negation(self):
        v1 = Vec3(1.0, -2.0, 3.0)
        v2 = Vec3(1.0, -2.0, 3.0)
        self.assertEqual(v1, v2)
        self.assertEqual(-v1, Vec3(-1.0, 2.0, -3.0))

    def test_indexing_and_iteration(self):
        v = Vec3(7.0, 8.0, 9.0)
        self.assertEqual(v[0], 7.0)
        self.assertEqual(v[1], 8.0)
        self.assertEqual(v[2], 9.0)
        with self.assertRaises(IndexError):
            _ = v[3]
        v[2] = 10.0
        self.assertEqual(v.z, 10.0)
        self.assertEqual(len(v), 3)
        self.assertEqual(list(v), [7.0, 8.0, 10.0])

    def test_arithmetic(self):
        v1 = Vec3(1.0, 2.0, 3.0)
        v2 = Vec3(4.0, 5.0, 6.0)
        self.assertEqual(v1 + v2, Vec3(5.0, 7.0, 9.0))
        self.assertEqual(v2 - v1, Vec3(3.0, 3.0, 3.0))
        self.assertEqual(v1 * 2.0, Vec3(2.0, 4.0, 6.0))
        self.assertEqual(v2 / 2.0, Vec3(2.0, 2.5, 3.0))

    def test_dot_cross_length_norm(self):
        v1 = Vec3(1.0, 0.0, 0.0)
        v2 = Vec3(0.0, 1.0, 0.0)
        cross = v1.cross(v2)
        self.assertEqual(cross, Vec3(0.0, 0.0, 1.0))
        self.assertAlmostEqual(v1.dot(v2), 0.0)
        v = Vec3(3.0, 4.0, 0.0)
        self.assertAlmostEqual(v.length(), 5.0)
        self.assertAlmostEqual(v.norm(), 5.0)
        self.assertEqual(v.normalized(), Vec3(0.6, 0.8, 0.0))


if __name__ == "__main__":
    unittest.main()
