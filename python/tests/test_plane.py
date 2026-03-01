import unittest

from c8_tracer.c8_tracer_ext import Vec3
from c8_tracer.c8_tracer_ext import Plane


class TestPlane(unittest.TestCase):
    def test_constructor(self):
        center = Vec3(0, 0, 0)
        normal = Vec3(0, 0, 1)
        plane = Plane(center, normal)

        self.assertEqual(normal, plane.getNormal())
        self.assertEqual(center, plane.getCenter())


if __name__ == "__main__":
    unittest.main()
