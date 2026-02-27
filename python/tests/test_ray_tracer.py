import unittest

from c8_tracer.c8_tracer_ext import environment
from c8_tracer.c8_tracer_ext import Plane
from c8_tracer.c8_tracer_ext import Vec3
from c8_tracer.c8_tracer_ext.tracer import RayTracer2D

from c8_tracer import logging

logging.logger_tracer.set_level(logging.LogLevel.ERROR)


class TestNRMCInterface(unittest.TestCase):
    def setUp(self):
        n_deep = 1.78
        delta_n = 0.423
        length_scale = 77.0
        up = Vec3(0, 0, 1)
        ref_point = Vec3(0, 0, 0)

        self.env_flat = environment.IsotropicEnvironment(n_deep + delta_n)
        self.env_expo = environment.CartesianSingleExponentialEnvironment(
            n_deep,
            delta_n,
            length_scale,
            up,
            ref_point,
        )

        self.axis_of_symmetry = up
        self.perp_dir = (Vec3(2.3, 1.1, 0.0).cross(up)).normalized()
        self.min_step = 0.00001
        self.max_step = 10.0
        self.tolerance = 1e-6
        self.nRays = 21

        self.dz_plane = 40.0
        self.dr_pos = 50.0

        self.start = Vec3(0, 0, 0)  # start (0,0,0)
        self.target_even = self.start + self.perp_dir * self.dr_pos
        self.target_general = (
            self.start + self.perp_dir * self.dr_pos + up * 0.5 * self.dz_plane
        )
        self.plane_height = self.start + up * self.dz_plane
        self.target_across_plane = (
            self.start + self.perp_dir * self.dr_pos + up * 1.2 * self.dz_plane
        )

    def test_good_setup(self):
        dot = self.axis_of_symmetry.dot(self.perp_dir)
        self.assertAlmostEqual(0.0, dot)

        def height(x):
            return self.axis_of_symmetry.dot(x)

        def range(a, b):
            dpos = a - b
            dr = dpos - self.axis_of_symmetry * dpos.dot(self.axis_of_symmetry)
            return dr.norm()

        self.assertGreater(self.plane_height.z, self.start.z)  # plane is above start
        self.assertGreater(
            height(self.plane_height), height(self.start)
        )  # plane is above
        self.assertGreater(height(self.plane_height), height(self.target_even))
        self.assertGreater(height(self.plane_height), height(self.target_general))
        self.assertGreater(height(self.target_across_plane), height(self.plane_height))
        self.assertAlmostEqual(height(self.start), height(self.target_even))

        self.assertAlmostEqual(range(self.start, self.target_even), self.dr_pos)
        self.assertAlmostEqual(range(self.start, self.target_general), self.dr_pos)
        self.assertAlmostEqual(range(self.start, self.target_across_plane), self.dr_pos)

    def get_tracer(self, with_mirror=False):
        tracer = RayTracer2D(
            self.axis_of_symmetry,
            self.min_step,
            self.max_step,
            self.tolerance,
            self.nRays,
        )

        if with_mirror:
            tracer.AddReflectionLayer(Plane(self.plane_height, self.axis_of_symmetry))

        return tracer

    def test_RayTracer2D(self):
        self.get_tracer()

        RayTracer2D(
            self.axis_of_symmetry * -1.0,
            self.min_step,
            self.max_step,
            self.tolerance,
            self.nRays,
        )

        # same step size
        RayTracer2D(
            self.axis_of_symmetry,
            self.max_step,
            self.max_step,
            self.tolerance,
            self.nRays,
        )

        with self.assertRaises(RuntimeError):
            # invert step size
            RayTracer2D(
                self.axis_of_symmetry,
                self.max_step,
                self.min_step,
                self.tolerance,
                self.nRays,
            )

        for n_rays_bad in range(-1, 2, 1):
            with self.assertRaises(RuntimeError):
                # invalid steps
                RayTracer2D(
                    self.axis_of_symmetry,
                    self.min_step,
                    self.max_step,
                    self.tolerance,
                    n_rays_bad,
                )

    def test_AddReflectionLayer(self):
        tracer = self.get_tracer()
        tracer.AddReflectionLayer(Plane(Vec3(0, 0, 0), self.axis_of_symmetry))

        with self.assertRaises(ValueError):  # not perp to grad
            tracer.AddReflectionLayer(Plane(Vec3(0, 0, 0), self.perp_dir))

    def test_PropToPoint(self):
        tracer = self.get_tracer(True)

        # uniform environment
        paths_even = tracer.PropagateToPoint(
            self.start, self.target_even, self.env_flat
        )
        self.assertEqual(len(paths_even), 2)
        self.assertNotAlmostEqual(paths_even[0].emit.dot(paths_even[1].emit), 1.0)

        paths_general = tracer.PropagateToPoint(
            self.start, self.target_general, self.env_flat
        )
        self.assertEqual(len(paths_general), 2)
        self.assertNotAlmostEqual(paths_general[0].emit.dot(paths_general[1].emit), 1.0)

        paths = tracer.PropagateToPoint(
            self.start, self.target_across_plane, self.env_flat
        )
        self.assertEqual(
            len(paths), 1
        )  # nothing to change direction to allow two for transmission

        # exponential environment
        paths = tracer.PropagateToPoint(self.start, self.target_even, self.env_expo)
        self.assertEqual(len(paths), 2)

        paths = tracer.PropagateToPoint(self.start, self.target_general, self.env_expo)
        self.assertEqual(len(paths), 2)
        paths = tracer.PropagateToPoint(
            self.start, self.target_across_plane, self.env_expo
        )
        self.assertEqual(len(paths), 2)

    def test_GetSignalPaths(self):
        tracer = self.get_tracer(True)

        # uniform environment
        paths = tracer.GetSignalPaths(self.start, self.target_even, self.env_flat)
        self.assertEqual(len(paths), 2)

        paths = tracer.GetSignalPaths(self.start, self.target_general, self.env_flat)
        self.assertEqual(len(paths), 2)

        paths = tracer.GetSignalPaths(
            self.start, self.target_across_plane, self.env_flat
        )  # only one solution due to being on far side of plane
        self.assertEqual(len(paths), 1)

        # exponential environment
        paths = tracer.GetSignalPaths(self.start, self.target_even, self.env_expo)
        self.assertEqual(len(paths), 2)

        paths = tracer.GetSignalPaths(self.start, self.target_general, self.env_expo)
        self.assertEqual(len(paths), 2)

        paths = tracer.GetSignalPaths(
            self.start, self.target_across_plane, self.env_expo
        )
        self.assertEqual(len(paths), 2)


if __name__ == "__main__":
    unittest.main()
