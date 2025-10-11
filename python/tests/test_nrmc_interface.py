import unittest

import numpy as np

from c8_tracer import nrmc_interface
from c8_tracer.c8_tracer_ext import Vec3


class FakeNRMCModel:
    def __init__(self):
        self.z_air_boundary = 0.0

    def get_index_of_refraction(self, _):
        return 1.33

    def get_gradient_of_index_of_refraction(self, _):
        return np.zeros(3)


class TestNRMCInterface(unittest.TestCase):
    def setUp(self):
        self.nrmc_model = FakeNRMCModel()

    def test_WrappedEnvironment(self):
        env = nrmc_interface.WrappedEnvironment(self.nrmc_model)

        true_n = self.nrmc_model.get_index_of_refraction([0, 0, 0])
        true_grad = Vec3(self.nrmc_model.get_gradient_of_index_of_refraction([0, 0, 0]))

        self.assertAlmostEqual(env.get_n(Vec3(1, 2, 3)), true_n)
        self.assertAlmostEqual(env.get_n(Vec3(0, 0, 0)), true_n)

        self.assertAlmostEqual((env.get_grad_n(Vec3(0, 0, 0)) - true_grad).norm(), 0.0)
        self.assertAlmostEqual((env.get_grad_n(Vec3(1, 2, 3)) - true_grad).norm(), 0.0)

    def test_CreateNRMCRayTracer(self):
        tracer = nrmc_interface.CreateNRMCRayTracer(
            self.nrmc_model,
            0.0001,
            1.234,
            1e-7,
            15,
        )

        tracer.PrintProfiling()

        # perp to axis
        proj = tracer.Get2DProjection(Vec3(0, 0, 0), Vec3(1, 0, 0))
        self.assertAlmostEqual(proj.norm(), 1.0)

        # perp and parallel to axis
        proj = tracer.Get2DProjection(Vec3(0, 0, 0), Vec3(1, 0, 1))
        self.assertAlmostEqual(proj.norm(), 1.0)

        # parallel to axis
        proj = tracer.Get2DProjection(Vec3(0, 0, 0), Vec3(0, 0, 1))
        self.assertAlmostEqual(proj.norm(), 0.0)

    def test_CreateNRMCWrappedRayTracer(self):
        trace_to_point, ray_tracer = nrmc_interface.CreateNRMCWrappedRayTracer(
            self.nrmc_model,
            0.0001,
            1.0,
            1e-8,
            13,
        )

        path_list = trace_to_point(
            np.array([1.0, 0.0, -100.0]), np.array([50.0, 0.0, -100.0])
        )
        self.assertEqual(len(path_list), 2)

        proj = ray_tracer.Get2DProjection(Vec3(0, 0, 0), Vec3(1, 0, 0))
        self.assertAlmostEqual(proj.norm(), 1.0)

    def test_CreateNRMCInterpolationTable(self):
        antenna_pos = np.array([0, 0, -100.0])
        get_path_to_antenna, ray_tracer = nrmc_interface.CreateNRMCInterpolationTable(
            self.nrmc_model,
            [antenna_pos],  # ant_position_list
            50.0,  # maxR
            -150,  # minZ
            -50,  # maxZ
            25,  # bin_length
            0.0001,  # min_step
            1.0,  # max_step
            1e-8,  # tolerance
            13,  # nRays
        )

        proj = ray_tracer.Get2DProjection(Vec3(0, 0, 0), Vec3(1, 0, 0))
        self.assertAlmostEqual(proj.norm(), 1.0)

        signal_paths = get_path_to_antenna(
            np.array([20.0, 20.0, -70.0]),
            antenna_pos,
        )
        self.assertEqual(len(signal_paths), 2)

        with self.assertRaises(RuntimeError):
            wrong_pos = np.array(antenna_pos)
            wrong_pos[2] += 1.0

            get_path_to_antenna(
                np.array([20.0, 20.0, -70.0]),
                wrong_pos,
            )
