import unittest
from unittest.mock import MagicMock
from c8_tracer_py import Vec3
from c8_tracer_py.integrators import CashKarpIntegrator
from c8_tracer_py.environment import EnvironmentBase


class TestCashKarpIntegrator(unittest.TestCase):
    def setUp(self):
        class TestEnvUniform(EnvironmentBase):
            def __init__(self):
                super().__init__()

            def get_n(self, position: Vec3) -> float:
                return 1.5

            def get_grad_n(self, position: Vec3) -> Vec3:
                return Vec3(0.0, 0.0, 0.0)

        class TestEnvNonUniform(EnvironmentBase):
            def __init__(self):
                super().__init__()

            def get_n(self, position: Vec3) -> float:
                return 1.5

            def get_grad_n(self, position: Vec3) -> Vec3:
                return Vec3(100.0, 100.0, 100.0)

        # Mock environment with constant refractive index and zero gradient
        self.env_uniform = TestEnvUniform()
        self.env_non = TestEnvNonUniform()

        self.start_pos = Vec3(0.0, 0.0, 0.0)
        self.start_dir = Vec3(1.0, 0.0, 0.0)
        minStep = 0.01
        self.maxStep = 1.0
        tolerance = 0.001
        self.integrator = CashKarpIntegrator(minStep, self.maxStep, tolerance)

    def test_step_basic(self):
        end_pos = Vec3(0.0, 0.0, 0.0)
        end_dir = Vec3(0.0, 0.0, 0.0)
        dir_error = Vec3(0.0, 0.0, 0.0)
        h = 0.1

        (end_pos, end_dir, dir_error) = self.integrator.Step(
            self.start_pos,
            self.start_dir,
            h,
            self.env_uniform,
        )

        # Check that end_pos and end_dir have changed
        self.assertNotEqual(end_pos, self.start_pos)
        self.assertEqual(end_dir, self.start_dir)
        self.assertTrue(dir_error.norm() >= 0.0)

        (end_pos, end_dir, dir_error) = self.integrator.Step(
            self.start_pos,
            self.start_dir,
            h,
            self.env_non,
        )

        # Check that end_pos and end_dir have changed
        self.assertNotEqual(end_pos, self.start_pos)
        self.assertNotEqual(end_dir, self.start_dir)
        self.assertTrue(dir_error.norm() >= 0.0)

    def test_adaptive_step_converges(self):
        h0 = 0.5

        end_pos, end_dir, h0 = self.integrator.AdaptiveStep(
            self.start_pos, self.start_dir, self.env_uniform, True
        )

        # Ensure output is updated
        self.assertNotEqual(end_pos, self.start_pos)
        self.assertEqual(end_dir, self.start_dir)
        self.assertTrue(0.01 <= h0 <= 1.0)

    def test_adaptive_step_min_step(self):
        h0 = 0.5

        _, _, h0 = self.integrator.AdaptiveStep(
            self.start_pos, self.start_dir, self.env_non, True
        )

        # Should clamp to minStep due to high error
        self.assertAlmostEqual(h0, self.maxStep)

    def test_constructor_parameters(self):
        integrator = CashKarpIntegrator(0.05, 2.0, 0.01)
        self.assertIsInstance(integrator, CashKarpIntegrator)


if __name__ == "__main__":
    unittest.main()
