import unittest

from c8_tracer_py import Vec3
from c8_tracer_py import environment


class TestEnvironment(unittest.TestCase):
    def test_base(self):
        n0 = 1.5
        gradient = Vec3(0.0, 0.0, 0.0)

        class TestEnvUniform(environment.EnvironmentBase):
            def __init__(self):
                super().__init__()

            def get_n(self, position: Vec3) -> float:
                return n0

            def get_grad_n(self, position: Vec3) -> Vec3:
                return gradient

        env = TestEnvUniform()

        self.assertEqual(env.get_n(Vec3(0, 0, 0)), n0)
        self.assertEqual(env.get_grad_n(Vec3(0, 0, 0)), gradient)

    def test_base_no_override_virtual(self):
        class NoOverride(environment.EnvironmentBase):
            def __init__(self):
                super().__init__()

        with self.assertRaises(RuntimeError):
            env = NoOverride()
            env.get_n(Vec3(0, 0, 0))

    def test_isotropic(self):
        n0 = 1.5
        env = environment.IsotropicEnvironment(n0)
        self.assertEqual(env.get_n(Vec3(0, 0, 0)), n0)
        self.assertEqual(env.get_n(Vec3(1000, 1000, 1000)), n0)

    def test_linear_radial(self):
        n0 = 1.5
        center = Vec3(0, 0, 1)
        dn_dr = 10.0
        env = environment.LinearRadialEnvironment(center, n0, dn_dr)

        # self.assertEqual(env.get_n(center), n0)
        self.assertNotEqual(env.get_n(Vec3(0.0, 0.0, 2.0)), n0)
        self.assertNotEqual(env.get_n(Vec3(1000, 1000, 1000)), n0)
        self.assertAlmostEqual(env.get_n(Vec3(0.0, 0.0, 2.0)), n0 + dn_dr)


if __name__ == "__main__":
    unittest.main()
