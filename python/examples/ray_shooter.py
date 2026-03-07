from c8_tracer.c8_tracer_ext import logging

from c8_tracer.c8_tracer_ext.tracer import RayTracer2D
from c8_tracer.c8_tracer_ext import Vec3, Plane
from c8_tracer.c8_tracer_ext.environment import CartesianSingleExponentialEnvironment
from c8_tracer.c8_tracer_ext.environment import EnvironmentBase

"""
Example script showing the how a custom environment can be set up and how the ray
tracer can be constructed to launch rays into the ice from that air. By default, this
uses the ShootToMinimumZ functions will just just propagate a ray until a minimum z
has been reached regardless of how close to the target it is
"""


class TestNRMCInterface():
    def __init__(self):
        from c8_tracer.c8_tracer_ext import environment
        logging.logger_tracer.set_level(logging.LogLevel.TRACE)


        n_deep = 1.78
        delta_n = 0.423
        length_scale = 77.0
        axis = Vec3(0, 0, 1)
        ref_point = Vec3(0, 0, 0)

        self.env_flat = environment.IsotropicEnvironment(n_deep + delta_n)
        self.env_expo = environment.CartesianSingleExponentialEnvironment(
            n_deep,
            delta_n,
            length_scale,
            axis,
            ref_point,
        )

        self.axis_of_symmetry = axis
        self.perp_dir = (Vec3(2.3, 1.1, 1.3).cross(axis)).normalized()
        self.min_step = 0.00001
        self.max_step = 10.0
        self.tolerance = 1e-9
        self.nRays = 13

        self.start = Vec3(0, 0, 0)
        self.target_even = self.start + self.perp_dir * 50.0
        self.target_general = self.target_even + self.axis_of_symmetry * 40.0
        self.plane_height = self.start + self.axis_of_symmetry * 100.0
        self.target_across_plane = self.plane_height + self.perp_dir * 50.0

        tracer = self.get_tracer(True)
        paths = tracer.GetSignalPaths(self.start, self.target_even, self.env_flat)
        print(len(paths))
        print(paths)

    def get_tracer(self, with_mirror=False):
        tracer = RayTracer2D(
            self.axis_of_symmetry,
            self.min_step,
            self.max_step,
            self.tolerance,
            self.nRays,
        )

        if with_mirror:
            tracer.AddReflectionLayer(
                Plane(self.plane_height, self.axis_of_symmetry)
            )

        return tracer

TestNRMCInterface()

exit()


class AirAndIce(EnvironmentBase):
    def __init__(self):
        super().__init__()
        n_deep = 1.78
        delta_n = 0.423
        length_scale = 77.0
        axis = Vec3(0, 0, 1)
        ref_point = Vec3(0, 0, 0)
        self.ice = CartesianSingleExponentialEnvironment(
            n_deep,
            delta_n,
            length_scale,
            axis,
            ref_point,
        )

    def get_n(self, pos: Vec3) -> float:
        if pos.z > 0.0:  # in air
            return 1.0001

        return self.ice.get_n(pos)  # in ice

    def get_grad_n(self, pos: Vec3) -> Vec3:
        if pos.z > 0.0:  # in air
            return Vec3(0, 0, 0)

        return self.ice.get_grad_n(pos)  # in ice


env = AirAndIce()

# make the ray tracer
axis_of_symmetry = Vec3(0, 0, 1)  # z_hat
min_step = 0.0001  # smallest step size that will be taken
max_step = 1.0  # largest step size that will be taken
tolerance = 1e-8  # relative error tolerance that defines the adaptive step size
ray_tracer = RayTracer2D(axis_of_symmetry, min_step, max_step)

# make a reflection plane at the origin facing downwards
center = Vec3(0, 0, 0)
normal = Vec3(0, 0, -1)
plane = Plane(center, normal)
ray_tracer.AddReflectionLayer(plane)

start = Vec3(0, 0, 50)  # Start in the air
target = Vec3(300, 0, -2000)
startDir = Vec3(1, 0, -1).normalized()  # Head down into the ice

testPos = Vec3(0, 0, 0)
testDir = Vec3(0, 0, 0)

logging.logger_tracer.set_level(logging.LogLevel.TRACE)

ray_tracer.ShootOneRayToMinimumZ(start, startDir, testPos, testDir, target, env)
print("Final pos:", testPos)
print(
    "Lateral distance to target", ray_tracer.Get2DRadialDistance(start, target, testPos)
)
print()

# try a better guess
startDir = Vec3(0.2509, 0, -1).normalized()
ray_tracer.ShootOneRayToMinimumZ(start, startDir, testPos, testDir, target, env)
print("Final pos:", testPos)
print(
    "Lateral distance to target", ray_tracer.Get2DRadialDistance(start, target, testPos)
)
