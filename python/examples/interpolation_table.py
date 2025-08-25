import numpy as np
import matplotlib.pyplot as plt

from c8_tracer_py import logging
from c8_tracer_py import Plane
from c8_tracer_py.tracer import RayTracer2D
from c8_tracer_py.tables import InterpolationTableGenerator2D
from c8_tracer_py import Vec3
from c8_tracer_py.environment import CartesianSingleExponentialEnvironment

logging.logger.set_level(logging.LogLevel.ERROR)

# set up the exponential profile
# n(z) = n_deep - delta_n * exp((z - z0) / length_scale)
n_deep = 1.78
delta_n = 0.423
length_scale = 77.0
axis = Vec3(0, 0, 1)  # points against the gradient (towards zenith for glacial ice)
ref_point = Vec3(0, 0, 0)
env = CartesianSingleExponentialEnvironment(
    n_deep,
    delta_n,
    length_scale,
    axis,
    ref_point,
)

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

table_gen = InterpolationTableGenerator2D(ray_tracer, env)

# define the location from which all the rays will be traced
# must be below the reflection plane
antenna_pos = Vec3(0, 0, -100)

min_r = 1
max_r = 300
n_bins_r = 13
min_z = -200  # negative values cover locations below the antenna
max_z = -antenna_pos.z - 1  # 1m below the reflection plane
n_bins_z = 9
tables = table_gen.GenerateTables(
    min_r, max_r, n_bins_r, min_z, max_z, n_bins_z, antenna_pos
)

for table in tables:
    table.Print()

target_pos = antenna_pos + Vec3(max_r / 2.0, 0.0, antenna_pos.z + max_z / 2.0)

print("\n\nTesting propagation for", target_pos)

for itable, table in enumerate(tables):
    if not table.ContainsPoint(target_pos):
        print("TABLE DOES NOT CONTAIN TARGET", target_pos)
        continue

    print(f"------- Ray solution {itable + 1} -------")
    path = table.GetSignalPath(target_pos, 1.0, 1.0)
    print("    Launch:   ", path.emit)
    print("    Receive:  ", path.receive)
    print("    Duration: ", path.propagation_time)
    print("    Length:   ", path.R_distance)
    print("    FresnelS: ", path.fresnelS)
    print("    FresnelP: ", path.fresnelP)
