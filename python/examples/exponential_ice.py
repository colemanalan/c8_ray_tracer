import numpy as np
import matplotlib.pyplot as plt

from c8_tracer_py import logging
from c8_tracer_py import Plane
from c8_tracer_py import RayTracer2D
from c8_tracer_py import Vec3
from c8_tracer_py.environment import CartesianSingleExponentialEnvironment
from c8_tracer_py.path import SignalPath

logging.logger.set_level(logging.LogLevel.INFO)

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

# propagate from start to finish
start = Vec3(0, 0, -50)
end = Vec3(300, 0, -200.0)
paths: list[SignalPath] = ray_tracer.PropagateToPoint(start, end, env)

for path in paths:
    dist = (path.getEnd() - end).norm()
    logging.logger.info(f"Distance to target {dist}m")

print()
ray_tracer.PrintProfiling()

# make a plot of the final path
fig = plt.figure(1)
ax = fig.add_subplot(1, 2, 1)
for ipath, path in enumerate(paths):
    steps = np.array([[loc.x, loc.z] for loc in path])
    ax.plot(steps.T[0], steps.T[1], label=f"sol {ipath + 1}")

ax.scatter([start.x], [start.z], color="k", marker="o", label="start")
ax.scatter([end.x], [end.z], color="r", marker="s", label="end")
ax.legend()
ax.set_aspect("equal")

# make a plot of the density profile
heights = np.linspace(*ax.get_ylim(), 100)
n_indexes = [env.get_n(Vec3(0, 0, z)) for z in heights]
ax = fig.add_subplot(1, 2, 2)
ax.plot(n_indexes, heights)
ax.set_xlabel("Refractive index")

fig.savefig("ExponentialIcePath.png")
