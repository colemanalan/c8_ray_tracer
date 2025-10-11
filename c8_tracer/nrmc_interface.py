import numpy as np

from c8_tracer.c8_tracer_ext import Vec3
from c8_tracer.c8_tracer_ext.environment import EnvironmentBase
from c8_tracer.c8_tracer_ext.path import SignalPath
from c8_tracer.c8_tracer_ext.tracer import RayTracer2D


class WrappedEnvironment(EnvironmentBase):
    """
    Creates a wrapped Environment class based a NRMC ice model

    :param nrmc_model: any NRMC ice model that inherits from `IceModel`
    """

    def __init__(self, nrmc_model):
        super().__init__()
        self.nrmc_model = nrmc_model

    def get_n(self, pos: Vec3) -> float:
        return self.nrmc_model.get_index_of_refraction([pos.x, pos.y, pos.z])

    def get_grad_n(self, pos: Vec3) -> Vec3:
        grad_nrmc = self.nrmc_model.get_gradient_of_index_of_refraction(
            [pos.x, pos.y, pos.z]
        )
        return Vec3(grad_nrmc[0], grad_nrmc[1], grad_nrmc[2])


def CreateNRMCRayTracer(
    nrmc_model,
    min_step: float = 0.0001,
    max_step: float = 1.0,
    tolerance: float = 1e-8,
    nRays: float = 13,
) -> RayTracer2D:
    """
    Creates a ray tracer that is initialized for use with NRMC simulations

    :param nrmc_model: any NRMC ice model that inherits from `IceModel`
    :param min_step: smallest step size that will be taken during propagation (meters)
    :param max_step: largest step size that will be taken during propagation (meters)
    :param tolerance: relative error tolerance that defines the adaptive step size
    :param nRays: defines how many rays are cast in the initial scan each time a solution
    is searched for (more will run slowed but will fail less often)

    :return: initialized instance of a C8 raytracer
    """
    axis_of_symmetry = Vec3(0, 0, 1)  # z_hat
    ray_tracer = RayTracer2D(axis_of_symmetry, min_step, max_step, tolerance, nRays)

    if nrmc_model.z_air_boundary is not None:
        from c8_tracer.c8_tracer_ext import Plane

        # make a reflection plane at the origin facing downwards
        ice_surface = Vec3(0, 0, nrmc_model.z_air_boundary)
        normal = Vec3(0, 0, -1)
        plane = Plane(ice_surface, normal)
        ray_tracer.AddReflectionLayer(plane)

    return ray_tracer


def CreateNRMCWrappedRayTracer(
    nrmc_model,
    min_step: float = 0.0001,
    max_step: float = 1.0,
    tolerance: float = 1e-8,
    nRays: float = 13,
):
    """
    Creates a ray tracer and wraps it in a function that takes start and end positions
    in the usual NRMC coord sys that can be used for finding solutions during simulations

    :param nrmc_model: any NRMC ice model that inherits from `IceModel`
    :param min_step: smallest step size that will be taken during propagation (meters)
    :param max_step: largest step size that will be taken during propagation (meters)
    :param tolerance: relative error tolerance that defines the adaptive step size
    :param nRays: defines how many rays are cast in the initial scan each time a solution
    is searched for (more will run slowed but will fail less often)

    :return: function that wraps an initialized ray tracer. The function takes two
    np.ndarrays as pos1 and pos2 and returns the direct and refracted solutions
    """

    env = WrappedEnvironment(nrmc_model)

    ray_tracer = CreateNRMCRayTracer(nrmc_model, min_step, max_step, tolerance, nRays)

    def trace_to_point(
        pos1: np.ndarray, pos2: np.ndarray
    ) -> list[SignalPath]:
        
        # Must be different altitudes for proper z-check
        if pos1[2] == pos2[2]:
            pos2[2] -= 1e-9

        start = Vec3(pos1[0], pos1[1], pos1[2])
        end = Vec3(pos2[0], pos2[1], pos2[2])
        paths: list[SignalPath] = ray_tracer.GetSignalPaths(start, end, env)

        return paths

    return trace_to_point, ray_tracer


def CreateNRMCInterpolationTable(
    nrmc_model,
    ant_position_list: list[np.ndarray],
    maxR: float,
    minZ: float,
    maxZ: float,
    bin_length: float,
    min_step: float = 0.0001,
    max_step: float = 1.0,
    tolerance: float = 1e-8,
    nRays: float = 13,
):
    """
    Creates a ray tracer and wraps it in a function that takes start and end positions
    in the usual NRMC coord sys that can be used for finding solutions during simulations

    :param nrmc_model: any NRMC ice model that inherits from `IceModel`

    :param ant_position_list: list of antenna positions, one ray tracing table will be
    created for each one
    :param maxR: maximum radius to be tabulated w.r.t. each antenna
    :param minZ: the minimum z to be tabulated
    :param maxZ: the maximum z to be tabulated
    :param bin_length: the side length of each bin (will be square in r/z)

    :param min_step: smallest step size that will be taken during propagation (meters)
    :param max_step: largest step size that will be taken during propagation (meters)
    :param tolerance: relative error tolerance that defines the adaptive step size
    :param nRays: defines how many rays are cast in the initial scan each time a solution
    is searched for (more will run slowed but will fail less often)

    :return get_path_to_antenna: function that can be queried, `get_path_to_antenna(start_pos, ant_pos)`
    and will return two signals paths. Signals paths that are `None` denote no solution
    :return ray_tracer: initialized ray tracer
    """

    import numpy as np
    from c8_tracer.c8_tracer_ext.tables import InterpolationTableGenerator2D

    env = WrappedEnvironment(nrmc_model)

    ray_tracer = CreateNRMCRayTracer(nrmc_model, min_step, max_step, tolerance, nRays)
    table_generator = InterpolationTableGenerator2D(ray_tracer, env)

    all_tables = {}

    for pos in ant_position_list:
        if _GetDictKey(pos) in all_tables.keys():
            msg = (
                f"Already have key {_GetDictKey(pos)} in dict. Cannot make duplicates at"
                + " this time. See: https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/merge_requests/671#note_34831"
            )
            raise RuntimeError(msg)

        # tables are defined relative to the antenna
        z_up = maxZ - pos[2]
        z_down = minZ - pos[2]
        z_bins = int(abs(z_up - z_down) / bin_length + 1)

        # cover the whole cylindrical volume
        r_min = 0.5  # this  should NEVER be exactly zero
        r_bins = int(abs(maxR - r_min) / bin_length + 1)

        vec_pos = Vec3(pos[0], pos[1], pos[2])
        tables = table_generator.GenerateTables(
            float(r_min),
            float(maxR),
            int(r_bins),
            float(z_down),
            float(z_up),
            int(z_bins),
            vec_pos,
        )

        all_tables[_GetDictKey(pos)] = tables

    # construct a function that can be queried in NRMC
    def get_path_to_antenna(
        start_pos: np.ndarray, antenna_pos: np.ndarray
    ) -> list[SignalPath]:
        # check that a table has been made for this antenna
        if _GetDictKey(antenna_pos) not in all_tables.keys():
            raise RuntimeError(f"No dict has been made for antenna at {antenna_pos}")

        tables = all_tables[_GetDictKey(antenna_pos)]

        r_i = Vec3(start_pos)
        r_f = Vec3(antenna_pos)

        signal_paths = []
        for table in tables:
            if not table.ContainsPoint(r_i):
                signal_paths.append(None)  # None for no solution paths

            signal_path: SignalPath = table.GetSignalPath(
                r_i, env.get_n(r_i), env.get_n(r_f)
            )
            signal_paths.append(signal_path)

        return signal_paths

    return get_path_to_antenna, ray_tracer


def ConvertSignalPath(signal_path: SignalPath):
    """
    Converter from C8 data class to bare numbers

    :param signal_path: path to be converted

    :return dt: time of flight (seconds)
    :return length: length of ray (meter)
    :return emit: emission unit vector
    :return receive: receive unit vector
    """

    if signal_path is None:
        return None, None, None, None

    dt = signal_path.propagation_time
    emit_ = signal_path.emit
    rec_ = signal_path.receive
    emit = np.array([emit_.x, emit_.y, emit_.z])
    receive = np.array([rec_.x, rec_.y, rec_.z])
    length = signal_path.R_distance

    return dt, length, emit, receive


def _GetDictKey(pos: np.ndarray):
    """
    Helper function to go from
    """
    return pos[2]
