import logging

import numpy as np

from c8_tracer.nrmc_interface import (
    CreateNRMCWrappedRayTracer,
    CreateNRMCInterpolationTable,
)

from NuRadioReco.utilities import units
from NuRadioMC.SignalProp.propagation_base_class import ray_tracing_base
from NuRadioMC.SignalProp.propagation import solution_types_revert


class C8RayTracerIndividual(ray_tracing_base):
    """
    C8 ray tracer in the style of the NRMC ray tracing base class
    This ray tracer finds the solutions each time they are queried
    and returns the path
    """

    def __init__(
        self,
        medium,
        attenuation_model=None,
        log_level=logging.NOTSET,
        n_frequencies_integration=None,
        n_reflections=None,
        config=None,
        detector=None,
        min_step: float = 0.0001,
        max_step: float = 1.0,
        tolerance: float = 1e-8,
        n_rays: int = 21,
    ):
        """
        C8 ray tracer in the style of the NRMC ray tracing base class
        This ray tracer finds the solutions each time they are queried
        and returns the path

        :param medium: class that inherits the `IceModel` class
        :param attenuation_model: unused by this class
        :param log_level: controls the granualrity of the print messages
        :param n_frequencies_integration: unused by this class
        :param config: unused by this class
        :param detector: unused by this class
        :param min_step: smallest step size that will be taken during propagation (meters)
        :param max_step: largest step size that will be taken during propagation (meters)
        :param tolerance: relative error tolerance that defines the adaptive step size

        """
        if n_reflections != 1:
            raise RuntimeError(
                "This tracer only works for a single reflection at the surface"
            )

        super().__init__(
            medium=medium,
            attenuation_model=attenuation_model,
            log_level=log_level,
            n_frequencies_integration=n_frequencies_integration,
            n_reflections=n_reflections,
            config=config,
            detector=detector,
        )

        # function that returns a signal path for each of ONLY 2 solutions
        self._trace_to_point, self._tracer = CreateNRMCWrappedRayTracer(
            medium,
            min_step,
            max_step,
            tolerance,
            nRays=n_rays,
        )

        self.set_config(config=config)

    def reset_solutions(self):
        self._X1 = None
        self._X2 = None
        self._results = None
        self._paths = None

    def find_solutions(self):
        if self._X1 is not None and self._X2 is not None:
            self._results = self._trace_to_point(self._X1, self._X2)
            self._paths = [None] * len(self._results)
            return

        self.__logger.fatal(
            "Asking to find solutions without first settings the end points"
        )
        raise RuntimeError()

    def has_solution(self):
        return self._results is not None and len(self._results) > 0

    def get_solution_type(self, iS):
        if iS == 0:
            return solution_types_revert["direct"]

        return solution_types_revert["refracted"]

    def get_path(self, iS, n_points=None):
        # ensure solutions exist
        if self._results is None:
            self.find_solutions()

        assert self._paths is not None
        assert self._results is not None

        # results are stored as Vec3 unless asked for, convert
        if self._paths[iS] is None:
            self._paths[iS] = np.array(
                [np.array([x.x, x.y, x.z]) for x in self._results[iS]]
            )

        # interpolate points if number specified
        if n_points is not None:
            x = np.linspace(0.0, 1.0, n_points)
            xp = np.linspace(0.0, 1.0, len(self._paths[iS]))
            return np.transpose(
                [
                    np.interp(x, xp, self._paths[iS, :, 0]), # type: ignore
                    np.interp(x, xp, self._paths[iS, :, 1]),
                    np.interp(x, xp, self._paths[iS, :, 2]),
                ]
            )

        # return the native path otherwise
        return self._paths[iS]

    def get_launch_vector(self, iS):
        assert self._results is not None
        return self._results[iS].emit

    def get_receive_vector(self, iS):
        assert self._results is not None
        return -self._results[iS].receive

    def get_reflection_angle(self, iS):
        if iS != 1:
            return None

        # find the point closest to z and ensure that it is plausibly
        # a reflection, note that the tracer will ensure that a path point
        # is with 5e-6 m of a reflection plane, so it won't be subtle
        path = self.get_path(iS)
        z_ref = self._medium.z_air_boundary

        dz = np.abs(path[:, 2] - z_ref)

        if not np.any(dz < 1e-2):
            return None

        imin = np.argmin(dz)
        dr = path[imin] - path[imin - 1]
        dr /= np.linalg.norm(dr)

        return abs(np.arccos(dr[2]))

    def get_path_length(self, iS, analytic=False):
        assert self._results is not None
        return self._results[iS].R_distance * units.m

    def get_travel_time(self, iS, analytic=False):
        assert self._results is not None
        return self._results[iS].propagation_time * units.s


class C8RayTracerTable(ray_tracing_base):
    def __init__(
        self,
        medium,
        antenna_positions: list[np.ndarray],
        maxR: float,
        minZ: float,
        maxZ: float,
        bin_length: float,
        min_step: float = 0.0001,
        max_step: float = 1.0,
        tolerance: float = 1e-8,
        n_rays: int = 21,
        attenuation_model=None,
        log_level=logging.NOTSET,
        n_frequencies_integration=None,
        n_reflections=None,
        config=None,
        detector=None,
    ):
        """
        C8 ray tracer in the style of the NRMC ray tracing base class
        This ray tracer pre-calculates solutions and stores them in tables.
        One table is made for each unique z-position.
        It is recommended that `maxZ` is not exactly at the surface

        :param medium: class that inherits the `IceModel` class

        :param antenna_positions: (N, 3) list of positions of the antennas, one table
        will be built for each unique z-coordinate
        :param maxR: maximum radius to create the table for (meters)
        :param minZ: minimum height to make the table
        :param maxZ: maximum height to make the table
        :param bin_length: the span between interpolation nodes
        :param min_step: smallest step size that will be taken during propagation (meters)
        :param max_step: largest step size that will be taken during propagation (meters)
        :param tolerance: relative error tolerance that defines the adaptive step size

        :param attenuation_model: unused by this class
        :param log_level: controls the granualrity of the print messages
        :param n_frequencies_integration: unused by this class
        :param config: unused by this class
        :param detector: unused by this class

        """
        if n_reflections != 1:
            raise RuntimeError(
                "This tracer only works for a single reflection at the surface"
            )

        super().__init__(
            medium=medium,
            attenuation_model=attenuation_model,
            log_level=log_level,
            n_frequencies_integration=n_frequencies_integration,
            n_reflections=n_reflections,
            config=config,
            detector=detector,
        )

        self._trace_to_point, self._tracer = CreateNRMCInterpolationTable(
            medium,
            antenna_positions,
            maxR,
            minZ,
            maxZ,
            bin_length,
            min_step,
            max_step,
            tolerance,
            n_rays,
        )

        self.set_config(config=config)

    def reset_solutions(self):
        self._X1 = None
        self._X2 = None
        self._results = None
        self._paths = None

    def find_solutions(self):
        if self._X1 is not None and self._X2 is not None:
            self._results = self._trace_to_point(self._X1, self._X2)
            self._paths = [None] * len(self._results)
            return

        self.__logger.fatal(
            "Asking to find solutions without first settings the end points"
        )
        raise RuntimeError()

    def has_solution(self):
        return self._results is not None and len(self._results) > 0

    def get_solution_type(self, iS):
        if iS == 0:
            return solution_types_revert["direct"]

        return solution_types_revert["refracted"]

    def get_path(self, iS, n_points=None):
        raise NotImplementedError("Paths are not stored in the interpolation table")

    def get_launch_vector(self, iS):
        assert self._results is not None
        return self._results[iS].emit

    def get_receive_vector(self, iS):
        assert self._results is not None
        return self._results[iS].receive

    def get_reflection_angle(self, iS):
        raise NotImplementedError("Paths are not stored in the interpolation table")

    def get_path_length(self, iS, analytic=False):
        assert self._results is not None
        return self._results[iS].R_distance * units.m

    def get_travel_time(self, iS, analytic=False):
        assert self._results is not None
        return self._results[iS].propagation_time * units.s
