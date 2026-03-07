import numpy as np

from c8_tracer.nrmc_interface import ConvertSignalPath, CreateNRMCInterpolationTable
from c8_tracer.c8_tracer_ext import logging
from c8_tracer.c8_tracer_ext.tracer import SolutionMethod

"""
This sketches out how this class could be used in the main NRMC loop for simulating
events by first building an interpolation table
"""


# Can be quite loud without this
logging.logger.set_level(logging.LogLevel.INFO)
logging.logger_tracer.set_level(logging.LogLevel.ERROR)


# this is a mock up of the NRMC `IceModel` base class
# so that it can be inherited below
class DummyIceModel:
    def __init__(self):
        self.z_air_boundary = 0.0

    def get_index_of_refraction(self, _):
        return 1.33

    def get_gradient_of_index_of_refraction(self, _):
        return np.zeros(3)


# would normally come from the NRMC simulation settings
ice_model = DummyIceModel()
antennas = [
    np.array([50.0, 0.0, -100.0]),
    np.array([50.0, 0.0, -150.0]),
    np.array([50.0, 0.0, -200.0]),
]
shower_pos = np.array([-50.0, 0.0, -150.0])

# This will create the tables and return a function that smoothly converts to/from NRMC
get_path_to_antenna, tracer = CreateNRMCInterpolationTable(
    ice_model,
    antennas,
    maxR=300.0,
    minZ=-300.0,  # for SP should be more like -3000
    maxZ=-0.1,
    bin_length=10.0,
    min_step=0.0001,
    max_step=1.0,
    tolerance=1e-8,
    nRays=10,
    method=SolutionMethod.Brent,
)

# would be the main loop in the NRMC code
for ant_pos in antennas:
    signal_paths = get_path_to_antenna(shower_pos, ant_pos)

    print("Paths from", shower_pos, "to", ant_pos)

    for signal_path in signal_paths:  # ordered [direct, refracted]
        dt, length, emit, receive = ConvertSignalPath(signal_path)
        print(f"\t dt={dt * 1e9:.1f}ns, L={length:.1f}m, emit={emit}, rec={receive}")

        # ...do whatever NRMC normally does
