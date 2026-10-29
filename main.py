""" Stripmap SAR simulation with an FMCW radar """

# %% Preamble
import matplotlib.pyplot as plt

from params import Params
from geometry import Geometry
from signal_generator import SignalGenerator
from rma import RangeMigrationAlgorithm


# %% Radar parameters
params = Params(fc=10e9, fs=10e6, R_max=3e3, bandwidth=150e6)


# %% Scene: point targets on the ground
# x = along-track position [m], y = ground range [m], z = height [m], amplitude = reflectivity [linear]
targets = [
    {'x':   0.0, 'y': 800.0, 'z': 0.0, 'amplitude': 1.0},  # Scene center
    {'x': -15.0, 'y': 800.0, 'z': 0.0, 'amplitude': 1.0},  # Same range, different azimuth
    {'x':  15.0, 'y': 800.0, 'z': 0.0, 'amplitude': 1.0},
    {'x':   0.0, 'y': 780.0, 'z': 0.0, 'amplitude': 1.0},  # Same azimuth, different range
    {'x':   0.0, 'y': 820.0, 'z': 0.0, 'amplitude': 1.0},
]


# %% Geometry: platform, antenna and scene
# beamwidth: 'wide' (12 deg), 'medium' (6 deg), 'narrow' (3 deg), or any number in degrees
geometry = Geometry(params, targets, beamwidth='medium', platform_speed=100.0, altitude=300.0)
geometry.plot_geometry()


# %% Raw data
# stop_and_go: True = platform frozen during each chirp; False = platform moving (realistic). Used by the RMA too
stop_and_go = True
sgen = SignalGenerator(params, geometry, stop_and_go=stop_and_go, snr_db=None)
raw_data = sgen.generate()
# sgen.plot_range_compressed(raw_data)  # Quick look at the raw data (the RMA plots it too)


# %% Image formation: Range Migration Algorithm
rma = RangeMigrationAlgorithm(params, geometry, stop_and_go=stop_and_go)
image = rma.focus(raw_data, plot=True)

plt.show()
print('here')