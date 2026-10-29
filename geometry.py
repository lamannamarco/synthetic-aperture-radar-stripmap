""" Class containing the stripmap SAR geometry: platform trajectory, antenna and scene

Coordinate system (right-handed, in meters):
    x : along-track (azimuth), the direction in which the platform flies
    y : cross-track ground range, the side the antenna looks at
    z : height above the ground

The platform flies in a straight line along x, at y = 0 and z = altitude, with
constant speed. The antenna is side-looking (it points along +y) and, being a
stripmap SAR, it never steers: the beam stays perpendicular to the flight path
(zero squint) and sweeps over the ground as the platform moves.
"""

import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
from scipy.constants import speed_of_light as c

class Geometry:

    # Azimuth 3 dB beamwidths [deg] of the predefined antennas
    BEAMWIDTH_PRESETS_DEG = {
        'wide': 12.0,
        'medium': 6.0,
        'narrow': 3.0,
    }


    def __init__(self, params, targets, beamwidth='medium', platform_speed=100.0, altitude=300.0):
        """ Initializing the class

        params          : Params object (radar parameters)
        targets         : list of dicts, one per point target, with keys 'x', 'y', 'z' [m] and 'amplitude' [linear]
        beamwidth       : 'wide', 'medium', 'narrow' or a number (azimuth 3 dB beamwidth in degrees)
        platform_speed  : [m/s] speed of the platform along x
        altitude        : [m] height of the platform above the ground
        """

        self.params = params
        self.targets = targets
        self.platform_speed = platform_speed
        self.altitude = altitude

        # --- Antenna ---
        # Accept either a preset name or a number in degrees
        if isinstance(beamwidth, str):
            self.beam_name = beamwidth
            beamwidth_deg = Geometry.BEAMWIDTH_PRESETS_DEG[beamwidth]
        else:
            self.beam_name = 'custom'
            beamwidth_deg = beamwidth
        self.beamwidth = np.deg2rad(beamwidth_deg)  # [rad] Azimuth 3 dB beamwidth (one-way, full width)

        # A Gaussian beam is never exactly zero, so we decide to ignore the target once it is
        # further than one full 3 dB beamwidth away from boresight. At that angle the two-way
        # gain is already -24 dB, so cutting the beam there barely changes the signal
        self.illumination_angle = self.beamwidth  # [rad] Half-angle of the illuminated cone

        # --- Slow time (one sample per chirp) ---
        # FMCW chirps are sent back to back, so the pulse repetition interval (PRI) is the chirp length
        self.PRI = params.Tc_full  # [s] Time between two consecutive chirps
        self.PRF = 1 / self.PRI  # [Hz] Pulse (chirp) repetition frequency
        self.azimuth_step = self.platform_speed * self.PRI  # [m] Distance flown between two chirps

        # --- Flight path ---
        # Every target must be seen over its whole illumination interval: the platform has to start
        # before the first target enters the beam and stop after the last target leaves it.
        # A target at closest-approach slant range R0 is inside the beam while the platform is
        # within +/- R0*tan(illumination_angle) of it, along x
        target_x = np.array([target['x'] for target in targets])
        target_R0 = np.array([self.closest_approach_range(target) for target in targets])
        self.half_illumination_length = target_R0.max() * np.tan(self.illumination_angle)  # [m] Worst case (farthest target)
        self.x_start = target_x.min() - self.half_illumination_length  # [m] Platform x at the first chirp
        self.x_end = target_x.max() + self.half_illumination_length  # [m] Platform x at the last chirp

        self.Nc = int(np.ceil((self.x_end - self.x_start) / self.azimuth_step)) + 1  # Chirps needed to fly from x_start to x_end
        self.t_slow = np.arange(self.Nc) * self.PRI  # [s] Slow time: start time of each chirp
        self.x_platform = self.x_start + self.platform_speed * self.t_slow  # [m] Platform x at the start of each chirp

        # --- Scene ---
        self.R0_min = target_R0.min()  # [m] Closest target (at closest approach)
        self.R0_max = target_R0.max()  # [m] Farthest target (at closest approach)
        self.R_max_scene = self.R0_max / np.cos(self.illumination_angle)  # [m] Largest range ever seen (farthest target, at the beam edge)

        # --- Figures of merit ---
        # Synthetic aperture length: distance flown while a target is inside the 3 dB beam
        self.L_sa = 2 * self.R0_max * np.tan(self.beamwidth / 2)  # [m]
        # Doppler of a target at angle theta from boresight: f_D = 2 * v * sin(theta) / lambda.
        # Over the illuminated cone, Doppler spans from -f_D(illumination_angle) to +f_D(illumination_angle)
        self.doppler_bandwidth = 4 * self.platform_speed * np.sin(self.illumination_angle) / params.lc  # [Hz]
        # Azimuth resolution of a stripmap SAR (it does not depend on range!)
        self.azimuth_resolution = params.lc / (4 * np.sin(self.beamwidth / 2))  # [m]

        Geometry.check_geometry(self)
        Geometry.print_geometry_to_console(self)


    def closest_approach_range(self, target):
        """ Slant range [m] between the target and the flight path, when the platform is right abeam of it """

        return np.sqrt(target['y']**2 + (target['z'] - self.altitude)**2)


    def platform_position(self, t):
        """ Platform position (x, y, z) [m] at time t [s] (t can be a scalar or an array) """

        x = self.x_start + self.platform_speed * t
        y = np.zeros_like(x)
        z = np.full_like(x, self.altitude)
        return x, y, z


    def antenna_gain(self, angle):
        """ Two-way antenna amplitude gain for a given azimuth angle [rad] from boresight

        One-way Gaussian power pattern: G(theta) = exp(-4 ln2 (theta / beamwidth)^2), which is 0.5 (-3 dB) at
        theta = beamwidth/2. The received signal amplitude goes through the antenna twice (transmit and
        receive), i.e. sqrt(G) * sqrt(G) = G. The elevation pattern is ignored: we assume the elevation beam
        is wide enough to illuminate the whole scene evenly.
        """

        gain = np.exp(-4 * np.log(2) * (angle / self.beamwidth)**2)
        gain = np.where(np.abs(angle) <= self.illumination_angle, gain, 0.0)  # Ignore the target outside the illuminated cone
        return gain


    def check_geometry(self):
        """ Warning the user about geometries that the radar cannot handle """

        if self.R_max_scene > self.params.R_max:
            print(f'\nWARNING: the scene extends to {self.R_max_scene:.1f} m, beyond the radar maximum range '
                  f'({self.params.R_max:.1f} m). Those echoes will alias in range.')

        if self.doppler_bandwidth > self.PRF:
            print(f'\nWARNING: the Doppler bandwidth ({self.doppler_bandwidth:.1f} Hz) exceeds the PRF ({self.PRF:.1f} Hz). '
                  f'The azimuth signal is undersampled and will alias: reduce the beamwidth or the platform speed.')


    def print_geometry_to_console(self):
        """ Printing to console all the geometry parameters """

        data = [
                ['Antenna', f'{self.beam_name}', '/', ''],
                ['Azimuth beamwidth', f'{np.rad2deg(self.beamwidth):.3f}', 'deg', '3 dB, one-way'],
                ['Platform speed', f'{self.platform_speed:.3f}', 'm/s', ''],
                ['Altitude', f'{self.altitude:.3f}', 'm', ''],
                ['PRF', f'{self.PRF:.3f}', 'Hz', '1/Tc (chirps back to back)'],
                ['Azimuth step', f'{self.azimuth_step*1e3:.3f}', 'mm', 'Distance flown between chirps'],
                ['Doppler bandwidth', f'{self.doppler_bandwidth:.3f}', 'Hz', 'Must be < PRF'],
                ['Closest range', f'{self.R0_min:.3f}', 'm', 'Slant range, closest approach'],
                ['Farthest range', f'{self.R0_max:.3f}', 'm', 'Slant range, closest approach'],
                ['Synthetic aperture', f'{self.L_sa:.3f}', 'm', 'Farthest target, 3 dB beam'],
                ['Flight path', f'{self.x_end - self.x_start:.3f}', 'm', ''],
                ['Nc', f'{self.Nc}', '/', 'Chirps over the flight path'],
                ['Acquisition time', f'{self.t_slow[-1]:.3f}', 's', ''],
                ['Azimuth resolution', f'{self.azimuth_resolution*1e2:.3f}', 'cm', 'Theoretical'],
        ]
        headers = ['Item', 'Value', 'Unit', 'Comments']
        table = tabulate(data, headers=headers, disable_numparse=True)
        print('\nGEOMETRY TABLE')
        print(table)


    def plot_geometry(self):
        """ Top view (x-y) and side view (y-z) of the platform, the beam and the targets """

        fig, (ax_top, ax_side) = plt.subplots(1, 2, figsize=(13, 5))
        target_x = [target['x'] for target in self.targets]
        target_y = [target['y'] for target in self.targets]
        target_z = [target['z'] for target in self.targets]

        # Top view: flight path, targets and the beam footprint at the start, middle and end of the acquisition
        ax_top.plot([self.x_start, self.x_end], [0, 0], 'k-', linewidth=2, label='Flight path')
        ax_top.scatter(target_x, target_y, marker='x', color='tab:red', zorder=3, label='Targets')
        y_far = 1.2 * max(target_y)
        for x_p, alpha in zip([self.x_start, (self.x_start + self.x_end) / 2, self.x_end], [0.15, 0.3, 0.15]):
            # The beam is a cone: at ground range y the slant range is sqrt(y^2 + altitude^2), and the 3 dB
            # footprint half-width is that slant range times tan(half beamwidth)
            half_width = np.sqrt(y_far**2 + self.altitude**2) * np.tan(self.beamwidth / 2)
            ax_top.fill([x_p, x_p - half_width, x_p + half_width], [0, y_far, y_far], color='tab:blue', alpha=alpha)
            ax_top.plot(x_p, 0, 'o', color='tab:blue')
        ax_top.set_xlabel('x, along-track [m]')
        ax_top.set_ylabel('y, ground range [m]')
        ax_top.set_title('Top view (shaded: 3 dB beam at start, middle, end)')
        ax_top.legend(loc='lower right')
        ax_top.grid(True)

        # Side view: platform altitude and line of sight to the targets
        ax_side.plot(0, self.altitude, 'o', color='tab:blue', label='Platform')
        ax_side.scatter(target_y, target_z, marker='x', color='tab:red', zorder=3, label='Targets')
        for y, z in zip(target_y, target_z):
            ax_side.plot([0, y], [self.altitude, z], '--', color='gray', linewidth=0.8)
        ax_side.set_xlabel('y, ground range [m]')
        ax_side.set_ylabel('z, height [m]')
        ax_side.set_title('Side view (dashed: line of sight at closest approach)')
        ax_side.legend(loc='upper right')
        ax_side.grid(True)

        fig.tight_layout()
