""" Class generating the raw (dechirped, complex I/Q) FMCW data of a stripmap SAR

FMCW signal model
-----------------
The radar transmits a chirp of length Tc whose frequency rises linearly at a rate alpha = B/Tc:
    s_tx(t) = exp(j 2pi (fc t + alpha/2 t^2)),    0 <= t < Tc

A point target at range R sends back a copy of it, delayed by the round-trip time tau = 2R/c:
    s_rx(t) = s_tx(t - tau)

The receiver mixes ("dechirps") the echo with the transmitted chirp, i.e. it multiplies the received signal by
the complex conjugate of the transmitted one. The output is called the beat (or IF) signal:
    s_beat(t) = s_rx(t) * conj(s_tx(t))
              = exp(-j 2pi (alpha tau t  +  fc tau  -  alpha/2 tau^2))
                            |                |             |
                            |                |             +-- residual video phase (RVP), small
                            |                +---------------- carrier phase: changes from chirp to chirp
                            |                                  as the platform moves -> this is what SAR uses
                            +--------------------------------- beat frequency f_b = -alpha tau (negative with this
                                                               convention): proportional to the range

The ADC samples s_beat(t) at fs, giving N_ADC complex samples per chirp ("fast time"). Repeating this
for every chirp along the flight path ("slow time") gives a 2D matrix of shape (Nc, N_ADC).
"""

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.fft import ifft
from scipy.signal.windows import hann
from scipy.constants import speed_of_light as c

class SignalGenerator:


    def __init__(self, params, geometry, stop_and_go=True, snr_db=None, seed=0):
        """ Initializing the class

        params      : Params object (radar parameters)
        geometry    : Geometry object (platform, antenna, targets)
        stop_and_go : if True, the platform is assumed frozen during each chirp (it "jumps" between chirps).
                      If False, the platform keeps moving during the chirp, which is what really happens.
                      With FMCW the chirps are long, so this can matter (more on that when we process the data)
        snr_db      : [dB] signal-to-noise ratio per sample of a target with amplitude 1. None = no noise
        seed        : seed of the random noise generator, to get repeatable results
        """

        self.params = params
        self.geometry = geometry
        self.stop_and_go = stop_and_go
        self.snr_db = snr_db
        self.rng = np.random.default_rng(seed)


    def generate(self):
        """ Generating the raw data matrix, of shape (Nc, N_ADC): one row per chirp, one column per ADC sample """

        params = self.params
        geometry = self.geometry
        data = np.zeros((geometry.Nc, params.N_ADC), dtype=np.complex64)  # complex64 halves the memory

        for n in tqdm(range(geometry.Nc), desc='Generating chirps'):

            # Time at each ADC sample of this chirp: the chirp starts at t_slow[n], then the fast time runs
            if self.stop_and_go:
                t = geometry.t_slow[n]  # One single time for the whole chirp (platform frozen)
            else:
                t = geometry.t_slow[n] + params.t_full  # One time per ADC sample (platform moving)

            # Where the platform is at that time
            x_p, y_p, z_p = geometry.platform_position(t)

            for target in geometry.targets:

                # Distance between platform and target
                dx = target['x'] - x_p
                dy = target['y'] - y_p
                dz = target['z'] - z_p
                R = np.sqrt(dx**2 + dy**2 + dz**2)

                # Azimuth angle of the target as seen from the antenna (0 = straight abeam).
                # If the target is outside the beam it does not contribute to this chirp
                angle = np.arcsin(dx / R)
                gain = geometry.antenna_gain(angle)
                if np.all(gain == 0):
                    continue

                # Round-trip delay and beat signal, rx * conj(tx) (see the module docstring for the formula).
                # We neglect the 1/R^2 propagation loss: the targets are at similar ranges anyway
                tau = 2 * R / c
                phase = -2 * np.pi * (params.alpha * tau * params.t_full + params.fc * tau - 0.5 * params.alpha * tau**2)
                data[n, :] += target['amplitude'] * gain * np.exp(1j * phase)

        # Complex white Gaussian noise: the power is split equally between I and Q
        if self.snr_db is not None:
            noise_power = 10**(-self.snr_db / 10)  # A target of amplitude 1 has power 1
            noise = np.sqrt(noise_power / 2) * (self.rng.standard_normal(data.shape) + 1j * self.rng.standard_normal(data.shape))
            data += noise.astype(np.complex64)

        return data


    def plot_range_compressed(self, data, margin=30.0):
        """ Quick look at the data after a range FFT (each row: one chirp, each column: one range bin)

        This is not SAR focusing yet: it only shows, for every platform position, at which range each target
        is. As the platform flies by, the range to a target first decreases then increases again, so each
        target traces a hyperbola ("range migration"). Focusing will collapse each hyperbola into a point.

        margin : [m] range shown before the closest and after the farthest target
        """

        params = self.params
        geometry = self.geometry

        # Range compression, chirp by chirp. The beat frequency is negative, -alpha tau, so we use the inverse FFT
        # (kernel exp(+j 2pi f t)): it puts each target at the positive bin k = alpha tau N_ADC / fs.
        # The Hann window lowers the sidelobes
        window = hann(params.N_ADC, sym=False).astype(np.float32)  # float32, so that the data stays complex64
        range_profiles = ifft(data * window, axis=1)

        # Bin k corresponds to |f_b| = k * fs / N_ADC, i.e. to range R = |f_b| * c / (2 alpha)
        f_beat = np.arange(params.N_ADC) * params.fs / params.N_ADC
        r_axis = f_beat * c / (2 * params.alpha)

        # Keep only the range bins around the targets
        keep = (r_axis >= geometry.R0_min - margin) & (r_axis <= geometry.R_max_scene + margin)
        image = 20 * np.log10(np.abs(range_profiles[:, keep]) + 1e-12)
        image -= image.max()  # Normalize: 0 dB = strongest pixel

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(image, aspect='auto', origin='lower', cmap='viridis', vmin=-50, vmax=0,
                       extent=[r_axis[keep][0], r_axis[keep][-1], geometry.x_platform[0], geometry.x_platform[-1]])
        ax.set_xlabel('Slant range [m]')
        ax.set_ylabel('Platform position x [m]')
        ax.set_title('Range-compressed raw data')
        fig.colorbar(im, ax=ax, label='[dB]')
        fig.tight_layout()
