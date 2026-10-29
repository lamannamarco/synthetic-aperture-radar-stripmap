""" Range Migration Algorithm (RMA, also called omega-k) for a stripmap FMCW SAR

The idea in one sentence
------------------------
In the 2D frequency (wavenumber) domain, a point target at closest-approach range R0 and along-track
position x_t has a phase that is *almost* linear, and a linear phase is exactly what a 2D inverse FFT
turns into a sharp point. The RMA applies a few multiplications and one change of variable (Stolt) to
make that phase exactly linear.

Notation
--------
    t   : fast time (inside one chirp), sampled at fs
    x   : platform position along the flight path (slow time), sampled every azimuth_step
    R(x): distance between platform and target = sqrt(R0^2 + (x - x_t)^2)
    Kr  : range wavenumber [rad/m]  = 4 pi f / c   (f = instantaneous transmitted frequency)
    Kx  : azimuth wavenumber [rad/m], the Fourier variable of x
    Ky  : the wavenumber "toward the scene", Ky = sqrt(Kr^2 - Kx^2)

Steps
-----
0. Raw data: the beat signal rx * conj(tx), see signal_generator.py
       s(t, x) = exp(-j 2pi (alpha tau t + fc tau - alpha/2 tau^2)),  tau = 2 R(x) / c

1. RVP removal (+ deskew), then a Hann window over fast time. Removes the +pi alpha tau^2 term. What is left is
       s(t, x) = exp(-j 4pi/c (fc + alpha t) R(x)) = exp(-j Kr R(x)),   Kr = 4pi (fc + alpha t) / c
   This is the key FMCW trick: the transmitted frequency sweeps linearly in time, so each fast-time
   sample "sees" the scene at a different wavenumber Kr. The dechirped data is ALREADY in the range
   frequency domain: no range FFT is needed.

2. Azimuth FFT (over x). By the principle of stationary phase, the hyperbola R(x) becomes
       S(Kx, Kr) = exp(-j sqrt(Kr^2 - Kx^2) R0 - j Kx x_t)

   (2b. Only if the platform moves during the chirp: compensate for that motion, see the method.)

3. Reference function multiply (bulk compression), at a reference range R_ref (scene center):
       S(Kx, Kr) * exp(+j sqrt(Kr^2 - Kx^2) R_ref) = exp(-j sqrt(Kr^2 - Kx^2) (R0 - R_ref) - j Kx x_t)
   A target at R0 = R_ref now has a linear phase (-Kx x_t only): it is perfectly focused.
   Targets at other ranges are not, because sqrt(Kr^2 - Kx^2) is not linear in Kr and Kx.

4. Stolt interpolation: change of variable Ky = sqrt(Kr^2 - Kx^2). Resample each Kx row so that it is
   uniformly sampled in Ky instead of Kr:
       S(Kx, Ky) = exp(-j Ky (R0 - R_ref) - j Kx x_t)
   The phase is now linear in Ky and Kx for EVERY target.

5. 2D inverse FFT: Ky -> range, Kx -> azimuth. Each target becomes a 2D sinc at (R0, x_t).

All the steps assume the "stop-and-go" approximation (platform frozen during each chirp), except step 2b.
"""

import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
from scipy.fft import fft, ifft, fftfreq, fftshift
from scipy.signal.windows import hann
from scipy.constants import speed_of_light as c

class RangeMigrationAlgorithm:


    def __init__(self, params, geometry, stop_and_go=True, R_ref=None, range_window=True):
        """ Initializing the class

        params      : Params object (radar parameters)
        geometry    : Geometry object (platform, antenna, targets)
        stop_and_go : must match what was used to generate the data. If False, step 2b is applied
        R_ref       : [m] reference slant range for step 3. Default: center of the scene
        range_window: if True, a Hann window is applied over fast time (i.e. over Kr) in step 1, to lower the range
                      sidelobes. Without it, each target has sinc sidelobes at -13 dB spreading over the image
        """

        self.params = params
        self.geometry = geometry
        self.stop_and_go = stop_and_go
        self.R_ref = R_ref if R_ref is not None else (geometry.R0_min + geometry.R0_max) / 2
        self.window = hann(params.N_ADC, sym=False).astype(np.float32) if range_window else np.ones(params.N_ADC, dtype=np.float32)

        # --- Range wavenumber: one value per fast-time sample (step 1) ---
        self.Kr = 4 * np.pi * (params.fc + params.alpha * params.t_full) / c  # [rad/m]
        self.dKr = self.Kr[1] - self.Kr[0]  # [rad/m] Spacing between two fast-time samples

        # --- Azimuth wavenumber: one value per chirp (step 2) ---
        # Kept in the FFT order (0, positive, then negative values), like the output of fft()
        self.Kx = 2 * np.pi * fftfreq(geometry.Nc, d=geometry.azimuth_step)  # [rad/m]

        # --- Ky grid for the Stolt interpolation (step 4) ---
        # A target seen at an angle theta from boresight ends up at Kx = Kr sin(theta), Ky = Kr cos(theta).
        # Over the illuminated cone, Ky therefore spans from Kr_min * cos(illumination_angle) up to Kr_max.
        # We sample this span uniformly with as many points as the original range samples
        Ky_min = self.Kr[0] * np.cos(geometry.illumination_angle)
        self.Ky = np.linspace(Ky_min, self.Kr[-1], params.N_ADC)  # [rad/m]
        self.dKy = self.Ky[1] - self.Ky[0]

        # --- Image axes (step 5) ---
        self.x_axis = geometry.x_platform  # [m] Azimuth axis: one pixel per platform position
        self.r_axis = self.range_axis(self.dKy)  # [m] Range axis after the Stolt interpolation


    def range_axis(self, dK):
        """ Slant range axis [m] after an FFT over a wavenumber axis with spacing dK, centered at R_ref

        A phase -K * dR, sampled every dK, is a "tone" whose inverse FFT peaks at dR = m * 2pi / (N dK), m = bin index.
        Bins in the upper half of the FFT correspond to negative dR (targets closer than R_ref): fftfreq()
        returns exactly these values, and fftshift() puts them in increasing order
        """

        return self.R_ref + fftshift(fftfreq(self.params.N_ADC, d=dK / (2 * np.pi)))


    # ------------------------------------------------------------------------------------------------
    # Processing steps
    # ------------------------------------------------------------------------------------------------

    def range_compress(self, data, apply_window=True):
        """ Initial image: range compression over fast time, as in a conventional FMCW radar (not part of the RMA)

        With the rx * conj(tx) convention the beat frequency is negative, -alpha tau. The inverse FFT (kernel
        exp(+j 2pi f t)) puts each target at the positive bin k = alpha tau N_ADC / fs, i.e. at range k * R_max / N_ADC
        """

        window = self.window if apply_window else 1.0
        image = ifft(data * window, axis=1)
        r_axis = np.arange(self.params.N_ADC) * self.params.R_max / self.params.N_ADC
        return image, r_axis


    def remove_rvp(self, data):
        """ Step 1: residual video phase (RVP) removal and deskew, then a Hann window over fast time

        The RVP, +pi alpha tau^2, depends on the delay tau, so it differs from target to target and cannot be
        removed directly in the time domain. But in the beat-frequency domain each target sits at its own
        frequency f = -alpha tau, so there the RVP is a known function of f:  +pi alpha tau^2 = +pi f^2 / alpha.
        We multiply the spectrum by exp(-j pi f^2 / alpha) to cancel it.

        As a side effect, this filter also advances each echo in time by its own delay tau ("deskew"), so that
        all the echoes line up. To avoid circular wrap-around we zero-pad to twice the length, then keep the
        first N_ADC samples.

        Finally, a Hann window over fast time (i.e. over Kr) lowers the range sidelobes of the final image.
        """

        params = self.params
        N = params.N_ADC

        spectrum = fft(data, n=2 * N, axis=1)  # Beat-frequency domain (zero-padded)

        # Beat frequency of each FFT bin. With the rx * conj(tx) convention all the beat frequencies are negative,
        # between -fs and 0 (complex I/Q samples). Bin k holds k * fs / (2N), which is the same as k * fs / (2N) - fs
        # (frequencies are periodic with period fs): we take the negative version, except for bin 0
        f = np.arange(2 * N) * params.fs / (2 * N)  # [Hz]
        f[1:] -= params.fs

        rvp_filter = np.exp(-1j * np.pi * f**2 / params.alpha).astype(np.complex64)
        data_no_rvp = ifft(spectrum * rvp_filter, axis=1)[:, :N]  # Back to fast time

        return data_no_rvp * self.window


    def azimuth_fft(self, data):
        """ Step 2: FFT along the flight path (over the chirps), from x to Kx """

        return fft(data, axis=0)


    def compensate_motion_during_chirp(self, spectrum):
        """ Step 2b: only needed if the platform moves during each chirp (stop_and_go=False)

        During a chirp the platform moves by v*t, so the sample at fast time t was really taken at position
        x + v*t instead of x. A shift in x by v*t is, after the azimuth FFT, a phase exp(+j Kx v t).
        Multiplying by exp(-j Kx v t) undoes it, and we are back to the stop-and-go case.
        """

        v = self.geometry.platform_speed
        t = self.params.t_full
        correction = np.exp(-1j * self.Kx[:, np.newaxis] * v * t[np.newaxis, :]).astype(np.complex64)
        return spectrum * correction


    def reference_function_multiply(self, spectrum):
        """ Step 3: bulk compression at the reference range R_ref

        exp(+j sqrt(Kr^2 - Kx^2) R_ref) cancels the phase of a target at R_ref exactly. For all the other
        targets, it leaves a residual phase -sqrt(Kr^2 - Kx^2) (R0 - R_ref), which is small because R0 - R_ref
        is small: that is why the next step can get away with simple linear interpolation.
        """

        Kr = self.Kr[np.newaxis, :]  # Shape (1, N_ADC): Kr changes along the columns
        Kx = self.Kx[:, np.newaxis]  # Shape (Nc, 1): Kx changes along the rows
        Ky = np.sqrt(Kr**2 - Kx**2)  # Shape (Nc, N_ADC)
        reference = np.exp(+1j * Ky * self.R_ref).astype(np.complex64)
        return spectrum * reference


    def stolt_interpolation(self, spectrum):
        """ Step 4: resample each Kx row from a uniform Kr grid to a uniform Ky grid

        For a given Kx, the sample at Kr really belongs to Ky = sqrt(Kr^2 - Kx^2). These Ky values are not
        uniformly spaced (and they are smaller than Kr, more so for large |Kx|), so we interpolate the data on
        the uniform grid self.Ky. Outside the measured Ky values there is no data: we set it to zero.
        """

        stolt = np.zeros((self.geometry.Nc, self.params.N_ADC), dtype=np.complex64)

        for i, Kx in enumerate(self.Kx):  # One row (one Kx value) at a time
            Ky_measured = np.sqrt(self.Kr**2 - Kx**2)  # Where the samples of this row actually are, in Ky
            # np.interp only works with real numbers: interpolate the real and imaginary parts separately
            real_part = np.interp(self.Ky, Ky_measured, spectrum[i, :].real, left=0, right=0)
            imag_part = np.interp(self.Ky, Ky_measured, spectrum[i, :].imag, left=0, right=0)
            stolt[i, :] = real_part + 1j * imag_part

        return stolt


    def to_image(self, spectrum):
        """ Step 5: 2D inverse FFT, from the wavenumber domain (Kx, K) to the image domain (azimuth, range)

        Range: the phase is -K (R0 - R_ref), and the inverse FFT over K turns it into a peak at R0 - R_ref.
        Azimuth: the phase is -Kx x_t, and the inverse FFT over Kx turns it into a peak at x_t.
        """

        image = ifft(spectrum, axis=1)  # Range
        image = fftshift(image, axes=1)  # Negative R0 - R_ref first, see range_axis()
        image = ifft(image, axis=0)  # Azimuth
        return image


    def focus(self, raw_data, plot=True):
        """ Running the whole RMA on the raw data and (optionally) plotting the intermediate results """

        print('\nRMA: step 1/5, RVP removal')
        data = self.remove_rvp(raw_data)

        print('RMA: step 2/5, azimuth FFT')
        spectrum = self.azimuth_fft(data)
        if not self.stop_and_go:
            print('RMA: step 2b,  compensation of the platform motion during each chirp')
            spectrum = self.compensate_motion_during_chirp(spectrum)

        print('RMA: step 3/5, reference function multiply')
        spectrum_ref = self.reference_function_multiply(spectrum)

        print('RMA: step 4/5, Stolt interpolation')
        spectrum_stolt = self.stolt_interpolation(spectrum_ref)

        print('RMA: step 5/5, 2D inverse FFT')
        image = self.to_image(spectrum_stolt)

        if plot:
            self.plot_steps(raw_data, data, spectrum, spectrum_ref, spectrum_stolt, image)

        self.print_target_table(image)
        return image


    # ------------------------------------------------------------------------------------------------
    # Figures and checks
    # ------------------------------------------------------------------------------------------------

    def to_db(self, image):
        """ Magnitude in dB, normalized so that the strongest pixel is 0 dB """

        image_db = 20 * np.log10(np.abs(image) + 1e-12)
        return image_db - image_db.max()


    def show(self, ax, image, r_axis, x_axis, r_lim, x_lim, title, vmin=-50):
        """ Showing the part of an (azimuth, range) image inside r_lim and x_lim, in dB """

        keep_r = (r_axis >= r_lim[0]) & (r_axis <= r_lim[1])
        keep_x = (x_axis >= x_lim[0]) & (x_axis <= x_lim[1])
        crop = self.to_db(image)[np.ix_(keep_x, keep_r)]
        im = ax.imshow(crop, aspect='auto', origin='lower', cmap='viridis', vmin=vmin, vmax=0,
                       extent=[r_axis[keep_r][0], r_axis[keep_r][-1], x_axis[keep_x][0], x_axis[keep_x][-1]])
        ax.set_xlabel('Slant range [m]')
        ax.set_ylabel('Azimuth x [m]')
        ax.set_title(title)
        return im


    def scene_limits(self, margin):
        """ Range and azimuth limits [m] that contain all the targets, plus a margin """

        target_x = [target['x'] for target in self.geometry.targets]
        r_lim = (self.geometry.R0_min - margin, self.geometry.R0_max + margin)
        x_lim = (min(target_x) - margin, max(target_x) + margin)
        return r_lim, x_lim


    def show_spectrum(self, ax, spectrum, K_axis, K_label, title):
        """ Showing the magnitude of a 2D spectrum (rows: Kx in FFT order, columns: K_axis), in dB """

        Kx_shifted = fftshift(self.Kx)  # From FFT order to increasing order
        Kx_max = 1.5 * self.Kr[-1] * np.sin(self.geometry.illumination_angle)  # Show a bit beyond the signal
        keep = np.abs(Kx_shifted) <= Kx_max
        spectrum_db = self.to_db(fftshift(spectrum, axes=0)[keep, :])
        im = ax.imshow(spectrum_db, aspect='auto', origin='lower', cmap='viridis', vmin=-40, vmax=0,
                       extent=[K_axis[0], K_axis[-1], Kx_shifted[keep][0], Kx_shifted[keep][-1]])
        ax.set_xlim(self.Ky[0], self.Kr[-1])  # Same horizontal scale for every spectrum
        ax.set_xlabel(f'{K_label} [rad/m]')
        ax.set_ylabel('Kx [rad/m]')
        ax.set_title(title)
        return im


    def plot_steps(self, raw_data, data_step1, spectrum, spectrum_ref, spectrum_stolt, image):
        """ One figure with the input of the RMA and the output of each step

        Steps 1 and 3 are not directly visible: step 1 changes the data in fast time, so we look at it after a
        range FFT, and step 3 changes only the phase of the spectrum, so we look at it after the 2D inverse
        transform of step 5 (skipping the Stolt interpolation)
        """

        fig, axes = plt.subplots(2, 3, figsize=(19, 11))
        r_lim_raw = (self.geometry.R0_min - 20, self.geometry.R_max_scene + 20)  # The hyperbolas extend to R_max_scene
        x_lim_raw = (self.x_axis[0], self.x_axis[-1])  # The whole flight path
        r_lim_scene, x_lim_scene = self.scene_limits(margin=10)

        # Input: what a conventional FMCW radar would show
        initial_image, r_axis_raw = self.range_compress(raw_data)
        im = self.show(axes[0, 0], initial_image, r_axis_raw, self.x_axis, r_lim_raw, x_lim_raw,
                       '(a) Input: range-compressed raw data\nEach target traces a hyperbola (range migration)')
        axes[0, 0].set_ylabel('Platform position x [m]')
        fig.colorbar(im, ax=axes[0, 0], label='[dB]')

        # Step 1, looked at after range compression (the Hann window was already applied in step 1)
        step1_image, _ = self.range_compress(data_step1, apply_window=False)
        im = self.show(axes[0, 1], step1_image, r_axis_raw, self.x_axis, r_lim_raw, x_lim_raw,
                       '(b) Step 1: RVP removal + deskew\n'
                       r'Removes $+\pi\alpha\tau^2$, data = exp(-j Kr R(x)), Kr = 4$\pi$(fc + $\alpha$t)/c')
        axes[0, 1].set_ylabel('Platform position x [m]')
        fig.colorbar(im, ax=axes[0, 1], label='[dB]')

        # Step 2 (and 2b)
        title_step2 = '(c) Step 2: azimuth FFT (x → Kx)' + ('' if self.stop_and_go else ' + 2b: motion during chirp')
        im = self.show_spectrum(axes[0, 2], spectrum, self.Kr, 'Kr',
                                title_step2 + '\nEach column (Kr) is one frequency of the chirp')
        fig.colorbar(im, ax=axes[0, 2], label='[dB]')

        # Step 3, looked at after the 2D inverse transform (no Stolt): only the targets at R_ref are focused
        image_bulk = self.to_image(spectrum_ref)
        r_axis_bulk = self.range_axis(self.dKr)  # Without Stolt, the range axis comes from the Kr spacing
        im = self.show(axes[1, 0], image_bulk, r_axis_bulk, self.x_axis, r_lim_scene, x_lim_scene,
                       f'(d) Step 3: reference function at R_ref = {self.R_ref:.1f} m (dotted)\n'
                       'Shown after a 2D IFFT: only targets at R_ref are focused')
        axes[1, 0].axvline(self.R_ref, color='w', linestyle=':', linewidth=1)
        fig.colorbar(im, ax=axes[1, 0], label='[dB]')

        # Step 4
        im = self.show_spectrum(axes[1, 1], spectrum_stolt, self.Ky, 'Ky',
                                '(e) Step 4: Stolt interpolation, Kr → Ky = sqrt(Kr² - Kx²)\n'
                                'The spectrum bends, the phase becomes linear for all targets')
        fig.colorbar(im, ax=axes[1, 1], label='[dB]')

        # Step 5
        im = self.show(axes[1, 2], image, self.r_axis, self.x_axis, r_lim_scene, x_lim_scene,
                       '(f) Step 5: 2D inverse FFT (Ky → range, Kx → azimuth)\nFinal image: every target is focused')
        fig.colorbar(im, ax=axes[1, 2], label='[dB]')

        fig.suptitle('Range Migration Algorithm, step by step', fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.97))  # Leave room at the top for the suptitle


    def print_target_table(self, image, half_size_r=4.0, half_size_x=1.0):
        """ Comparing, for each target, the true position with the position of the brightest pixel nearby """

        image_db = self.to_db(image)
        data = []
        for target in self.geometry.targets:
            R0 = self.geometry.closest_approach_range(target)
            keep_r = np.flatnonzero(np.abs(self.r_axis - R0) <= half_size_r)
            keep_x = np.flatnonzero(np.abs(self.x_axis - target['x']) <= half_size_x)
            crop = image_db[np.ix_(keep_x, keep_r)]
            i_x, i_r = np.unravel_index(np.argmax(crop), crop.shape)
            r_peak = self.r_axis[keep_r[i_r]]
            x_peak = self.x_axis[keep_x[i_x]]
            data.append([f'{target["x"]:.3f}', f'{x_peak:.3f}', f'{R0:.3f}', f'{r_peak:.3f}', f'{crop.max():.2f}'])

        headers = ['x true [m]', 'x peak [m]', 'R0 true [m]', 'R0 peak [m]', 'Peak [dB]']
        print('\nRMA: TARGETS (peak = brightest pixel; pixel size '
              f'{abs(self.r_axis[1] - self.r_axis[0]):.3f} m in range, {self.geometry.azimuth_step*1e3:.1f} mm in azimuth)')
        print(tabulate(data, headers=headers, disable_numparse=True))
