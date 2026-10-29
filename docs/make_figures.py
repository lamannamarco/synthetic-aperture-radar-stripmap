""" Generating the figures used in README.md

Run from the repository root:  python docs/make_figures.py
The figures are saved in docs/figures/
"""

# %% Preamble
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # No window: we only save the figures
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light as c
from scipy.fft import fftshift

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from params import Params
from geometry import Geometry
from signal_generator import SignalGenerator
from rma import RangeMigrationAlgorithm

FIGURES_DIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# Colors
BLUE = '#2a78d6'
ORANGE = '#eb6834'
GRAY = '#8a8984'


# %% Same setup as main.py
params = Params(fc=10e9, fs=10e6, R_max=3e3, bandwidth=150e6)
targets = [
    {'x':   0.0, 'y': 800.0, 'z': 0.0, 'amplitude': 1.0},
    {'x': -15.0, 'y': 800.0, 'z': 0.0, 'amplitude': 1.0},
    {'x':  15.0, 'y': 800.0, 'z': 0.0, 'amplitude': 1.0},
    {'x':   0.0, 'y': 780.0, 'z': 0.0, 'amplitude': 1.0},
    {'x':   0.0, 'y': 820.0, 'z': 0.0, 'amplitude': 1.0},
]
geometry = Geometry(params, targets, beamwidth='medium', platform_speed=100.0, altitude=300.0)


# %% Figure: geometry
geometry.plot_geometry()
plt.gcf().savefig(os.path.join(FIGURES_DIR, 'geometry.png'), dpi=110)
plt.close('all')


# %% Figure: the RMA step by step (real simulated data)
raw_data = SignalGenerator(params, geometry, stop_and_go=True).generate()
rma = RangeMigrationAlgorithm(params, geometry, stop_and_go=True)
rma.focus(raw_data, plot=True)
plt.gcf().savefig(os.path.join(FIGURES_DIR, 'rma_steps.png'), dpi=90)
plt.close('all')


# %% Figure: why Stolt is needed, seen on the phase of one target
# A single target at x_t = 0 and R0 = R_ref + dR. After the reference function (step 3) its spectrum is
#     exp(-j sqrt(Kr^2 - Kx^2) dR)              (before Stolt)
# and after the Stolt change of variable
#     exp(-j Ky dR)                             (after Stolt)
# Only cos(phase) is shown, which does not depend on the sign
dR = 19.0  # [m] Same offset as the near/far targets of the simulation
Kr = np.linspace(rma.Kr[0], rma.Kr[-1], 600)
Ky = np.linspace(rma.Ky[0], rma.Ky[-1], 600)
Kx_max = rma.Kr[-1] * np.sin(geometry.illumination_angle)
Kx = np.linspace(-Kx_max, Kx_max, 400)

KR, KX = np.meshgrid(Kr, Kx)
phase_before = -np.sqrt(KR**2 - KX**2) * dR

KY, KX2 = np.meshgrid(Ky, Kx)
phase_after = -KY * dR
KR_of_sample = np.sqrt(KY**2 + KX2**2)  # Which Kr each (Ky, Kx) point came from
phase_after[(KR_of_sample < rma.Kr[0]) | (KR_of_sample > rma.Kr[-1])] = np.nan  # No measurement there

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
for ax, phase, K, label, title in [
        (axes[0], phase_before, Kr, 'Kr', 'What we have (pre-Stolt)'),
        (axes[1], phase_after, Ky, 'Ky', 'What we want (post-Stolt)')]:
    im = ax.imshow(np.cos(phase), aspect='auto', origin='lower', cmap='RdBu_r', vmin=-1, vmax=1,
                   extent=[K[0], K[-1], Kx[0], Kx[-1]])
    ax.set_xlim(rma.Ky[0], rma.Kr[-1])
    ax.set_xlabel(f'{label} [rad/m]')
    ax.set_title(title, fontsize=10)
axes[0].set_ylabel('Kx [rad/m]')
fig.colorbar(im, ax=axes, label='cos(phase)')
fig.savefig(os.path.join(FIGURES_DIR, 'stolt_phase.png'), dpi=110, bbox_inches='tight')
plt.close('all')


# %% Figure: what the Stolt interpolation does to the samples
# Wide beam (12 deg) to make the bending easier to see. Only a few samples are drawn
wide = np.deg2rad(12.0)
Kr_min, Kr_max = rma.Kr[0], rma.Kr[-1]
Kr_samples = np.linspace(Kr_min, Kr_max, 9)
Kx_rows = np.linspace(-Kr_max * np.sin(wide), Kr_max * np.sin(wide), 9)
highlight = -1  # Row drawn in orange (largest Kx)

fig, axes = plt.subplots(1, 3, figsize=(17, 5), gridspec_kw={'width_ratios': [1, 1, 1.2]})

# (a) Before: a rectangular grid, uniform in Kr
ax = axes[0]
for i, Kx_row in enumerate(Kx_rows):
    color = ORANGE if i == len(Kx_rows) + highlight else BLUE
    ax.plot(Kr_samples, np.full_like(Kr_samples, Kx_row), 'o', color=color, markersize=6)
ax.set_xlabel('Kr [rad/m]')
ax.set_ylabel('Kx [rad/m]')
ax.set_title('(a) Measured samples, after step 3\nUniform grid in (Kr, Kx)', fontsize=10)

# (b) After: each sample moves to Ky = sqrt(Kr^2 - Kx^2), along a circle of radius Kr
ax = axes[1]
angles = np.linspace(-1.15 * wide, 1.15 * wide, 200)
for Kr_value, style in [(Kr_min, '--'), (Kr_max, '-')]:  # The two circles bounding the measured band
    ax.plot(Kr_value * np.cos(angles), Kr_value * np.sin(angles), style, color=GRAY, linewidth=1,
            label=f'Circle Kr = {Kr_value:.1f} rad/m')
for i, Kx_row in enumerate(Kx_rows):
    color = ORANGE if i == len(Kx_rows) + highlight else BLUE
    ax.plot(np.sqrt(Kr_samples**2 - Kx_row**2), np.full_like(Kr_samples, Kx_row), 'o', color=color, markersize=6)
ax.legend(loc='center left', fontsize=8)
ax.set_xlabel('Ky [rad/m]')
ax.set_title('(b) Same samples at their true (Ky, Kx)\nThey fill an annular sector: not a grid anymore', fontsize=10)
ax.set_xlim(axes[0].get_xlim()[0] - 9.5, axes[0].get_xlim()[1])
axes[0].set_xlim(ax.get_xlim())
for a in axes[:2]:
    a.set_ylim(-1.2 * Kr_max * np.sin(wide), 1.2 * Kr_max * np.sin(wide))

# (c) One row: measured samples (non-uniform in Ky) are interpolated onto the uniform Ky grid
ax = axes[2]
Kx_row = Kx_rows[highlight]
Kr_fine = np.linspace(Kr_min, Kr_max, 40)  # A coarse version of the real 3003 samples
Ky_measured = np.sqrt(Kr_fine**2 - Kx_row**2)
Ky_grid = np.linspace(Kr_min * np.cos(wide), Kr_max, 30)
signal = lambda Ky_value: np.cos(-Ky_value * 2.0)  # Real part of exp(-j Ky dR), dR = 2 m to keep it readable
Ky_curve = np.linspace(Ky_measured[0], Ky_measured[-1], 800)  # Only where this row was measured
ax.plot(Ky_curve, signal(Ky_curve), '-', color=GRAY, linewidth=1, label='True spectrum (unknown)')
ax.plot(Ky_measured, signal(Ky_measured), 'o', color=ORANGE, markersize=7, label='Measured samples')
interpolated = np.interp(Ky_grid, Ky_measured, signal(Ky_measured), left=0, right=0)
ax.plot(Ky_grid, interpolated, 'x', color=BLUE, markersize=8, markeredgewidth=2, label='Interpolated on uniform Ky grid')
for Ky_value in Ky_grid:
    ax.axvline(Ky_value, color=GRAY, linewidth=0.4, alpha=0.5)
ax.set_xlabel('Ky [rad/m]')
ax.set_ylabel('Re{spectrum}')
ax.set_title(f'(c) Row Kx = {Kx_row:.1f} rad/m (orange in a, b)\nResample on the uniform Ky grid, zero where there is no data', fontsize=10)
ax.legend(loc='upper right', fontsize=8)

fig.suptitle('Stolt interpolation (wide 12° beam, few samples, for visibility)', fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'stolt_mapping.png'), dpi=110)
plt.close('all')

# %% Figure: 3D view of the acquisition (Section 2.2)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

h = geometry.altitude
half_beam = np.deg2rad(12.0) / 2  # Wide preset, to make the footprint easier to see
y_near, y_far = 650.0, 950.0  # Ground range covered by the elevation beam (the swath)
x_min, x_max = -300.0, 300.0  # Part of the flight path drawn
x_p = 0.0  # Platform position shown

fig = plt.figure(figsize=(12, 7.5))
ax = fig.add_subplot(111, projection='3d')

# Imaged strip on the ground: the swath, swept along x as the platform flies
strip = [[(x_min, y_near, 0), (x_max, y_near, 0), (x_max, y_far, 0), (x_min, y_far, 0)]]
ax.add_collection3d(Poly3DCollection(strip, facecolor=GRAY, alpha=0.12, edgecolor='none'))
ax.text(x_max - 40, y_far - 10, 0, 'Imaged strip (swath)', color='#52514e', fontsize=9)

# Flight path, a few past/future platform positions, and the direction of flight
ax.plot([x_min, x_max], [0, 0], [h, h], 'k-', linewidth=2)
for x_ghost in np.arange(x_min + 50, x_max, 100):
    if x_ghost != x_p:
        ax.scatter(x_ghost, 0, h, color=GRAY, s=15)
ax.quiver(x_max - 130, 0, h + 30, 110, 0, 0, color='k', arrow_length_ratio=0.25, linewidth=1.5)
ax.text(x_max - 150, 0, h + 50, 'Flight direction, speed v', fontsize=9)
ax.scatter(x_p, 0, h, color=BLUE, s=60, depthshade=False)
ax.text(x_p - 60, 0, h + 30, 'Radar', color=BLUE, fontsize=10)

# Beam: from the antenna to its footprint on the ground. The azimuth half-width grows with range
def half_width(y):
    return np.sqrt(y**2 + h**2) * np.tan(half_beam)
corners = [(x_p - half_width(y_near), y_near, 0), (x_p + half_width(y_near), y_near, 0),
           (x_p + half_width(y_far), y_far, 0), (x_p - half_width(y_far), y_far, 0)]
ax.add_collection3d(Poly3DCollection([corners], facecolor=BLUE, alpha=0.3, edgecolor=BLUE))
apex = (x_p, 0, h)
for k in range(4):  # Sides of the beam
    side = [apex, corners[k], corners[(k + 1) % 4]]
    ax.add_collection3d(Poly3DCollection([side], facecolor=BLUE, alpha=0.06, edgecolor=BLUE, linewidth=0.5))
ax.text(x_p + half_width(y_far) + 15, y_far - 20, 0, 'Beam footprint', color=BLUE, fontsize=9)

# Targets
for target in targets:
    ax.scatter(target['x'], target['y'], target['z'], color=ORANGE, marker='x', s=40, depthshade=False)
ax.text(-60, 740, 0, 'Targets', color=ORANGE, fontsize=9, ha='right')

# Slant range, altitude and ground range to the scene center
y_c = 800.0
ax.plot([x_p, 0], [0, y_c], [h, 0], '--', color=ORANGE, linewidth=1.5)
ax.text(15, y_c / 2, h / 2 + 15, 'Slant range $R_0$', color=ORANGE, fontsize=9)
ax.plot([x_p, x_p], [0, 0], [h, 0], ':', color='#52514e', linewidth=1.2)
ax.text(x_p + 10, 0, h / 2, 'Altitude h', color='#52514e', fontsize=9)
ax.plot([x_p, 0], [0, y_c], [0, 0], ':', color='#52514e', linewidth=1.2)
ax.text(x_p + 10, y_c / 2, 0, 'Ground range', color='#52514e', fontsize=9)

ax.set_xlim(x_min, x_max)
ax.set_ylim(0, 1000)
ax.set_zlim(0, 380)
ax.set_box_aspect((1.3, 2.0, 1.0))
ax.set_xlabel('x, along-track [m]')
ax.set_ylabel('y, ground range [m]')
ax.set_zlabel('z, height [m]', labelpad=2)
ax.view_init(elev=20, azim=-38)
ax.set_title('Stripmap acquisition: the radar flies along x and looks sideways (along y)', fontsize=11)

fig.subplots_adjust(left=0, right=0.95, bottom=0.03, top=0.95)
fig.savefig(os.path.join(FIGURES_DIR, 'acquisition_3d.png'), dpi=110, bbox_inches='tight', pad_inches=0.3)
plt.close('all')


# %% Figure: the synthetic aperture (Section 2.2)
# Five platform positions, 0 to 4, while the target crosses the beam. At each position, theta_i is the look angle:
# the angle between broadside (dotted) and the line of sight to the target. The line of sight is drawn thicker
# the closer the target is to the beam center (i.e. the more strongly it is illuminated)
fig, ax = plt.subplots(figsize=(11, 6.2))
half_beam = np.deg2rad(28)  # Beam made wide on purpose, to make the picture readable
y_t = 8.0
L = 2 * y_t * np.tan(half_beam)  # Distance flown while the target is inside the beam
y_top = y_t + 2.5
positions = np.linspace(-L / 2, L / 2, 5)

ax.plot([-11, 11], [0, 0], 'k-', linewidth=2)
ax.text(-10.8, 0.3, 'Flight path (x)', fontsize=9)
for i, x_p in enumerate(positions):
    # Beam
    edge = y_top * np.tan(half_beam)
    ax.fill([x_p, x_p - edge, x_p + edge], [0, y_top, y_top], color=BLUE, alpha=0.07)
    ax.plot([x_p, x_p - edge], [0, y_top], color=BLUE, linewidth=0.7, alpha=0.6)
    ax.plot([x_p, x_p + edge], [0, y_top], color=BLUE, linewidth=0.7, alpha=0.6)

    # Line of sight, thicker when the target is closer to the beam center
    look_angle = np.arctan2(-x_p, y_t)  # Angle of the target from broadside
    gain = np.exp(-4 * np.log(2) * (look_angle / (2 * half_beam))**2)  # Gaussian beam, as in geometry.py
    ax.plot([x_p, 0], [0, y_t], '-', color=ORANGE, linewidth=0.5 + 3.5 * gain, alpha=0.9)

    # Broadside direction and look angle theta_i
    ax.plot([x_p, x_p], [0, 3.6], ':', color='#52514e', linewidth=1)
    if abs(look_angle) > 1e-9:
        los_deg = np.rad2deg(np.arctan2(y_t, -x_p))  # Direction of the line of sight, from the +x axis
        ax.add_patch(matplotlib.patches.Arc((x_p, 0), 5.0, 5.0, theta1=min(90, los_deg), theta2=max(90, los_deg),
                                            color='#52514e', linewidth=1.2))
        mid = np.deg2rad((90 + los_deg) / 2)
        ax.text(x_p + 3.1 * np.cos(mid), 3.1 * np.sin(mid), f'$\\theta_{i}$', ha='center', va='center', fontsize=12)
    else:
        ax.text(x_p, 4.1, f'$\\theta_{i} = 0$', ha='center', va='center', fontsize=12,
                bbox=dict(facecolor='white', edgecolor='none', pad=1.5))

    ax.plot(x_p, 0, 'o', color=BLUE, markersize=9)
    ax.text(x_p, -0.5, f'{i}', ha='center', va='top', fontsize=12, fontweight='bold')

ax.plot(0, y_t, 'x', color=ORANGE, markersize=14, markeredgewidth=3, zorder=5)
ax.text(0.4, y_t + 0.35, 'Target', color=ORANGE, fontsize=10)
ax.set_xlim(-11, 11)
ax.set_ylim(-1.6, y_top + 0.3)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('Positions 0 to 4 see the target at look angles $\\theta_0$ to $\\theta_4$ (from broadside, dotted)\n'
             'Orange: line of sight, thicker = target closer to the beam center', fontsize=10)

fig.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'synthetic_aperture.png'), dpi=110)
plt.close('all')


# %% Figure: range migration (Section 2.4)
fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={'width_ratios': [1.1, 1]})

# (a) Sketch: lines of sight from several platform positions to one target
ax = axes[0]
y_t = 10.0
ax.plot([-9, 9], [0, 0], 'k-', linewidth=2)
ax.plot(0, y_t, 'x', color=ORANGE, markersize=12, markeredgewidth=3)
ax.text(0.4, y_t + 0.35, 'Target', color=ORANGE, fontsize=10)
for x_p in [-7, -3.5, 0, 3.5, 7]:
    R = np.hypot(x_p, y_t)
    ax.plot([x_p, 0], [0, y_t], '--' if x_p != 0 else '-', color=BLUE, linewidth=1.2)
    ax.plot(x_p, 0, 'o', color=BLUE, markersize=8)
    ax.text(x_p, -0.6, f'R = {R:.1f}', ha='center', va='top', fontsize=9)
ax.text(0.3, y_t / 2, 'R₀ (closest\napproach)', fontsize=9, color='#52514e')
ax.annotate('', xy=(8.5, -2.0), xytext=(-8.5, -2.0), arrowprops=dict(arrowstyle='-|>', color='#52514e'))
ax.text(0, -2.3, 'Platform moving along x', ha='center', va='top', fontsize=9)
ax.set_xlim(-9.5, 9.5)
ax.set_ylim(-3.5, y_t + 1.5)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('(a) The distance R(x) first decreases, then increases\n(sketch, arbitrary units)', fontsize=10)

# (b) Real numbers: R(x) for a target at R0 = 854 m, medium beam, against the 1 m range bins
ax = axes[1]
R0 = geometry.closest_approach_range(targets[0])
x_half = geometry.half_illumination_length
x_p = np.linspace(-x_half, x_half, 1001)
R = np.sqrt(R0**2 + x_p**2)
bin_size = c / (2 * params.bandwidth)
for k in range(-1, 7):  # Range bin edges around R0
    ax.axvline(R0 + (k - 0.5) * bin_size, color=GRAY, linewidth=0.6, alpha=0.6)
ax.plot(R, x_p, color=BLUE, linewidth=2.5, label='R(x) = √(R₀² + x²)')
ax.plot(R0, 0, 'o', color=ORANGE, markersize=8, label=f'Closest approach R₀ = {R0:.1f} m')
ax.text(R0 + 5.6 * bin_size, 0.2 * x_half, f'Over the illuminated\nflight path, R changes\nby {R.max() - R0:.1f} m, i.e. about\n{(R.max() - R0) / bin_size:.0f} range bins',
        ha='right', fontsize=9, color='#52514e')
ax.set_xlim(R0 - 1.5 * bin_size, R0 + 6.5 * bin_size)
ax.set_xlabel('Slant range [m]  (thin lines: 1 m range bins)')
ax.set_ylabel('Platform position x relative to the target [m]')
ax.set_title('(b) After a range FFT, the target traces a hyperbola\nacross several range bins (medium beam)', fontsize=10)
ax.legend(loc='upper left', fontsize=9)

fig.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'range_migration.png'), dpi=110)
plt.close('all')


# %% Figure: one target through all the domains, from hyperbola to point (Section 2.5)
# A single target 19 m beyond R_ref (R_ref = scene center of main.py), so that steps 3 and 4 make a visible difference
one_target = [{'x': 0.0, 'y': 820.0, 'z': 0.0, 'amplitude': 1.0}]
geometry_one = Geometry(params, one_target, beamwidth='medium', platform_speed=100.0, altitude=300.0)
raw_one = SignalGenerator(params, geometry_one, stop_and_go=True).generate()
R_ref_center = np.sqrt(800.0**2 + 300.0**2)
rma_one = RangeMigrationAlgorithm(params, geometry_one, stop_and_go=True, R_ref=R_ref_center)

# Running the steps one by one, to keep every intermediate result
range_compressed, r_axis_raw = rma_one.range_compress(raw_one)
data_one = rma_one.remove_rvp(raw_one)  # Includes the Hann window
spectrum_ref_one = rma_one.reference_function_multiply(rma_one.azimuth_fft(data_one))
spectrum_stolt_one = rma_one.stolt_interpolation(spectrum_ref_one)
image_one = rma_one.to_image(spectrum_stolt_one)

R0_one = geometry_one.closest_approach_range(one_target[0])
Kx_shifted = fftshift(rma_one.Kx)
keep_kx = np.abs(Kx_shifted) <= rma_one.Kr[-1] * np.sin(geometry_one.illumination_angle)

fig, axes = plt.subplots(1, 4, figsize=(19, 4.8))

# (a) Range-compressed data: the hyperbola
rma_one.show(axes[0], range_compressed, r_axis_raw, rma_one.x_axis, (R0_one - 4, R0_one + 9),
             (rma_one.x_axis[0], rma_one.x_axis[-1]), '(a) Raw data after range compression\ndomain (x, R): a hyperbola')
axes[0].set_ylabel('Platform position x [m]')

# (b) and (c) Real part of the spectrum, before and after Stolt
for ax, spectrum, K_axis, K_label, title in [
        (axes[1], spectrum_ref_one, rma_one.Kr, 'Kr', '(b) After steps 1 to 3, Re{S}\ndomain (Kx, Kr): curved fringes'),
        (axes[2], spectrum_stolt_one, rma_one.Ky, 'Ky', '(c) After step 4 (Stolt), Re{S}\ndomain (Kx, Ky): straight fringes')]:
    # The azimuth FFT counts x from the first platform position, x_start, so a target at x_t = 0 carries a
    # phase +Kx*|x_start| that wraps hundreds of times along Kx. Removing it (display only) leaves the
    # part of the phase that matters here
    display = spectrum * np.exp(-1j * rma_one.Kx[:, np.newaxis] * geometry_one.x_start)
    crop = fftshift(display, axes=0)[keep_kx, :]
    ax.imshow(crop.real / np.abs(crop).max(), aspect='auto', origin='lower', cmap='RdBu_r', vmin=-1, vmax=1,
              extent=[K_axis[0], K_axis[-1], Kx_shifted[keep_kx][0], Kx_shifted[keep_kx][-1]])
    ax.set_xlim(rma_one.Ky[0], rma_one.Kr[-1])
    ax.set_xlabel(f'{K_label} [rad/m]')
    ax.set_ylabel('Kx [rad/m]')
    ax.set_title(title, fontsize=10)

# (d) Final image, zoomed on the target
im = rma_one.show(axes[3], image_one, rma_one.r_axis, rma_one.x_axis, (R0_one - 5, R0_one + 5), (-1.0, 1.0),
                  '(d) After step 5 (2D inverse FFT)\ndomain (x, R): a point', vmin=-40)
axes[3].plot(R0_one, 0, 'w+', markersize=14)
fig.colorbar(im, ax=axes[3], label='[dB]')

fig.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'one_target_domains.png'), dpi=100)
plt.close('all')


print(f'\nFigures saved in {FIGURES_DIR}')
