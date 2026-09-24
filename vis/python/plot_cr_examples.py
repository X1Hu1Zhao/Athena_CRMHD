"""Plot one-dimensional CR benchmark outputs against their reference solutions."""

import argparse
import glob
import math
from pathlib import Path

import athena_read
import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_TIMES = {
    "streaming": (0.0, 0.15, 0.30),
    "landau": (0.0, 0.05, 0.10),
    "growing": (0.0, 1.0, 2.0),
}
DEFAULT_FILES = {
    "streaming": ("CR_streaming.out2.*.athdf", "streaming_comparison.png"),
    "landau": ("CR_streaming_nll.out2.*.athdf", "landau_comparison.png"),
    "growing": ("CRMHD_growing.out2.*.athdf", "growing_comparison.png"),
}
VMAX = {"streaming": 100.0, "landau": 100.0, "growing": 1000.0}


def select_snapshots(pattern, target_times):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise ValueError(f"No HDF5 snapshots match {pattern!r}")

    times = []
    for path in paths:
        with h5py.File(path, "r") as snapshot:
            times.append(float(snapshot.attrs["Time"]))

    selected = []
    for target_time in target_times:
        nearest = int(np.argmin(np.abs(np.asarray(times) - target_time)))
        tolerance = max(1e-3, 0.02 * abs(target_time))
        if abs(times[nearest] - target_time) > tolerance:
            raise ValueError(
                f"No snapshot near t={target_time:g}; nearest is t={times[nearest]:g}. "
                "Check the run's final time and HDF5 output cadence."
            )
        if nearest in selected:
            raise ValueError("The requested times select the same snapshot")
        selected.append(nearest)

    return [(paths[index], times[index]) for index in selected]


def streaming_solution(positions, time):
    speed = 5.0 / 3.0
    plateau_edge = math.sqrt((speed * time) ** 2 + 4.0 * speed * time)
    energy = np.maximum(2.0 + speed * time - np.maximum(np.abs(positions), plateau_edge), 0.0)
    flux = np.where(np.abs(positions) > plateau_edge, np.sign(positions) * speed * energy, 0.0)
    return energy, flux


def landau_solution(positions, time):
    speed = 11.0 / 6.0
    initial_integral = 2.0 / 3.0 * math.log(16.0 / 10.0)

    def balance(edge):
        denominator = 10.0 - 3.0 * speed * time + 3.0 * edge
        return (2.0 / 3.0 * math.log(16.0 / denominator)
                + 2.0 * edge / denominator - initial_integral)

    if time == 0.0:
        plateau_edge = 0.0
    else:
        lower, upper = 0.0, max(1.0, float(np.max(np.abs(positions))))
        while balance(upper) > 0.0:
            upper *= 2.0
        for _ in range(60):
            midpoint = (lower + upper) / 2.0
            if balance(midpoint) > 0.0:
                lower = midpoint
            else:
                upper = midpoint
        plateau_edge = (lower + upper) / 2.0

    energy = 1.0 / (10.0 + 3.0 * (np.maximum(np.abs(positions), plateau_edge) - speed * time))
    flux = np.where(np.abs(positions) < plateau_edge, 0.0, np.sign(positions) * speed * energy)
    return energy, flux


def growing_solution(positions, time, domain_length):
    alfven_speed = 0.1
    sound_speed = 0.01
    vmax = VMAX["growing"]
    frequency = 0.3218887947075146 + 0.14203248351814002j
    wavenumber = 2.0 * math.pi / domain_length
    density = (10.0e-3
               * np.exp(1j * (wavenumber * positions - frequency.real * wavenumber
                              * alfven_speed * time))
               * np.exp(frequency.imag * wavenumber * alfven_speed * time))
    flux = (alfven_speed ** 3 * (3.0 * frequency + 1.0)
            * (frequency ** 2 - (sound_speed / alfven_speed) ** 2) * density
            / (1.0 - 3.0 * frequency ** 2 * alfven_speed ** 2 / vmax ** 2))
    pressure = ((1.0 + frequency * alfven_speed ** 2 / vmax ** 2) * flux
                / (alfven_speed * (3.0 * frequency + 1.0)))
    return 3.0 * pressure.real, flux.real


def load_profiles(path, mode):
    snapshot = athena_read.athdf(path)
    positions = np.asarray(snapshot["x1v"])
    built_in_fields = ("cr_1_energy_dens", "cr_1_flux1", "rho")
    if all(field in snapshot for field in built_in_fields):
        energy = snapshot["cr_1_energy_dens"]
        flux = snapshot["cr_1_flux1"] * VMAX[mode]
        density = snapshot["rho"]
    else:
        legacy_fields = ("user_out_var0", "user_out_var1", "user_out_var4")
        missing = [field for field in legacy_fields if field not in snapshot]
        if missing:
            raise ValueError(
                "Output must contain cr_1_energy_dens, cr_1_flux1, and rho; "
                f"missing legacy fallback fields: {', '.join(missing)}"
            )
        energy = snapshot["user_out_var0"]
        flux = snapshot["user_out_var1"] * VMAX[mode]
        density = snapshot["user_out_var4"]

    expected_shape = (1, 1, len(positions))
    if energy.shape != expected_shape or flux.shape != expected_shape \
            or density.shape != expected_shape:
        raise ValueError("This plotter expects one-dimensional output")
    length = float(snapshot["x1f"][-1] - snapshot["x1f"][0])
    return (positions, energy[0, 0, :], flux[0, 0, :], density[0, 0, :], length)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=DEFAULT_TIMES)
    parser.add_argument(
        "snapshots", nargs="?",
        help="Quoted .athdf glob (default: the selected example's output files)")
    parser.add_argument(
        "output", nargs="?",
        help="Output PNG filename (default: the selected example's comparison plot)")
    parser.add_argument("--times", type=float, nargs=3,
                        help="Three output times (default: benchmark times)")
    args = parser.parse_args()
    target_times = args.times if args.times is not None else DEFAULT_TIMES[args.mode]
    default_snapshots, default_output = DEFAULT_FILES[args.mode]
    snapshots = args.snapshots if args.snapshots is not None else default_snapshots
    output = args.output if args.output is not None else default_output
    selected = select_snapshots(snapshots, target_times)

    if args.mode == "growing":
        figure, axis = plt.subplots(figsize=(8, 4), constrained_layout=True)
        axes = (axis,)
    else:
        figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                   constrained_layout=True)
    colors = ("tab:blue", "tab:orange", "tab:green")
    analytical_labeled = [False] * len(axes)
    for color, (path, time) in zip(colors, selected):
        positions, energy, flux, _, length = load_profiles(path, args.mode)
        if args.mode == "growing":
            actual = (energy - np.mean(energy),)
            predicted = (growing_solution(positions, time, length)[0],)
        else:
            actual = (energy, flux)
            predicted = (streaming_solution(positions, time) if args.mode == "streaming"
                         else landau_solution(positions, time))
        if time == 0.0:
            relative_initial_error = (np.linalg.norm(actual[0] - predicted[0])
                                      / np.linalg.norm(predicted[0]))
            if relative_initial_error > 0.01:
                raise ValueError(
                    f"Initial CR energy differs from the {args.mode} reference "
                    f"by {relative_initial_error:.1%}; do not interpret this "
                    "plot as an analytical validation."
                )

        for axis_index, axis in enumerate(axes):
            axis.plot(positions, actual[axis_index], color=color, alpha=0.8,
                      label=f"Simulation t={time:.4g}")
            if args.mode == "growing" or time > 0.0 or axis_index == 0:
                axis.plot(positions, predicted[axis_index], color="black", ls="--",
                          lw=1.3, label="Analytical" if not analytical_labeled[axis_index]
                          else "_nolegend_")
                analytical_labeled[axis_index] = True

    if args.mode == "growing":
        labels = (r"$\mathcal{E}_{\rm cr}-\langle\mathcal{E}_{\rm cr}\rangle$",)
    else:
        labels = (r"$\mathcal{E}_{\rm cr}$", r"$F_{{\rm cr},1}$")
    for axis, label in zip(axes, labels):
        axis.set_ylabel(label)
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    axes[-1].set_xlabel("$x$")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=170)
    plt.close(figure)
    print(f"Saved {output}")


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, ValueError) as error:
        raise SystemExit(f"Cannot plot CR example: {error}") from error
