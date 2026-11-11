import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from field_io import list_grid_files, load_fields, stream_function

MARKERS = {0: "o", 50: "s", 100: "^"}


def select(files, **criteria):
    chosen = {}
    for path, params in files:
        if all(float(params[key]) == value for key, value in criteria.items()):
            chosen[int(float(params["Ha"]))] = path
    return dict(sorted(chosen.items()))


def plot_stream_function(files, Ra, output_dir):
    selected = select(files, Ra=Ra, N=0.0, Pr=0.054)
    if not selected:
        return
    fig, axes = plt.subplots(1, len(selected), figsize=(4.4 * len(selected), 4.3), dpi=150, squeeze=False)
    for ax, (Ha, path) in zip(axes[0], selected.items()):
        fields = load_fields(path)
        psi = stream_function(fields["u"])
        filled = ax.contourf(fields["X"], fields["Y"], psi, levels=25, cmap="coolwarm")
        ax.contour(fields["X"], fields["Y"], psi, levels=12, colors="k", linewidths=0.4)
        ax.set_title(f"Ra = {Ra:g}, Ha = {Ha}  (|\u03c8|max = {np.abs(psi).max():.2f})")
        ax.set_xlabel("X (lattice unit)")
        ax.set_ylabel("Y (lattice unit)")
        ax.set_aspect("equal")
        fig.colorbar(filled, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(output_dir / f"streamfunction_Ra{int(Ra)}.png", bbox_inches="tight")
    plt.close(fig)


def plot_midheight_profiles(files, Ra, output_dir):
    selected = select(files, Ra=Ra, N=0.0, Pr=0.054)
    if not selected:
        return
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), dpi=150)
    for Ha, path in selected.items():
        fields = load_fields(path)
        n = fields["T"].shape[1]
        row = fields["T"].shape[0] // 2
        x = np.arange(n) / (n - 1)
        axes[0].plot(x, fields["T"][row, :], label=f"Ha = {Ha}")
        axes[1].plot(x, fields["C"][row, :], label=f"Ha = {Ha}")
        axes[2].plot(x, fields["v"][row, :], label=f"Ha = {Ha}")
    for ax, title in zip(axes, ("Temperature at mid-height", "Concentration at mid-height", "Vertical velocity v at mid-height")):
        ax.set_title(title)
        ax.set_xlabel("x / L")
        ax.grid(alpha=0.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / f"midheight_profiles_Ra{int(Ra)}.png", bbox_inches="tight")
    plt.close(fig)


def plot_transfer_rates(summary_path, output_dir):
    summary = pd.read_csv(summary_path)
    summary = summary[(summary["N"] == 0) & (summary["Pr"] == 0.054)]
    fig, (ax_nu, ax_sh) = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150)
    for Ha, marker in MARKERS.items():
        subset = summary[summary["Ha"] == Ha].sort_values("Ra")
        if subset.empty:
            continue
        ax_nu.plot(subset["Ra"], subset["avgNu_left"], marker=marker, label=f"Ha = {Ha}")
        ax_sh.plot(subset["Ra"], subset["avgSh_left"].abs(), marker=marker, label=f"Ha = {Ha}")
    for ax, title in ((ax_nu, "Average Nusselt number (hot wall)"), (ax_sh, "Average Sherwood number (left wall)")):
        ax.set_xscale("log")
        ax.set_xlabel("Rayleigh number, Ra")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend()
    ax_nu.set_ylabel("Nu")
    ax_sh.set_ylabel("Sh")
    fig.tight_layout()
    fig.savefig(output_dir / "Nu_Sh_vs_Ra.png", bbox_inches="tight")
    plt.close(fig)


def plot_hartmann_reduction(summary_path, output_dir):
    summary = pd.read_csv(summary_path)
    summary = summary[(summary["N"] == 0) & (summary["Pr"] == 0.054)]
    fig, (ax_nu, ax_sh) = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150)
    for Ra in sorted(summary["Ra"].unique()):
        subset = summary[summary["Ra"] == Ra].sort_values("Ha")
        base = subset[subset["Ha"] == 0]
        if base.empty:
            continue
        ax_nu.plot(subset["Ha"], subset["avgNu_left"] / base["avgNu_left"].iloc[0], marker="o", label=f"Ra = {Ra:g}")
        ax_sh.plot(subset["Ha"], subset["avgSh_left"] / base["avgSh_left"].iloc[0], marker="o", label=f"Ra = {Ra:g}")
    for ax, title in ((ax_nu, "Nu / Nu(Ha = 0)"), (ax_sh, "Sh / Sh(Ha = 0)")):
        ax.set_xlabel("Hartmann number, Ha")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "Nu_Sh_vs_Ha.png", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir")
    parser.add_argument("output_dir", nargs="?", default="output_plots")
    parser.add_argument("--summary", default="summary_results_Ra_Ha_study.csv")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    files = list_grid_files(args.data_dir)

    for Ra in (1000.0, 10000.0, 100000.0):
        plot_stream_function(files, Ra, output_dir)
        plot_midheight_profiles(files, Ra, output_dir)

    summary_path = Path(args.summary) if Path(args.summary).exists() else Path(args.data_dir) / args.summary
    if summary_path.exists():
        plot_transfer_rates(summary_path, output_dir)
        plot_hartmann_reduction(summary_path, output_dir)


if __name__ == "__main__":
    main()
