import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from field_io import list_grid_files, load_fields

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.titlesize": 12,
    }
)


def render(path, params, output_dir):
    fields = load_fields(path)
    X, Y = fields["X"], fields["Y"]
    title = f"Ra = {params['Ra']}  |  Ha = {params['Ha']}  |  N = {params['N']}  |  Pr = {params['Pr']}"

    fig, (ax_t, ax_c) = plt.subplots(1, 2, figsize=(10, 4.2), dpi=120)
    fig.suptitle(title, y=0.98, weight="bold")

    isotherms = ax_t.contourf(X, Y, fields["T"], levels=25, cmap="inferno")
    fig.colorbar(isotherms, ax=ax_t, fraction=0.046, pad=0.04).set_label("Temperature (T)", fontsize=9)
    ax_t.set_title("Isotherms (T)")
    ax_t.set_xlabel("X (Lattice Unit)")
    ax_t.set_ylabel("Y (Lattice Unit)")
    ax_t.set_aspect("equal", "box")

    concentration = ax_c.contourf(X, Y, fields["C"], levels=25, cmap="viridis")
    fig.colorbar(concentration, ax=ax_c, fraction=0.046, pad=0.04).set_label("Concentration (C)", fontsize=9)

    speed = np.hypot(fields["u"], fields["v"])
    peak = speed.max() if speed.max() > 0 else 1.0
    ax_c.streamplot(
        X[0, :],
        Y[:, 0],
        fields["u"],
        fields["v"],
        color="white",
        linewidth=0.8 + 1.2 * speed / peak,
        density=1.0,
        arrowsize=0.7,
    )
    ax_c.set_title("Concentration & Streamlines")
    ax_c.set_xlabel("X (Lattice Unit)")
    ax_c.set_ylabel("Y (Lattice Unit)")
    ax_c.set_aspect("equal", "box")

    fig.tight_layout()
    fig.savefig(output_dir / f"{path.stem}_plots.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir")
    parser.add_argument("output_dir", nargs="?", default="output_plots")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = list_grid_files(args.data_dir)
    for index, (path, params) in enumerate(files, 1):
        render(path, params, output_dir)
        print(f"[{index}/{len(files)}] {path.name}")


if __name__ == "__main__":
    main()
