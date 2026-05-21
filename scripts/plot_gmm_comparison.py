#!/usr/bin/env python3
"""Plot GMM comparison: MC scatter + GMM ellipses for 6x6 state subspace."""

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import chi2
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJ_ROOT, "OPM", "L1Halo_0noise", "gmm_analysis")
OUT_DIR = DATA_DIR

LABELS = ["X (km)", "Y (km)", "Z (km)", "VX (km/s)", "VY (km/s)", "VZ (km/s)"]
SHORT_LABELS = ["X", "Y", "Z", "VX", "VY", "VZ"]
COMP_LIST = [1, 3, 5]
METHODS = ["MC", "DA5"]
METHOD_COLORS = {"MC": "#1f77b4", "DA5": "#d62728"}
SCATTER_ALPHA = 0.15
SCATTER_COLOR = "#888888"


def load_csv(path):
    return np.loadtxt(path, delimiter=",", skiprows=1).T  # (6, N)


def load_json(path):
    with open(path) as f:
        return json.load(f)


def cov_ellipse(cov_2d, mean, n_std, **kwargs):
    """Draw confidence ellipse from 2x2 covariance."""
    if np.any(~np.isfinite(cov_2d)) or np.max(np.abs(cov_2d)) < 1e-30:
        return
    vals, vecs = np.linalg.eigh(cov_2d)
    if np.any(vals <= 0):
        return
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    scale = chi2.ppf(1 - 2 * (1 - 0.6827), 2) if n_std == 1 else chi2.ppf(0.95, 2)
    width, height = 2 * np.sqrt(vals * scale)
    e = Ellipse(xy=mean, width=width, height=height, angle=angle,
                fill=False, **kwargs)
    ax = plt.gca()
    ax.add_patch(e)


def plot_one_figure(n_comp, mc_particles, mc_json, da5_json):
    """Draw 6x6 subplot grid for one n_comp."""
    fig, axes = plt.subplots(6, 6, figsize=(24, 22))
    fig.suptitle(f"GMM Comparison: {n_comp} Component(s) — MC scatter (grey), "
                 f"MC-GMM (blue), DA5-GMM (red)",
                 fontsize=14, fontweight="bold")

    for row in range(6):
        for col in range(6):
            ax = axes[row, col]
            # Scatter MC particles
            ax.scatter(mc_particles[col, :], mc_particles[row, :],
                       s=1, c=SCATTER_COLOR, alpha=SCATTER_ALPHA, rasterized=True)

            # MC-GMM ellipses
            for gmm_data in [mc_json, da5_json]:
                color = METHOD_COLORS[gmm_data["method"]]
                for comp in gmm_data["components"]:
                    w = comp["weight"]
                    if w < 0.01:
                        continue
                    mean_2d = np.array([comp["mean"][col], comp["mean"][row]])
                    cov_full = np.array(comp["cov"])  # 6x6 2D array
                    cov_2d = cov_full[np.ix_([row, col], [row, col])]
                    alpha_lw = max(0.3, min(1.0, w * 2))
                    cov_ellipse(cov_2d, mean_2d, 1, ec=color,
                                linewidth=1.5 * alpha_lw, alpha=alpha_lw,
                                linestyle="-")
                    cov_ellipse(cov_2d, mean_2d, 2, ec=color,
                                linewidth=0.8 * alpha_lw, alpha=alpha_lw * 0.5,
                                linestyle="--")

            # Also plot GMM global means
            for gmm_data, color in [(mc_json, METHOD_COLORS["MC"]),
                                     (da5_json, METHOD_COLORS["DA5"])]:
                gm = gmm_data["global_mean"]
                ax.plot(gm[col], gm[row], marker='x', color=color,
                        markersize=8, mew=2)

            if row == 5:
                ax.set_xlabel(LABELS[col], fontsize=8)
            else:
                ax.set_xticklabels([])
            if col == 0:
                ax.set_ylabel(LABELS[row], fontsize=8)
            else:
                ax.set_yticklabels([])

            ax.tick_params(labelsize=7)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=SCATTER_COLOR,
               markersize=6, label='MC particles'),
        Line2D([0], [0], color=METHOD_COLORS["MC"], linewidth=1.5, label='MC-GMM 1σ'),
        Line2D([0], [0], color=METHOD_COLORS["MC"], linewidth=0.8,
               linestyle='--', label='MC-GMM 2σ'),
        Line2D([0], [0], color=METHOD_COLORS["DA5"], linewidth=1.5, label='DA5-GMM 1σ'),
        Line2D([0], [0], color=METHOD_COLORS["DA5"], linewidth=0.8,
               linestyle='--', label='DA5-GMM 2σ'),
        Line2D([0], [0], marker='x', color=METHOD_COLORS["MC"],
               markersize=8, mew=2, label='MC global mean', linestyle='None'),
        Line2D([0], [0], marker='x', color=METHOD_COLORS["DA5"],
               markersize=8, mew=2, label='DA5 global mean', linestyle='None'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4,
               fontsize=9, frameon=True)

    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    out_path = os.path.join(OUT_DIR, f"gmm_comp{n_comp}_comparison.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_component_summary(mc_particles, all_jsons):
    """Single summary figure: 3 rows (n_comp) x 1 col showing key 2D projection."""
    # Pick the most informative 2D projection: X-Y position
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("GMM Fit Comparison: X-Y Projection (MC scatter + GMM ellipses)",
                 fontsize=14, fontweight="bold")

    for idx, n_comp in enumerate(COMP_LIST):
        ax = axes[idx]
        ax.set_title(f"{n_comp} Component(s)", fontsize=12)
        ax.set_xlabel("X (km)")
        ax.set_ylabel("Y (km)")

        x_idx, y_idx = 0, 1
        ax.scatter(mc_particles[x_idx, :], mc_particles[y_idx, :],
                   s=1, c=SCATTER_COLOR, alpha=SCATTER_ALPHA, rasterized=True)

        mc_json = all_jsons[f"MC_{n_comp}"]
        da5_json = all_jsons[f"DA5_{n_comp}"]

        for gmm_data in [mc_json, da5_json]:
            color = METHOD_COLORS[gmm_data["method"]]
            for comp in gmm_data["components"]:
                w = comp["weight"]
                if w < 0.01:
                    continue
                mean_2d = np.array([comp["mean"][x_idx], comp["mean"][y_idx]])
                cov_full = np.array(comp["cov"])  # 6x6 2D array
                cov_2d = cov_full[np.ix_([y_idx, x_idx], [y_idx, x_idx])]
                alpha_lw = max(0.3, min(1.0, w * 2))
                cov_ellipse(cov_2d, mean_2d, 1, ec=color,
                            linewidth=1.5 * alpha_lw, alpha=alpha_lw, linestyle="-")
                cov_ellipse(cov_2d, mean_2d, 2, ec=color,
                            linewidth=0.8 * alpha_lw, alpha=alpha_lw * 0.5,
                            linestyle="--")

            gm = gmm_data["global_mean"]
            ax.plot(gm[x_idx], gm[y_idx], marker='x', color=color,
                    markersize=10, mew=2.5)

        ax.set_aspect('auto')

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=SCATTER_COLOR,
               markersize=6, label='MC particles'),
        Line2D([0], [0], color=METHOD_COLORS["MC"], linewidth=1.5, label='MC-GMM 1σ'),
        Line2D([0], [0], color=METHOD_COLORS["MC"], linewidth=0.8,
               linestyle='--', label='MC-GMM 2σ'),
        Line2D([0], [0], color=METHOD_COLORS["DA5"], linewidth=1.5, label='DA5-GMM 1σ'),
        Line2D([0], [0], color=METHOD_COLORS["DA5"], linewidth=0.8,
               linestyle='--', label='DA5-GMM 2σ'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=5,
               fontsize=9, frameon=True)
    plt.tight_layout(rect=[0, 0.06, 1, 0.94])
    out_path = os.path.join(OUT_DIR, "gmm_summary_xy.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def main():
    print("Loading MC particles...")
    mc_particles = load_csv(os.path.join(DATA_DIR, "mc_particles.csv"))
    print(f"  MC particles shape: {mc_particles.shape}")

    # Load all JSONs
    all_jsons = {}
    for n_comp in COMP_LIST:
        for method in METHODS:
            key = f"{method}_{n_comp}"
            fname = f"{method.lower()}_gmm_comp{n_comp}.json"
            all_jsons[key] = load_json(os.path.join(DATA_DIR, fname))
            g = all_jsons[key]
            print(f"  Loaded {fname}: method={g['method']}, "
                  f"n_comp={g['n_components']}, "
                  f"components: {[f'{c[\"weight\"]:.3f}' for c in g['components']]}")

    # Plot 6x6 figures for each n_comp
    for n_comp in COMP_LIST:
        print(f"\nPlotting 6x6 figure for {n_comp} component(s)...")
        mc_json = all_jsons[f"MC_{n_comp}"]
        da5_json = all_jsons[f"DA5_{n_comp}"]
        plot_one_figure(n_comp, mc_particles, mc_json, da5_json)

    # Plot summary figure
    print("\nPlotting summary figure...")
    plot_component_summary(mc_particles, all_jsons)

    print("\nDone! All figures saved to", OUT_DIR)


if __name__ == "__main__":
    main()
