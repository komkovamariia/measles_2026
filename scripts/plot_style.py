"""Shared publication-quality plotting style for the measles TRB analysis."""
from __future__ import annotations

import numpy as np
import matplotlib as mpl


COLORS = {
    "navy": "#24435A",
    "blue": "#4C78A8",
    "teal": "#2A9D8F",
    "orange": "#E07A5F",
    "gold": "#E9C46A",
    "purple": "#8064A2",
    "gray": "#7A858F",
    "light_gray": "#D9E0E6",
    "ink": "#22313A",
    "paper": "#FFFFFF",
}

FATE_COLORS = {
    "Retained expanded": COLORS["teal"],
    "Detected below threshold": COLORS["gold"],
    "Not detected": COLORS["orange"],
    "Not detected later": COLORS["orange"],
}

TIME_COLORS = {"T1": COLORS["blue"], "T2": COLORS["orange"]}
AGE_COLORS = {"child": COLORS["teal"], "adult": COLORS["purple"]}
BRANCH_COLORS = {
    "pairwise_rarefied": COLORS["teal"],
    "simpson_no_downsampling": COLORS["blue"],
}


def configure_publication_style() -> None:
    mpl.rcParams.update({
        "figure.facecolor": COLORS["paper"],
        "axes.facecolor": COLORS["paper"],
        "savefig.facecolor": COLORS["paper"],
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.labelsize": 9.5,
        "axes.labelcolor": COLORS["ink"],
        "axes.edgecolor": COLORS["light_gray"],
        "axes.linewidth": 0.8,
        "xtick.color": COLORS["ink"],
        "ytick.color": COLORS["ink"],
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "grid.color": COLORS["light_gray"],
        "grid.linewidth": 0.7,
        "grid.alpha": 0.65,
        "lines.linewidth": 1.8,
        "lines.markersize": 5,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    })


def finish_axis(ax, grid: str | None = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(True, axis=grid)
        ax.set_axisbelow(True)


def panel_label(ax, label: str) -> None:
    ax.text(-0.10, 1.04, label, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=COLORS["ink"], va="bottom")


def annotate_ranked_points(ax, x, y, labels, indices, color=None) -> None:
    """Label up to ten points in two fixed columns with non-overlapping leaders."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    labels = np.asarray(labels, dtype=object)
    indices = list(indices)[:10]
    if not indices:
        return
    left = indices[::2]
    right = indices[1::2]
    for side, chosen in [("left", left), ("right", right)]:
        if not chosen:
            continue
        y_positions = np.linspace(0.91, 0.17, len(chosen))
        x_text = 0.02 if side == "left" else 0.98
        ha = "left" if side == "left" else "right"
        for idx, y_text in zip(chosen, y_positions):
            ax.annotate(
                str(labels[idx]), xy=(x[idx], y[idx]), xycoords="data",
                xytext=(x_text, y_text), textcoords="axes fraction",
                ha=ha, va="center", fontsize=7.1, color=COLORS["ink"],
                arrowprops={"arrowstyle": "-", "color": color or COLORS["gray"],
                            "lw": 0.7, "alpha": 0.85},
                bbox={"boxstyle": "round,pad=0.18", "facecolor": "white",
                      "edgecolor": COLORS["light_gray"], "alpha": 0.92},
                zorder=5,
            )


configure_publication_style()
