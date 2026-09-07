#!/usr/bin/env python3
"""Apply the extended expansion-threshold grid to final tables and figures.

This runs after analysis_15, which has already curated the final output tree and moved
root CSV caches into results/cache. It overwrites only the threshold-sensitivity table
and the two sensitivity panels, so no duplicate figures are introduced.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FINAL = RESULTS / "final"
CACHE = RESULTS / "cache"
FIG = ROOT / "figures"

THRESHOLD_PCT = [0.05, 0.10] + [round(x / 10, 2) for x in range(2, 16)]
source = CACHE / "expanded_threshold_sensitivity_extended.csv"
master_path = FINAL / "01_clone_master.csv"
if not source.exists():
    raise FileNotFoundError(source)
if not master_path.exists():
    raise FileNotFoundError(master_path)

sens = pd.read_csv(source)
master = pd.read_csv(master_path)
observed = {round(float(x), 8) for x in sens.threshold_percent.dropna().unique()}
if observed != set(THRESHOLD_PCT):
    raise RuntimeError(f"Unexpected extended threshold set: {sorted(observed)}")
sens.to_csv(FINAL / "04_threshold_sensitivity.csv", index=False)


def truthy(s):
    return s.astype(str).str.lower().isin(["true", "1", "yes"])

primary = pd.to_numeric(master["T1_fraction"], errors="coerce").ge(0.001)

# Overwrite rarefied 02: same sampling-null panel, extended threshold grid on the right.
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), constrained_layout=True)
z = master[primary].copy()
pdet = pd.to_numeric(z.get("P_detection_if_unchanged_empirical", np.nan), errors="coerce")
axes[0].scatter(100 * pd.to_numeric(z.T1_fraction, errors="coerce"), pdet, s=14, alpha=0.35)
if "sampling_robust_not_detected" in z.columns:
    robust = truthy(z["sampling_robust_not_detected"])
    axes[0].scatter(
        100 * pd.to_numeric(z.loc[robust, "T1_fraction"], errors="coerce"),
        pdet.loc[robust], marker="*", s=60, label="Sampling-robust not detected"
    )
axes[0].set_xlabel("T1 repertoire fraction, percent")
axes[0].set_ylabel("P detection if unchanged")
axes[0].set_ylim(0, 1.02)
axes[0].set_title("Internal T1 sampling null")
axes[0].grid(alpha=0.2)
if "sampling_robust_not_detected" in z.columns:
    axes[0].legend(fontsize=8)

rare = sens[sens.branch.eq("pairwise_rarefied")].copy()
for age, label in [("child", "Children"), ("adult", "Adults")]:
    q = rare[rare.age_group.eq(age)]
    g = q.groupby("threshold_percent")["not_detected_mass_percent_median"].mean().reindex(THRESHOLD_PCT)
    axes[1].plot(g.index, g.values, marker="o", label=label)
axes[1].set_xlabel("Expanded threshold, percent")
axes[1].set_ylabel("Mean T1 mass not detected at T2, percent")
axes[1].set_title("Threshold sensitivity from 0.05 to 1.50 percent")
axes[1].set_xticks(THRESHOLD_PCT)
axes[1].tick_params(axis="x", rotation=45)
axes[1].legend()
axes[1].grid(alpha=0.2)
fig.savefig(FIG / "rarefied" / "02_sampling_and_sensitivity.png", dpi=220, bbox_inches="tight")
plt.close(fig)

# Overwrite Simpson raw 02: shared-clone representation plus the same full threshold grid.
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), constrained_layout=True)
z = master[primary & pd.to_numeric(master.T2_fraction, errors="coerce").gt(0)].copy()
x = 100 * pd.to_numeric(z.T1_fraction, errors="coerce")
y = 100 * pd.to_numeric(z.T2_fraction, errors="coerce")
axes[0].scatter(x, y, s=14, alpha=0.35)
lim = max(float(np.nanmax(x)) if len(x) else 1.0, float(np.nanmax(y)) if len(y) else 1.0)
axes[0].plot([0, lim], [0, lim], linestyle="--", linewidth=1)
axes[0].set_xlabel("T1 repertoire fraction, percent")
axes[0].set_ylabel("T2 repertoire fraction, percent")
axes[0].set_title("Shared T1-expanded clonotypes without downsampling")
axes[0].grid(alpha=0.2)

raw = sens[sens.branch.eq("simpson_no_downsampling")].copy()
for age, label in [("child", "Children"), ("adult", "Adults")]:
    q = raw[raw.age_group.eq(age)]
    g = q.groupby("threshold_percent")["not_detected_mass_percent_median"].mean().reindex(THRESHOLD_PCT)
    axes[1].plot(g.index, g.values, marker="o", label=label)
axes[1].set_xlabel("Expanded threshold, percent")
axes[1].set_ylabel("Mean T1 mass not detected at T2, percent")
axes[1].set_title("Threshold sensitivity without downsampling")
axes[1].set_xticks(THRESHOLD_PCT)
axes[1].tick_params(axis="x", rotation=45)
axes[1].legend()
axes[1].grid(alpha=0.2)
fig.savefig(FIG / "simpson_raw" / "02_shared_and_sensitivity.png", dpi=220, bbox_inches="tight")
plt.close(fig)

# Keep the figure index authoritative and human-readable.
index = FIG / "FIGURE_INDEX.txt"
if index.exists():
    text = index.read_text(encoding="utf-8")
    lines = []
    for line in text.splitlines():
        if line.lower().startswith("sensitivity thresholds"):
            line = "Sensitivity thresholds: 0.05, 0.10, then 0.20 to 1.50 percent in 0.10 percentage-point steps"
        lines.append(line)
    index.write_text("\n".join(lines) + "\n", encoding="utf-8")

print("Extended threshold grid applied to final output")
print("Thresholds percent:", THRESHOLD_PCT)
