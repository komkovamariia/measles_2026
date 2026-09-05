#!/usr/bin/env python3
"""Cross sectional repertoire analysis

Idea one uses every accepted T1 observation, every accepted T2 observation and
all accepted healthy controls. T1 observations are retained even when a patient
has no T2 sample. Independent tests compare healthy controls with T1 and healthy
controls with T2. T1 and T2 are not treated as independent groups for a direct
between group test because some observations belong to the same patient.
"""
from pathlib import Path
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FIG = ROOT / "figures" / "current"
FIG.mkdir(parents=True, exist_ok=True)
CURRENT.mkdir(parents=True, exist_ok=True)

metrics = pd.read_csv(RESULTS / "diversity_metrics.csv")
meta = pd.read_csv(CURRENT / "sample_metadata_current.csv")

for c in ["analysis_include_final", "is_control"]:
    if c in meta.columns:
        meta[c] = meta[c].astype(str).str.lower().isin(["true", "1", "yes"])

if meta["clinical_timepoint"].astype(str).eq("T0").any():
    raise RuntimeError("T0 does not exist in this study")

# Cohort counts use patient timepoint observations, not sequencing libraries
count_rows = []
for age in ["child", "adult"]:
    z = metrics[metrics["age_group"].eq(age)]
    for point in ["CONTROL", "T1", "T2", "T3"]:
        p = z[z["clinical_timepoint"].eq(point)]
        count_rows.append({
            "age_group": age,
            "clinical_group": point,
            "n_patient_observations": int(p["patient_id"].nunique()),
        })
counts = pd.DataFrame(count_rows)
counts.to_csv(CURRENT / "cross_sectional_cohort_counts.csv", index=False)

expected = {
    ("child", "CONTROL"): 4,
    ("child", "T1"): 27,
    ("child", "T2"): 8,
    ("adult", "CONTROL"): 13,
    ("adult", "T1"): 29,
    ("adult", "T2"): 19,
}
for key, n in expected.items():
    got = int(counts.loc[(counts.age_group == key[0]) & (counts.clinical_group == key[1]), "n_patient_observations"].iloc[0])
    if got != n:
        raise RuntimeError(f"Unexpected cohort size for {key}: expected {n}, found {got}")

metric_specs = [
    ("clonality", "Clonality"),
    ("shannon", "Shannon diversity"),
    ("expanded_read_fraction", "Expanded clone read share"),
    ("n_expanded_clonotypes", "Number of expanded clonotypes"),
]

stats = []
for age in ["child", "adult"]:
    z = metrics[metrics["age_group"].eq(age)].copy()
    for metric, label in metric_specs:
        for point in ["T1", "T2"]:
            x = pd.to_numeric(z.loc[z.clinical_timepoint.eq("CONTROL"), metric], errors="coerce").dropna().to_numpy(float)
            y = pd.to_numeric(z.loc[z.clinical_timepoint.eq(point), metric], errors="coerce").dropna().to_numpy(float)
            p = mannwhitneyu(x, y, alternative="two-sided").pvalue if len(x) >= 2 and len(y) >= 2 else np.nan
            stats.append({
                "age_group": age,
                "comparison": f"CONTROL compared with {point}",
                "metric": metric,
                "n_control": len(x),
                "n_case": len(y),
                "median_control": np.median(x) if len(x) else np.nan,
                "median_case": np.median(y) if len(y) else np.nan,
                "median_difference_case_minus_control": np.median(y) - np.median(x) if len(x) and len(y) else np.nan,
                "p_value": p,
            })
stat = pd.DataFrame(stats)
stat["fdr"] = np.nan
for age in stat.age_group.unique():
    idx = stat.index[(stat.age_group == age) & stat.p_value.notna()]
    if len(idx):
        stat.loc[idx, "fdr"] = multipletests(stat.loc[idx, "p_value"], method="fdr_bh")[1]
stat.to_csv(CURRENT / "cross_sectional_statistics.csv", index=False)


def jitter(n, center, width=0.10):
    if n <= 1:
        return np.full(n, center, dtype=float)
    return np.linspace(center - width, center + width, n)


def plot_metric(age, metric, label, filename):
    z = metrics[metrics["age_group"].eq(age)].copy()
    groups = ["CONTROL", "T1", "T2"]
    display = ["Healthy controls", "T1 early infection", "T2 recovery"]
    arrays = [pd.to_numeric(z.loc[z.clinical_timepoint.eq(g), metric], errors="coerce").dropna().to_numpy(float) for g in groups]
    fig, ax = plt.subplots(figsize=(8.4, 6.2))
    ax.boxplot(arrays, positions=np.arange(1, 4), widths=0.52, showfliers=False)
    for pos, arr in enumerate(arrays, start=1):
        ax.scatter(jitter(len(arr), pos), arr, s=24, alpha=0.72)
    ax.set_xticks(np.arange(1, 4), [f"{name}\nn {len(arr)}" for name, arr in zip(display, arrays)])
    ax.set_ylabel(label)
    age_title = "Children younger than 18 years" if age == "child" else "Adults 18 years and older"
    ax.set_title(f"{age_title}: {label.lower()} across clinical groups")
    ax.grid(axis="y", alpha=0.2)
    fig.text(
        0.5, 0.015,
        "All accepted T1 observations are included even when T2 is unavailable. Independent tests compare controls with T1 and controls with T2",
        ha="center", fontsize=8.5,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(FIG / f"{filename}.png", dpi=210, bbox_inches="tight")
    fig.savefig(FIG / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)

for age, suffix in [("child", "children"), ("adult", "adults")]:
    plot_metric(age, "clonality", "Clonality", f"01_cross_sectional_clonality_{suffix}")
    plot_metric(age, "shannon", "Shannon diversity", f"02_cross_sectional_shannon_{suffix}")
    plot_metric(age, "expanded_read_fraction", "Expanded clone read share", f"03_cross_sectional_expanded_share_{suffix}")

print("Cross sectional analysis completed")
print(counts.to_string(index=False))
