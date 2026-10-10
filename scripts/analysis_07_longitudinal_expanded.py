#!/usr/bin/env python3
"""Compare expanded aaV clonotypes in patients with accepted T1 and T2 samples.

Clonality and Shannon diversity use the complete repertoire. Overlap uses
clonotypes at or above 0.10 percent of productive TRB reads. Bulk TRB sequencing
does not identify naive T cells by phenotype."""

from pathlib import Path
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, mannwhitneyu
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FIG = ROOT / "figures" / "current"
FIG.mkdir(parents=True, exist_ok=True)
CURRENT.mkdir(parents=True, exist_ok=True)

EXPANDED_FRACTION = 0.001
metrics = pd.read_csv(RESULTS / "diversity_metrics.csv")
clones = pd.read_csv(RESULTS / "combined_trb_clonotypes.csv")
meta = pd.read_csv(CURRENT / "sample_metadata_current.csv")
clones["fraction"] = pd.to_numeric(clones["fraction"], errors="coerce").fillna(0.0)
clones["expanded"] = clones["fraction"] >= EXPANDED_FRACTION
for c in ["analysis_include_final", "is_control"]:
    if c in meta.columns:
        meta[c] = meta[c].astype(str).str.lower().isin(["true", "1", "yes"])

case_metrics = metrics[metrics["disease_status"].eq("measles")].copy()

# Exact specificity matches are optional and are used only for new T2 expansions
specificity = {}
match_path = RESULTS / "antigen_specific_clonotypes.csv"
if match_path.exists() and match_path.stat().st_size > 0:
    m = pd.read_csv(match_path)
    if {"sample_id", "cdr3aa", "pathogen"}.issubset(m.columns):
        m = m.dropna(subset=["sample_id", "cdr3aa", "pathogen"])
        for key, z in m.groupby(["sample_id", "cdr3aa"]):
            specificity[key] = "; ".join(sorted(set(z["pathogen"].astype(str))))


def point_table(patient, point):
    return clones[
        (clones.patient_id == patient) & (clones.clinical_timepoint == point)
    ].copy()


def clone_key_frame(z):
    x = z[["cdr3aa", "v_gene", "fraction", "expanded"]].copy()
    x["clone_key"] = list(zip(x.cdr3aa.astype(str), x.v_gene.astype(str)))
    return x.set_index("clone_key")


pair_rows = []
fate_rows = []
new_rows = []

for age in ["child", "adult"]:
    z = case_metrics[case_metrics.age_group.eq(age)]
    t1 = set(z.loc[z.clinical_timepoint.eq("T1"), "patient_id"])
    t2 = set(z.loc[z.clinical_timepoint.eq("T2"), "patient_id"])
    patients = sorted(t1 & t2)
    expected = 6 if age == "child" else 16
    if len(patients) != expected:
        raise RuntimeError(
            f"Unexpected paired T1 and T2 count for {age}: expected {expected}, found {len(patients)}"
        )

    for patient in patients:
        p1 = z[(z.patient_id == patient) & z.clinical_timepoint.eq("T1")].iloc[0]
        p2 = z[(z.patient_id == patient) & z.clinical_timepoint.eq("T2")].iloc[0]
        a = clone_key_frame(point_table(patient, "T1"))
        b = clone_key_frame(point_table(patient, "T2"))
        exp1 = set(a.index[a.expanded])
        exp2 = set(b.index[b.expanded])
        union = exp1 | exp2
        inter = exp1 & exp2
        jaccard = len(inter) / len(union) if union else np.nan

        # Batch information is kept for sensitivity analysis
        pm = meta[
            meta.patient_id.eq(patient)
            & meta.analysis_include_final
            & meta.clinical_timepoint.isin(["T1", "T2"])
        ]
        batches_t1 = sorted(
            set(pm.loc[pm.clinical_timepoint.eq("T1"), "batch"].astype(str))
        )
        batches_t2 = sorted(
            set(pm.loc[pm.clinical_timepoint.eq("T2"), "batch"].astype(str))
        )
        cross_batch = set(batches_t1) != set(batches_t2)

        retained = 0
        below = 0
        undetected = 0
        retained_mass = 0.0
        below_mass = 0.0
        undetected_mass = 0.0
        t1_expanded_mass = float(a.loc[list(exp1), "fraction"].sum()) if exp1 else 0.0

        for key in sorted(exp1):
            t1f = float(a.at[key, "fraction"])
            t2f = float(b.at[key, "fraction"]) if key in b.index else 0.0
            if key in exp2:
                fate = "retained expanded"
                retained += 1
                retained_mass += t1f
            elif key in b.index:
                fate = "detected below expansion threshold"
                below += 1
                below_mass += t1f
            else:
                fate = "not detected at T2"
                undetected += 1
                undetected_mass += t1f
            fate_rows.append(
                {
                    "patient_id": patient,
                    "age_group": age,
                    "cdr3aa": key[0],
                    "v_gene": key[1],
                    "T1_fraction": t1f,
                    "T2_fraction": t2f,
                    "fate": fate,
                }
            )

        for key in sorted(exp2 - exp1):
            t2f = float(b.at[key, "fraction"])
            t1f = float(a.at[key, "fraction"]) if key in a.index else 0.0
            origin = (
                "present below threshold at T1"
                if key in a.index
                else "not detected at T1"
            )
            sample_id = f"{patient}__T2"
            new_rows.append(
                {
                    "patient_id": patient,
                    "age_group": age,
                    "cdr3aa": key[0],
                    "v_gene": key[1],
                    "T1_fraction": t1f,
                    "T2_fraction": t2f,
                    "new_expansion_origin": origin,
                    "exact_reference_match": specificity.get((sample_id, key[0]), ""),
                }
            )

        pair_rows.append(
            {
                "patient_id": patient,
                "age_group": age,
                "cross_batch_pair": cross_batch,
                "n_T1_libraries": int(p1.n_replicate_libraries),
                "n_T2_libraries": int(p2.n_replicate_libraries),
                "T1_clonality": p1.clonality,
                "T2_clonality": p2.clonality,
                "clonality_change": p2.clonality - p1.clonality,
                "T1_shannon": p1.shannon,
                "T2_shannon": p2.shannon,
                "shannon_change": p2.shannon - p1.shannon,
                "T1_expanded_clonotypes": len(exp1),
                "T2_expanded_clonotypes": len(exp2),
                "expanded_jaccard": jaccard,
                "T1_expanded_read_share": t1_expanded_mass,
                "retained_expanded_count": retained,
                "detected_below_threshold_count": below,
                "not_detected_T2_count": undetected,
                "retained_T1_expanded_read_share": retained_mass,
                "below_threshold_T1_expanded_read_share": below_mass,
                "not_detected_T1_expanded_read_share": undetected_mass,
                "new_T2_expanded_count": len(exp2 - exp1),
            }
        )

pairs = pd.DataFrame(pair_rows)
fates = pd.DataFrame(fate_rows)
new_expansions = pd.DataFrame(new_rows)
pairs.to_csv(CURRENT / "longitudinal_expanded_clone_summary.csv", index=False)
fates.to_csv(CURRENT / "expanded_clone_fates.csv", index=False)
new_expansions.to_csv(CURRENT / "new_T2_expansions.csv", index=False)

stats = []
for age in ["child", "adult"]:
    p = pairs[pairs.age_group.eq(age)]
    controls = metrics[
        (metrics.age_group == age) & metrics.clinical_timepoint.eq("CONTROL")
    ]
    for metric, a_col, b_col in [
        ("clonality", "T1_clonality", "T2_clonality"),
        ("shannon", "T1_shannon", "T2_shannon"),
    ]:
        x = pd.to_numeric(p[a_col], errors="coerce").to_numpy(float)
        y = pd.to_numeric(p[b_col], errors="coerce").to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
        pw = (
            wilcoxon(y, x, alternative="two-sided").pvalue
            if len(x) >= 3 and np.any(np.abs(y - x) > 1e-12)
            else np.nan
        )
        c = pd.to_numeric(controls[metric], errors="coerce").dropna().to_numpy(float)
        p_c1 = (
            mannwhitneyu(c, x, alternative="two-sided").pvalue
            if len(c) >= 2 and len(x) >= 2
            else np.nan
        )
        p_c2 = (
            mannwhitneyu(c, y, alternative="two-sided").pvalue
            if len(c) >= 2 and len(y) >= 2
            else np.nan
        )
        stats.extend(
            [
                {
                    "age_group": age,
                    "metric": metric,
                    "comparison": "paired T1 compared with T2",
                    "n_a": len(x),
                    "n_b": len(y),
                    "p_value": pw,
                },
                {
                    "age_group": age,
                    "metric": metric,
                    "comparison": "controls compared with paired T1",
                    "n_a": len(c),
                    "n_b": len(x),
                    "p_value": p_c1,
                },
                {
                    "age_group": age,
                    "metric": metric,
                    "comparison": "controls compared with paired T2",
                    "n_a": len(c),
                    "n_b": len(y),
                    "p_value": p_c2,
                },
            ]
        )
stat = pd.DataFrame(stats)
stat["fdr"] = np.nan
for age in stat.age_group.unique():
    idx = stat.index[(stat.age_group == age) & stat.p_value.notna()]
    if len(idx):
        stat.loc[idx, "fdr"] = multipletests(stat.loc[idx, "p_value"], method="fdr_bh")[
            1
        ]
stat.to_csv(CURRENT / "longitudinal_statistics.csv", index=False)


def paired_metric_plot(age, metric, ylab, filename):
    p = pairs[pairs.age_group.eq(age)].copy()
    controls = metrics[
        (metrics.age_group == age) & metrics.clinical_timepoint.eq("CONTROL")
    ]
    c = pd.to_numeric(controls[metric], errors="coerce").dropna().to_numpy(float)
    a = p[f"T1_{metric}"].to_numpy(float)
    b = p[f"T2_{metric}"].to_numpy(float)
    fig, ax = plt.subplots(figsize=(8.8, 6.4))
    if len(c):
        ax.scatter(
            np.full(len(c), 1.0) + np.linspace(-0.08, 0.08, len(c)),
            c,
            s=30,
            alpha=0.75,
            label="Healthy controls",
        )
    for i in range(len(p)):
        ax.plot([2, 3], [a[i], b[i]], linewidth=1.1, alpha=0.55)
        ax.scatter([2, 3], [a[i], b[i]], s=30)
    ax.set_xticks(
        [1, 2, 3],
        [
            f"Healthy controls\nn {len(c)}",
            f"Paired T1\nn {len(p)}",
            f"Paired T2\nn {len(p)}",
        ],
    )
    ax.set_ylabel(ylab)
    age_title = (
        "Children younger than 18 years"
        if age == "child"
        else "Adults 18 years and older"
    )
    ax.set_title(f"{age_title}: paired change in {ylab.lower()}")
    ax.grid(axis="y", alpha=0.2)
    fig.text(
        0.5,
        0.015,
        "Lines connect the same patient from T1 to T2. Controls are shown as an external reference",
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(FIG / f"{filename}.png", dpi=210, bbox_inches="tight")
    fig.savefig(FIG / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)


def fate_plot(age, filename):
    p = pairs[pairs.age_group.eq(age)].copy()
    denom = (
        p.retained_expanded_count
        + p.detected_below_threshold_count
        + p.not_detected_T2_count
    ).replace(0, np.nan)
    vals = pd.DataFrame(
        {
            "Retained expanded": 100 * p.retained_expanded_count / denom,
            "Detected below threshold": 100 * p.detected_below_threshold_count / denom,
            "Not detected at T2": 100 * p.not_detected_T2_count / denom,
        }
    )
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    x = np.arange(len(vals))
    bottom = np.zeros(len(vals))
    for col in vals.columns:
        y = vals[col].fillna(0).to_numpy(float)
        ax.bar(x, y, bottom=bottom, label=col)
        bottom += y
    ax.set_xticks(x, p.patient_id.tolist(), rotation=55, ha="right")
    ax.set_ylabel("Percent of expanded T1 clonotypes")
    age_title = (
        "Children younger than 18 years"
        if age == "child"
        else "Adults 18 years and older"
    )
    ax.set_title(f"{age_title}: fate of expanded T1 clonotypes after recovery")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 100)
    fig.text(
        0.5,
        0.015,
        "Overlap is evaluated only among clonotypes expanded to at least 0.10 percent of productive TRB reads",
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(FIG / f"{filename}.png", dpi=210, bbox_inches="tight")
    fig.savefig(FIG / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)


for age, suffix in [("child", "children"), ("adult", "adults")]:
    paired_metric_plot(age, "clonality", "Clonality", f"04_paired_clonality_{suffix}")
    paired_metric_plot(
        age, "shannon", "Shannon diversity", f"05_paired_shannon_{suffix}"
    )
    fate_plot(age, f"06_expanded_clone_fate_{suffix}")

print("Paired expanded clone analysis completed")
print(pairs.groupby("age_group").patient_id.nunique())
if "БА1_071" in set(pairs.patient_id):
    r = pairs[pairs.patient_id.eq("БА1_071")].iloc[0]
    print(f"БА1_071 T2 uses {int(r.n_T2_libraries)} averaged libraries")
