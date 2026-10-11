#!/usr/bin/env python3
"""Export final clone tables and figures for rarefied and observed-frequency analyses."""

from __future__ import annotations

from pathlib import Path
import math
import os
import shutil
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from plot_style import (
    AGE_COLORS,
    BRANCH_COLORS,
    COLORS,
    FATE_COLORS,
    TIME_COLORS,
    annotate_ranked_points,
    finish_axis,
    panel_label,
)

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FINAL = RESULTS / "final"
FIG = ROOT / "figures"

PRIMARY = 0.001
THRESHOLD_PCT = [0.05, 0.10] + [round(x / 10, 2) for x in range(2, 16)]

for p in [
    CURRENT / "paired_T1_clone_master_nt.csv",
    CURRENT / "pairwise_threshold_sensitivity.csv",
    CURRENT / "simpson_raw_patient_metrics.csv",
    CURRENT / "simpson_raw_threshold_sensitivity.csv",
    CURRENT / "simpson_raw_primary_candidates.csv",
    CURRENT / "downsampling_iteration_convergence.csv",
]:
    if not p.exists():
        raise FileNotFoundError(f"Required input is missing: {p}")

master = pd.read_csv(CURRENT / "paired_T1_clone_master_nt.csv")
rare_sens = pd.read_csv(CURRENT / "pairwise_threshold_sensitivity.csv")
rare_patient = pd.read_csv(CURRENT / "pairwise_rarefaction_summary.csv")
raw_metrics = pd.read_csv(CURRENT / "simpson_raw_patient_metrics.csv")
raw_sens = pd.read_csv(CURRENT / "simpson_raw_threshold_sensitivity.csv")
raw_candidates = pd.read_csv(CURRENT / "simpson_raw_primary_candidates.csv")
convergence = pd.read_csv(CURRENT / "downsampling_iteration_convergence.csv")

for c in [
    "T1_fraction",
    "T2_fraction",
    "T1_fraction_rarefied_median",
    "T2_fraction_rarefied_median",
    "P_detection_if_unchanged_empirical",
]:
    if c in master.columns:
        master[c] = pd.to_numeric(master[c], errors="coerce")


def truthy(s):
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def short_seq(s, n=18):
    s = str(s)
    if len(s) <= n:
        return s
    return s[:8] + "..." + s[-6:]


def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


known_col = find_col(
    master, ["known_antigen", "antigen", "matched_antigen", "specificity", "pathogen"]
)
if known_col is None:
    master["known_antigen_final"] = ""
else:
    master["known_antigen_final"] = master[known_col].fillna("").astype(str)


def specificity_class(x):
    antigen = str(x).strip()
    antigen_lower = antigen.lower()
    if not antigen or antigen_lower in {
        "nan",
        "none",
        "unknown",
        "no match",
        "not available",
    }:
        return "unknown"
    if "measles" in antigen_lower or "rubeola" in antigen_lower:
        return "measles-associated"
    return "known non-measles memory or other specificity"


master["specificity_class"] = master["known_antigen_final"].map(specificity_class)
master["literature_query"] = master.apply(
    lambda r: f'"{r.get("cdr3aa", "")}" "{r.get("TRBV", "")}" TCR viral infection'
    + (
        f' "{r.get("known_antigen_final", "")}"'
        if str(r.get("known_antigen_final", "")).strip()
        else ""
    ),
    axis=1,
)
master["literature_status"] = np.where(
    master.specificity_class.eq("unknown"),
    "requires PubMed review after candidate selection",
    "exact database annotation available; confirm in literature",
)

primary_mask = (
    truthy(master["expanded_T1_0_10_percent"])
    if "expanded_T1_0_10_percent" in master.columns
    else master.T1_fraction.ge(PRIMARY)
)
master["rarefied_log2_T2_over_T1"] = np.where(
    master.T1_fraction_rarefied_median.gt(0) & master.T2_fraction_rarefied_median.gt(0),
    np.log2(master.T2_fraction_rarefied_median / master.T1_fraction_rarefied_median),
    np.nan,
)


def rare_direction(r):
    if not bool(primary_mask.loc[r.name]):
        return "not primary T1-expanded"
    t2_fraction = r.get("T2_fraction_rarefied_median", np.nan)
    log2_change = r.get("rarefied_log2_T2_over_T1", np.nan)
    if not np.isfinite(t2_fraction) or t2_fraction <= 0:
        return "not detected"
    if log2_change >= 1:
        return "increased"
    if log2_change <= -1:
        return "decreased"
    return "stable"


master["rarefied_direction_0_10_percent"] = master.apply(rare_direction, axis=1)

# Merge raw branch fields back into the one-row-per-T1-clonotype master.
raw_keep_cols = [
    c
    for c in [
        "patient",
        "clone_id_nt",
        "raw_log2_T2_over_T1",
        "raw_absolute_delta",
        "raw_fold_contraction",
        "raw_T1_expanded_mass_share",
        "raw_direction_0_10_percent",
        "raw_strong_contraction",
    ]
    if c in raw_candidates.columns
]
master = master.merge(
    raw_candidates[raw_keep_cols], on=["patient", "clone_id_nt"], how="left"
)


def candidate_class(r, branch="rarefied"):
    if branch == "rarefied":
        direction = str(r.get("rarefied_direction_0_10_percent", ""))
        robust = str(r.get("sampling_robust_not_detected", "")).lower() in {
            "true",
            "1",
            "yes",
        }
        strong = str(r.get("sampling_candidate_class", "")) == "strong contraction"
    else:
        direction = str(r.get("raw_direction_0_10_percent", ""))
        robust = False
        strong = bool(r.get("raw_strong_contraction", False))
    specificity = str(r.get("specificity_class", "unknown"))
    if direction == "increased":
        return "candidate developing or persisting measles response"
    if direction == "stable":
        return "stable circulating clone"
    if direction in {"decreased", "not detected"} or robust or strong:
        if specificity == "measles-associated":
            return "possible contraction of measles response"
        if specificity.startswith("known non-measles"):
            return "candidate immune amnesia"
        return "contracted clone with unknown specificity"
    return "other"


master["candidate_class_rarefied"] = master.apply(
    candidate_class, axis=1, branch="rarefied"
)
master["candidate_class_simpson_raw"] = master.apply(
    candidate_class, axis=1, branch="raw"
)


def enrichment_for_branch(df, branch):
    # Reindex the fixed T1-expanded pool for each subgroup so enrichment
    # cannot include clones outside the primary pool.
    subgroup_primary = primary_mask.reindex(df.index).fillna(False).astype(bool)
    clones = df.loc[subgroup_primary].copy()
    if branch == "rarefied":
        affected = (
            clones.get("sampling_candidate_class", "")
            .astype(str)
            .isin(["sampling-robust not detected", "strong contraction"])
        )
        retained = (
            clones.get("sampling_candidate_class", "")
            .astype(str)
            .eq("retained expanded")
        )
    else:
        affected = clones.get("raw_direction_0_10_percent", "").astype(str).eq(
            "not detected"
        ) | clones.get("raw_strong_contraction", False).fillna(False).astype(bool)
        retained = clones.T2_fraction.ge(PRIMARY) & ~affected
    clones = clones[affected | retained].copy()
    if clones.empty or affected.loc[clones.index].nunique() < 2:
        return pd.DataFrame()
    clones["affected"] = affected.loc[clones.index].astype(int)
    clones["log10_T1_fraction"] = np.log10(
        pd.to_numeric(clones.T1_fraction, errors="coerce").clip(lower=1e-8)
    )

    rows = []

    def fit(kind, feature, indicator):
        if (
            int(indicator.sum()) < 5
            or clones.loc[indicator, "patient"].nunique() < 3
            or indicator.nunique() < 2
        ):
            return
        X = pd.DataFrame(
            {
                "const": 1.0,
                "feature": indicator.astype(float).to_numpy(),
                "log10_T1_fraction": clones.log10_T1_fraction.to_numpy(float),
            }
        )
        try:
            m = sm.GEE(
                clones.affected.to_numpy(float),
                X,
                groups=clones.patient.astype(str).to_numpy(),
                family=sm.families.Binomial(),
            ).fit(maxiter=100)
            coef = float(m.params["feature"])
            se = float(m.bse["feature"])
            p = float(m.pvalues["feature"])
        except Exception:
            coef = se = p = np.nan
        rows.append(
            {
                "branch": branch,
                "feature_type": kind,
                "feature": feature,
                "n_candidate_clones": int(indicator.sum()),
                "n_patients": int(clones.loc[indicator, "patient"].nunique()),
                "log_odds_affected_vs_retained": coef,
                "standard_error": se,
                "odds_ratio_affected_vs_retained": (
                    math.exp(coef) if np.isfinite(coef) else np.nan
                ),
                "p_value": p,
            }
        )

    for col, kind in [("TRBV", "TRBV"), ("TRBJ", "TRBJ")]:
        if col not in clones.columns:
            continue
        for value in sorted(v for v in clones[col].dropna().astype(str).unique() if v):
            fit(kind, value, clones[col].astype(str).eq(value))
    motifs = {}
    for idx, aa in clones.cdr3aa.astype(str).items():
        for kmer in set(aa[i : i + 3] for i in range(max(0, len(aa) - 2))):
            motifs.setdefault(kmer, set()).add(idx)
    for motif, indices in motifs.items():
        fit(
            "CDR3_3mer",
            motif,
            pd.Series(clones.index.isin(indices), index=clones.index),
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr"] = np.nan
    for kind in out.feature_type.unique():
        idx = out.index[(out.feature_type.eq(kind)) & out.p_value.notna()]
        if len(idx):
            out.loc[idx, "fdr"] = multipletests(
                out.loc[idx, "p_value"], method="fdr_bh"
            )[1]
    return out.sort_values(["fdr", "p_value"], na_position="last")


rare_enrich_path = CURRENT / "lost_vs_retained_feature_enrichment.csv"
rare_enrich = (
    pd.read_csv(rare_enrich_path)
    if rare_enrich_path.exists() and rare_enrich_path.stat().st_size > 1
    else pd.DataFrame()
)
if not rare_enrich.empty:
    rare_enrich = rare_enrich.rename(
        columns={"log_odds_lost_vs_retained": "log_odds_affected_vs_retained"}
    )
    rare_enrich["branch"] = "pairwise_rarefied"
raw_enrich = enrichment_for_branch(master, "simpson_no_downsampling")
enrichment = pd.concat([rare_enrich, raw_enrich], ignore_index=True, sort=False)

SUPP = RESULTS / "supplementary_tables"
for outdir in [FINAL, SUPP]:
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

paired_patients = set(rare_patient["patient"].astype(str))
master = master[master["patient"].astype(str).isin(paired_patients)].copy()
primary_mask = primary_mask.loc[master.index]
candidates = master[primary_mask].copy()
for source, target in [
    ("T1_fraction", "T1_fraction_percent"),
    ("T2_fraction", "T2_fraction_percent"),
    ("T1_fraction_rarefied_median", "T1_rarefied_median_percent"),
    ("T2_fraction_rarefied_median", "T2_rarefied_median_percent"),
]:
    candidates[target] = 100 * pd.to_numeric(candidates[source], errors="coerce")
candidates = candidates.rename(
    columns={
        "rarefied_direction_0_10_percent": "direction_pairwise_rarefaction",
        "raw_direction_0_10_percent": "direction_observed_frequency",
        "candidate_class_rarefied": "candidate_interpretation",
        "T1_raw_count_sum": "T1_count",
        "T2_raw_count_sum": "T2_count",
    }
)
t1_fraction = pd.to_numeric(candidates["T1_fraction"], errors="coerce")
t2_fraction = pd.to_numeric(candidates["T2_fraction"], errors="coerce")
observed_fold = (t1_fraction / t2_fraction).where(t2_fraction.gt(0))
if "fold_contraction_lower_bound_if_T2_zero" in candidates:
    lower_bound_fold = pd.to_numeric(
        candidates["fold_contraction_lower_bound_if_T2_zero"], errors="coerce"
    )
else:
    zero_upper = (
        pd.to_numeric(candidates["T2_zero_upper_95_pair_depth"], errors="coerce")
        if "T2_zero_upper_95_pair_depth" in candidates
        else pd.Series(np.nan, index=candidates.index)
    )
    lower_bound_fold = (t1_fraction / zero_upper).where(
        t2_fraction.le(0) & zero_upper.gt(0)
    )
candidates["fold_contraction_or_lower_bound"] = observed_fold.where(
    observed_fold.notna(), lower_bound_fold
)
candidates["fold_contraction_is_lower_bound"] = (
    observed_fold.isna() & lower_bound_fold.notna()
)


# Keep the specified candidate schema; calculation fields stay in the
# supplementary full clone master.
def as_boolean(value):
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def human_fold_contraction(value, is_lower_bound):
    if pd.isna(value):
        return ""
    rendered = f"{float(value):.7g}"
    return f"≥{rendered}" if as_boolean(is_lower_bound) else rendered


candidate_export = pd.DataFrame(
    {
        "Patient": candidates["patient"].astype(str),
        "Nucleotide clonotype ID": candidates["clone_id_nt"].astype(str),
        "CDR3 amino acid sequence": candidates["cdr3aa"].astype(str),
        "TRBV gene": candidates["TRBV"].astype(str),
        "TRBJ gene": candidates["TRBJ"].astype(str),
        "T1 frequency (%)": pd.to_numeric(
            candidates["T1_fraction_percent"], errors="coerce"
        ),
        "T2 frequency (%)": pd.to_numeric(
            candidates["T2_fraction_percent"], errors="coerce"
        ),
        "Fate at 0.10%": candidates["fate_0_10_percent"].replace(
            {
                "not detected": "Not detected at T2",
                "retained expanded": "Retained at or above 0.10% in T2",
                "detected below 0.10 percent": "Detected below 0.10% in T2",
            }
        ),
        "T1 to T2 fold contraction": [
            human_fold_contraction(value, is_lower_bound)
            for value, is_lower_bound in zip(
                candidates["fold_contraction_or_lower_bound"],
                candidates["fold_contraction_is_lower_bound"],
            )
        ],
        "Loss robust to sampling": candidates["sampling_robust_not_detected"].map(
            lambda value: "Yes" if as_boolean(value) else "No"
        ),
        "Public or private": candidates["public_private"].astype(str).str.capitalize(),
        "Control prevalence (%)": 100
        * pd.to_numeric(candidates["control_prevalence"], errors="coerce"),
        "Known antigen": candidates["known_antigen_final"],
        "Motif cluster": candidates["motif_cluster"].astype(str),
        "T3 fate": candidates["T3_fate"].replace(
            {
                "not available": "T3 not available",
                "expanded at T3": "At or above 0.10% at T3",
                "detected below 0.10 percent at T3": "Detected below 0.10% at T3",
                "not detected at T3": "Not detected at T3",
            }
        ),
        "Candidate class": candidates["candidate_interpretation"].replace(
            {
                "contracted clone with unknown specificity": "Contracted clone with unknown antigen specificity",
                "stable circulating clone": "Stable circulating clone",
                "candidate developing or persisting measles response": (
                    "Candidate for a developing or persistent measles response"
                ),
            }
        ),
    }
)
candidate_export.sort_values(
    ["Patient", "T1 frequency (%)"], ascending=[True, False]
).to_csv(FINAL / "candidate_clones.csv", index=False)

master.to_csv(SUPP / "01_full_clone_master.csv", index=False)
raw_wide = raw_metrics.pivot(
    index=["patient", "age_group"],
    columns="timepoint",
    values=["simpson_D", "simpson_clonality", "normalized_richness", "top_gt1_mass"],
)
raw_wide.columns = [f"{metric}_{point}" for metric, point in raw_wide.columns]
patient_summary = rare_patient.merge(
    raw_wide.reset_index(), on=["patient", "age_group"], how="left"
)
patient_summary.to_csv(SUPP / "02_patient_summary.csv", index=False)

extended_path = RESULTS / "expanded_threshold_sensitivity_extended.csv"
if not extended_path.exists():
    raise FileNotFoundError(extended_path)
threshold_final = pd.read_csv(extended_path)
threshold_final.to_csv(SUPP / "03_threshold_sensitivity.csv", index=False)
enrichment.to_csv(SUPP / "04_feature_enrichment.csv", index=False)

validation_parts = []
for label, filename in [
    ("identity", "repertoire_identity_qc.csv"),
    ("T3", "extended_T3_clone_validation.csv"),
    ("GA1_003_batch", "GA1_003_cross_batch_sensitivity.csv"),
]:
    p = CURRENT / filename
    if p.exists() and p.stat().st_size > 1:
        x = pd.read_csv(p)
        x.insert(0, "validation_type", label)
        validation_parts.append(x)
validation = (
    pd.concat(validation_parts, ignore_index=True, sort=False)
    if validation_parts
    else pd.DataFrame()
)
validation.to_csv(SUPP / "05_validation.csv", index=False)
convergence.to_csv(SUPP / "06_downsampling_iteration_convergence.csv", index=False)

# Publish PNG figures in the final cohort directories.
for d in [
    "group",
    "main",
    "downsampling_pairwise",
    "sensitivity_thresholds",
    "supplementary",
    "rarefied",
    "simpson_raw",
    "validation",
    "patients",
    "adults",
    "children",
]:
    p = FIG / d
    if p.exists():
        shutil.rmtree(p)
for d in [FIG / "patients", FIG / "rarefied", FIG / "simpson_raw", FIG / "validation"]:
    d.mkdir(parents=True, exist_ok=True)
for cohort in ["adults", "children"]:
    for section in ["patients", "rarefied", "simpson_raw", "validation"]:
        (FIG / cohort / section).mkdir(parents=True, exist_ok=True)
patient_src = FIG / "current" / "patients"
if patient_src.exists():
    for p in sorted(patient_src.glob("*.png")):
        shutil.copy2(p, FIG / "patients" / p.name)
        patient = next(
            (
                value
                for value in sorted(
                    master.patient.astype(str).unique(), key=len, reverse=True
                )
                if p.stem.startswith(value)
            ),
            None,
        )
        if patient is None:
            continue
        ages = (
            master.loc[master.patient.astype(str).eq(patient), "age_group"]
            .dropna()
            .astype(str)
        )
        if not ages.empty and ages.iloc[0] in {"adult", "child"}:
            cohort = "adults" if ages.iloc[0] == "adult" else "children"
            shutil.copy2(p, FIG / cohort / "patients" / p.name)


def label_for_clone(row):
    return f"{row.patient}  {short_seq(row.cdr3aa, 16)}  {row.TRBV}"


def mean_by_age(df, value, age):
    return float(
        pd.to_numeric(df.loc[df.age_group.eq(age), value], errors="coerce").mean()
    )


def stacked_fate(ax, sens, title, suffix, age_groups=("child", "adult")):
    z = sens[
        np.isclose(pd.to_numeric(sens.threshold_percent, errors="coerce"), 0.10)
    ].copy()
    ages = list(age_groups)
    xs, bottom = np.arange(len(ages)), np.zeros(len(ages))
    specs = [
        (
            "retained_mass_percent" + suffix,
            "Retained expanded\nT2 ≥ selected threshold",
            "Retained expanded",
        ),
        (
            "below_mass_percent" + suffix,
            "Detected below threshold\n0 < T2 < selected threshold",
            "Detected below threshold",
        ),
        (
            "not_detected_mass_percent" + suffix,
            "Not detected\nzero observed T2 reads",
            "Not detected",
        ),
    ]
    for col, legend_label, color_key in specs:
        values = np.array([mean_by_age(z, col, age) for age in ages], dtype=float)
        bars = ax.bar(
            xs,
            values,
            bottom=bottom,
            width=0.62,
            label=legend_label,
            color=FATE_COLORS[color_key],
            edgecolor="white",
            linewidth=0.8,
        )
        for bar, value, base in zip(bars, values, bottom):
            if np.isfinite(value) and value >= 5:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    base + value / 2,
                    f"{value:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=8,
                )
        bottom += np.nan_to_num(values)
    counts = [z.loc[z.age_group.eq(age), "patient"].nunique() for age in ages]
    age_names = {"child": "Children", "adult": "Adults"}
    ax.set_xticks(
        xs,
        [
            f"{age_names.get(age, age.title())}\n(n={count} paired)"
            for age, count in zip(ages, counts)
        ],
    )
    ax.set_ylim(0, 100)
    ax.set_ylabel("Mean share of each patient's T1-expanded mass (%)")
    ax.set_title(title)
    finish_axis(ax, "y")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=3,
        title="Fate of clones expanded at T1",
    )


def shared_scatter(ax, xcol, ycol, title, dircol, age_group=None):
    keep_age = (
        pd.Series(True, index=master.index)
        if age_group is None
        else master.age_group.eq(age_group)
    )
    z = master[
        primary_mask & keep_age & pd.to_numeric(master[ycol], errors="coerce").gt(0)
    ].copy()
    x = 100 * pd.to_numeric(z[xcol], errors="coerce").to_numpy(float)
    y = 100 * pd.to_numeric(z[ycol], errors="coerce").to_numpy(float)
    direction = z[dircol].astype(str)
    colors = {
        "increased": COLORS["teal"],
        "stable": COLORS["blue"],
        "decreased": COLORS["orange"],
    }
    legend_labels = {
        "increased": "Increased: T2/T1 ≥ 2",
        "stable": "Stable: 0.5 < T2/T1 < 2",
        "decreased": "Decreased: T2/T1 ≤ 0.5",
    }
    for label in ["increased", "stable", "decreased"]:
        keep = direction.eq(label).to_numpy()
        ax.scatter(
            x[keep],
            y[keep],
            s=25,
            alpha=0.56,
            color=colors[label],
            edgecolor="white",
            linewidth=0.35,
            label=legend_labels[label],
        )
    positive = np.concatenate([x[x > 0], y[y > 0]])
    lo = max(float(np.nanmin(positive)) * 0.72, 1e-4) if len(positive) else 0.01
    hi = float(np.nanmax(positive)) * 1.35 if len(positive) else 1.0
    ax.plot([lo, hi], [lo, hi], linestyle="--", color=COLORS["gray"], linewidth=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(
        xlim=(lo, hi),
        ylim=(lo, hi),
        xlabel="T1 repertoire fraction (%)",
        ylabel="T2 repertoire fraction (%)",
        title=title,
    )
    finish_axis(ax, "both")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=3,
        title="Direction among clones detected at both time points",
    )
    labels = np.array(
        [label_for_clone(r) for r in z.itertuples(index=False)], dtype=object
    )
    annotate_ranked_points(ax, x, y, labels, np.argsort(-np.maximum(x, y))[:10])


def plot_threshold_lines(ax, sens, branches, age_groups=("child", "adult")):
    names = {
        "pairwise_rarefied": "Repeated rarefaction",
        "simpson_no_downsampling": "Observed frequency",
    }
    for branch in branches:
        q_branch = sens[sens.branch.eq(branch)].copy()
        groups = (
            [
                (age, {"child": "Children", "adult": "Adults"}.get(age, age.title()))
                for age in age_groups
            ]
            if age_groups
            else [(None, "All paired patients")]
        )
        for age, age_label in groups:
            q = q_branch if age is None else q_branch[q_branch.age_group.eq(age)]
            g = (
                q.groupby("threshold_percent")["not_detected_mass_percent_median"]
                .mean()
                .reindex(THRESHOLD_PCT)
            )
            branch_color = AGE_COLORS[age] if age else BRANCH_COLORS[branch]
            ax.plot(
                g.index,
                g.values,
                marker="o" if branch == "pairwise_rarefied" else "s",
                ms=3.5,
                linestyle="-" if branch == "pairwise_rarefied" else "--",
                markerfacecolor=(
                    branch_color if branch == "pairwise_rarefied" else "white"
                ),
                markeredgecolor=branch_color,
                color=branch_color,
                zorder=3 if branch == "simpson_no_downsampling" else 2,
                label=f"{names[branch]}: {age_label}" if age_groups else names[branch],
            )
            missing = g[g.isna()].index
            if len(missing):
                ax.scatter(
                    missing,
                    np.zeros(len(missing)),
                    marker="x",
                    s=26,
                    color=branch_color,
                    zorder=4,
                )
                first = float(missing[0])
                ax.annotate(
                    "n=0 eligible",
                    (first, 0),
                    xytext=(3, 7),
                    textcoords="offset points",
                    fontsize=7,
                    color=branch_color,
                )
    ax.axvline(0.10, color=COLORS["ink"], linestyle=":", linewidth=1.2)
    ax.text(
        0.10,
        0.98,
        " selected 0.10%",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=8,
    )
    ax.set_xticks(THRESHOLD_PCT)
    ax.tick_params(axis="x", rotation=45)
    ax.set_xlabel("T1 expansion threshold (%)")
    ax.set_ylabel("Mean T1 mass not detected at T2 (%)")
    finish_axis(ax, "both")


# Rarefied 01: shared clones and cohort-mean fate.
fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.1), constrained_layout=True)
shared_scatter(
    axes[0],
    "T1_fraction_rarefied_median",
    "T2_fraction_rarefied_median",
    "Shared T1-expanded clones after repeated rarefaction\nTop 10 by maximum T1/T2 fraction are labelled",
    "rarefied_direction_0_10_percent",
)
stacked_fate(
    axes[1],
    threshold_final[threshold_final.branch.eq("pairwise_rarefied")],
    "Cohort mean clone fate at the selected 0.10% threshold\nEach patient contributes one T1–T2 repertoire pair",
    "_median",
)
panel_label(axes[0], "A")
panel_label(axes[1], "B")
fig.suptitle(
    "Primary sampling-controlled analysis",
    fontsize=15,
    fontweight="bold",
    color=COLORS["ink"],
)
fig.savefig(FIG / "rarefied" / "01_shared_and_fate.png", bbox_inches="tight")
plt.close(fig)

# Rarefied 02: expected versus observed detection, plus threshold sensitivity.
fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.8), constrained_layout=True)
z = master[primary_mask].copy()
x = pd.to_numeric(z["P_detection_if_unchanged_empirical"], errors="coerce").to_numpy(
    float
)
y = pd.to_numeric(
    z["P_detection_in_observed_T2_rarefactions"], errors="coerce"
).to_numpy(float)
robust = (
    truthy(z["sampling_robust_not_detected"])
    if "sampling_robust_not_detected" in z.columns
    else pd.Series(False, index=z.index)
)
axes[0].fill_between([0.8, 1.0], 0, 0.2, color=COLORS["orange"], alpha=0.10)
axes[0].scatter(
    x[~robust.to_numpy()],
    y[~robust.to_numpy()],
    s=24,
    alpha=0.45,
    color=COLORS["blue"],
    edgecolor="white",
    linewidth=0.3,
    label="Other T1-expanded clones",
)
axes[0].scatter(
    x[robust.to_numpy()],
    y[robust.to_numpy()],
    marker="*",
    s=85,
    color=COLORS["orange"],
    edgecolor=COLORS["ink"],
    linewidth=0.35,
    label="Sampling-robust T2 absence:\nexpected detection ≥0.99, observed 0",
)
axes[0].plot([0, 1], [0, 1], linestyle="--", color=COLORS["gray"], linewidth=1)
axes[0].set(
    xlim=(0, 1.01),
    ylim=(0, 1.01),
    xlabel="Expected detection probability if T1 abundance were unchanged",
    ylabel="Observed T2 detection probability across rarefactions",
    title="Internal sampling control\nBottom-right clones disappear more often than depth alone predicts",
)
finish_axis(axes[0], "both")
axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
labels = np.array([label_for_clone(r) for r in z.itertuples(index=False)], dtype=object)
annotate_ranked_points(
    axes[0],
    x,
    y,
    labels,
    np.argsort(-np.nan_to_num(x - y, nan=-np.inf))[:10],
    color=COLORS["orange"],
)
plot_threshold_lines(axes[1], threshold_final, ["pairwise_rarefied"])
axes[1].set_title(
    "Sensitivity to the expansion threshold\n× at zero means no eligible patient, not zero clone loss"
)
axes[1].legend(loc="best")
panel_label(axes[0], "A")
panel_label(axes[1], "B")
fig.suptitle(
    "Sampling robustness and threshold sensitivity",
    fontsize=15,
    fontweight="bold",
    color=COLORS["ink"],
)
fig.savefig(FIG / "rarefied" / "02_sampling_and_sensitivity.png", bbox_inches="tight")
plt.close(fig)


def dumbbell_and_trbv(branch, xcol, ycol, dircol, outfile, age_group=None):
    keep_age = (
        pd.Series(True, index=master.index)
        if age_group is None
        else master.age_group.eq(age_group)
    )
    paired = master[
        primary_mask
        & master.patient.astype(str).isin(paired_patients)
        & keep_age
        & pd.to_numeric(master[ycol], errors="coerce").gt(0)
    ].copy()
    log_col = (
        "rarefied_log2_T2_over_T1" if branch == "rarefied" else "raw_log2_T2_over_T1"
    )
    blocks = []
    for category in ["increased", "stable", "decreased"]:
        q = paired[paired[dircol].astype(str).eq(category)].copy()
        q[log_col] = pd.to_numeric(q[log_col], errors="coerce")
        if category == "increased":
            q = q.sort_values(log_col, ascending=False).head(5)
        elif category == "decreased":
            q = q.sort_values(log_col, ascending=True).head(5)
        else:
            q = q.assign(abs_change=q[log_col].abs()).sort_values("abs_change").head(5)
        q = q.reset_index(drop=True)
        q["rank_within_direction"] = np.arange(1, len(q) + 1)
        q["direction_block"] = category
        blocks.append(q)
    top = pd.concat(blocks, ignore_index=True) if blocks else pd.DataFrame()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15.5, 7.8),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.62, 1]},
    )
    if top.empty:
        axes[0].text(
            0.5, 0.5, "No nonzero paired T1–T2 candidates", ha="center", va="center"
        )
        axes[0].axis("off")
    else:
        y_pos = np.arange(len(top))
        t1 = 100 * pd.to_numeric(top[xcol], errors="coerce").to_numpy(float)
        t2 = 100 * pd.to_numeric(top[ycol], errors="coerce").to_numpy(float)
        for yi, a, b in zip(y_pos, t1, t2):
            axes[0].plot(
                [a, b], [yi, yi], color=COLORS["light_gray"], linewidth=2.2, zorder=1
            )
        axes[0].scatter(t1, y_pos, color=TIME_COLORS["T1"], s=34, label="T1", zorder=3)
        axes[0].scatter(t2, y_pos, color=TIME_COLORS["T2"], s=34, label="T2", zorder=3)
        labels = [
            f"{r.direction_block.title()} {int(r.rank_within_direction)}  |  "
            f"{r.patient}  {short_seq(r.cdr3aa, 15)}  {r.TRBV}"
            for r in top.itertuples(index=False)
        ]
        axes[0].set_yticks(y_pos, labels)
        axes[0].invert_yaxis()
        axes[0].set_xlabel("Repertoire fraction (%)")
        axes[0].set_title(
            "Top five paired clones in each direction\nRanks run from 1 to 5 within each block"
        )
        for boundary in [4.5, 9.5]:
            if boundary < len(top):
                axes[0].axhline(
                    boundary, color=COLORS["gray"], linestyle="--", linewidth=0.9
                )
        for start, color in zip(
            [0, 5, 10], [COLORS["teal"], COLORS["blue"], COLORS["orange"]]
        ):
            if start < len(top):
                axes[0].axhspan(
                    start - 0.5,
                    min(start + 4.5, len(top) - 0.5),
                    color=color,
                    alpha=0.045,
                    zorder=0,
                )
        axes[0].legend(loc="lower right", title="Paired time points")
        finish_axis(axes[0], "x")

    decreased = paired[paired[dircol].astype(str).eq("decreased")].copy()
    if decreased.empty:
        axes[1].text(
            0.5, 0.5, "No decreased nonzero clonotypes", ha="center", va="center"
        )
        axes[1].axis("off")
    else:
        g = (
            decreased.groupby("TRBV")[xcol]
            .sum()
            .sort_values(ascending=False)
            .head(10)
            .sort_values()
        )
        axes[1].barh(g.index, 100 * g.values, color=COLORS["orange"], alpha=0.88)
        axes[1].set_xlabel("T1 repertoire fraction carried by decreased clones (%)")
        axes[1].set_title(
            "TRBV context of the decreased clone pool\nTop 10 genes by cumulative T1 fraction"
        )
        finish_axis(axes[1], "x")
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    name = (
        "Repeated rarefaction"
        if branch == "rarefied"
        else "Observed frequencies without downsampling"
    )
    fig.suptitle(
        f"Candidate trajectories: {name}",
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.text(
        0.5,
        -0.01,
        "Direction definitions: increased = T2/T1 ≥ 2; stable = 0.5 < T2/T1 < 2; "
        "decreased = T2/T1 ≤ 0.5 but detected at T2. Not-detected clones have zero T2 reads "
        "and are not shown because a T2 dumbbell endpoint does not exist.",
        ha="center",
        va="top",
        fontsize=8.0,
        color=COLORS["gray"],
        wrap=True,
    )
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


dumbbell_and_trbv(
    "rarefied",
    "T1_fraction_rarefied_median",
    "T2_fraction_rarefied_median",
    "rarefied_direction_0_10_percent",
    FIG / "rarefied" / "03_candidate_dumbbell_and_trbv.png",
)


def enrichment_plot(
    branch, outfile, enrichment_frame=None, cohort_label="All paired patients"
):
    source = enrichment if enrichment_frame is None else enrichment_frame
    e = (
        source[source.branch.astype(str).eq(branch)].copy()
        if not source.empty and "branch" in source.columns
        else pd.DataFrame()
    )
    if not e.empty:
        e = e[pd.to_numeric(e.get("fdr", np.nan), errors="coerce").le(0.05)].copy()
    fig, ax = plt.subplots(figsize=(10.8, 6.9) if not e.empty else (10.0, 4.8))
    fig.subplots_adjust(left=0.18, right=0.98, top=0.78, bottom=0.27)
    if e.empty:
        ax.text(
            0.5,
            0.54,
            "No TRBV, TRBJ or CDR3 3-mer feature passes BH FDR ≤ 0.05",
            ha="center",
            va="center",
            fontsize=11,
            transform=ax.transAxes,
        )
        ax.text(
            0.5,
            0.45,
            "The empty result is retained because it is the statistical conclusion.",
            ha="center",
            va="center",
            fontsize=9,
            color=COLORS["gray"],
            transform=ax.transAxes,
        )
        ax.set_axis_off()
    else:
        e = e.assign(
            abscoef=pd.to_numeric(
                e.log_odds_affected_vs_retained, errors="coerce"
            ).abs()
        )
        e = (
            e.sort_values("abscoef", ascending=False)
            .head(12)
            .sort_values("log_odds_affected_vs_retained")
        )
        y_pos = np.arange(len(e))
        coef = pd.to_numeric(e.log_odds_affected_vs_retained, errors="coerce").to_numpy(
            float
        )
        se = pd.to_numeric(e.get("standard_error", np.nan), errors="coerce").to_numpy(
            float
        )
        ax.errorbar(
            coef,
            y_pos,
            xerr=1.96 * se,
            fmt="none",
            ecolor=COLORS["gray"],
            elinewidth=1.2,
            capsize=2.5,
        )
        ax.scatter(
            coef,
            y_pos,
            c=np.where(coef >= 0, COLORS["orange"], COLORS["teal"]),
            s=42,
            zorder=3,
        )
        ax.axvline(0, color=COLORS["ink"], linewidth=0.9)
        feature_names = {
            "CDR3_3mer": "CDR3 aa 3-mer",
            "TRBV": "TRBV gene",
            "TRBJ": "TRBJ gene",
        }
        ax.set_yticks(
            y_pos,
            [
                f"{feature_names.get(a, a)}: {b}"
                for a, b in zip(e.feature_type, e.feature)
            ],
        )
        ax.set_xlabel("Adjusted log odds ratio for contracted or lost fate")
        finish_axis(ax, "x")
        legend_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=COLORS["orange"],
                markeredgecolor="none",
                markersize=7,
                label="Enriched in contracted or lost clones",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=COLORS["teal"],
                markeredgecolor="none",
                markersize=7,
                label="Enriched among retained clones",
            ),
            Line2D(
                [0],
                [0],
                color=COLORS["gray"],
                marker="o",
                markerfacecolor="white",
                markeredgecolor=COLORS["gray"],
                linewidth=1.2,
                markersize=5,
                label="Adjusted log odds ratio and 95% CI",
            ),
        ]
        ax.legend(
            handles=legend_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.015),
            ncol=3,
            fontsize=8.0,
            columnspacing=1.4,
            handlelength=2.0,
        )
    method = (
        "Patient-clustered GEE after repeated rarefaction"
        if branch == "pairwise_rarefied"
        else "Patient-clustered GEE using observed frequencies"
    )
    outcome_definition = (
        "Contracted: ≥5-fold decrease from T1 to T2.\n"
        "Lost: absent at T2 despite ≥0.99 expected detection if unchanged.\n"
        "Retained: detected at T2 at ≥0.10%."
        if branch == "pairwise_rarefied"
        else "Contracted: ≥5-fold decrease from T1 to T2.\n"
        "Lost: zero observed reads at T2.\n"
        "Retained: detected at T2 at ≥0.10%."
    )
    fig.suptitle(
        f"Sequence-feature enrichment by clone fate: {cohort_label}\n{method}",
        y=0.975,
        fontsize=14.5,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.add_artist(
        Line2D(
            [0.18, 0.98],
            [0.205, 0.205],
            transform=fig.transFigure,
            color=COLORS["light_gray"],
            linewidth=0.8,
        )
    )
    fig.text(
        0.18,
        0.185,
        "Outcome definitions",
        ha="left",
        va="top",
        fontsize=8.2,
        fontweight="semibold",
        color=COLORS["ink"],
    )
    fig.text(
        0.18,
        0.155,
        outcome_definition,
        ha="left",
        va="top",
        fontsize=7.7,
        color=COLORS["gray"],
        linespacing=1.35,
    )
    fig.text(
        0.60,
        0.185,
        "Statistical analysis",
        ha="left",
        va="top",
        fontsize=8.2,
        fontweight="semibold",
        color=COLORS["ink"],
    )
    fig.text(
        0.60,
        0.155,
        "Points show adjusted log odds ratios; horizontal bars show 95% confidence intervals.\n"
        "Binomial GEE clusters clones by patient and adjusts for log10-transformed\n"
        "T1 clone frequency. Up to 12 features with BH FDR ≤0.05 are shown.",
        ha="left",
        va="top",
        fontsize=7.7,
        color=COLORS["gray"],
        linespacing=1.35,
    )
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


enrichment_plot(
    "pairwise_rarefied", FIG / "rarefied" / "04_enrichment_and_specificity.png"
)


def simpson_and_fate(
    outfile, age_groups=("child", "adult"), cohort_label="All paired patients"
):
    """Keep only the non-redundant Simpson D panel plus clone fate."""
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.8), constrained_layout=True)
    q_metrics = raw_metrics[raw_metrics.age_group.isin(age_groups)].copy()
    extreme = []
    for patient, q in q_metrics.groupby("patient"):
        if not {"T1", "T2"}.issubset(set(q.timepoint)):
            continue
        a = float(q.loc[q.timepoint.eq("T1"), "simpson_D"].iloc[0])
        b = float(q.loc[q.timepoint.eq("T2"), "simpson_D"].iloc[0])
        age = str(q.age_group.iloc[0])
        color = AGE_COLORS.get(age, COLORS["gray"])
        axes[0].plot([1, 2], [a, b], color=color, alpha=0.38, linewidth=1.1)
        axes[0].scatter([1, 2], [a, b], color=color, s=19, alpha=0.78)
        extreme.append((abs(b - a), str(patient), b))
    for rank, (_, patient, b) in enumerate(sorted(extreme, reverse=True)[:3]):
        axes[0].annotate(
            patient,
            xy=(2, b),
            xycoords="data",
            xytext=(0.74, 0.88 - 0.11 * rank),
            textcoords="axes fraction",
            fontsize=7.5,
            ha="left",
            va="center",
            arrowprops={"arrowstyle": "-", "color": COLORS["gray"], "lw": 0.7},
            bbox={
                "facecolor": "white",
                "edgecolor": COLORS["light_gray"],
                "boxstyle": "round,pad=0.15",
                "alpha": 0.9,
            },
        )
    axes[0].set_xticks([1, 2], ["T1\nearly infection", "T2\nrecovery"])
    axes[0].set_ylabel("Simpson D = Σ(clone fraction²)")
    axes[0].set_title(
        "Repertoire concentration without downsampling\nHigher Simpson D means stronger dominance by large clones"
    )
    finish_axis(axes[0], "y")
    delta = q_metrics.pivot(
        index="patient", columns="timepoint", values="simpson_D"
    ).dropna()
    med = float((delta["T2"] - delta["T1"]).median()) if not delta.empty else np.nan
    axes[0].text(
        0.02,
        0.04,
        f"Median paired change: {med:+.3g}",
        transform=axes[0].transAxes,
        fontsize=8,
    )
    handles = [
        Line2D(
            [0],
            [0],
            color=AGE_COLORS[age],
            marker="o",
            label={"child": "Children", "adult": "Adults"}[age],
        )
        for age in age_groups
    ]
    axes[0].legend(handles=handles, loc="best")
    stacked_fate(
        axes[1],
        threshold_final[threshold_final.branch.eq("simpson_no_downsampling")],
        "Observed-frequency clone fate at 0.10%\nMean of patients; full libraries, no downsampling",
        "_median",
        age_groups=age_groups,
    )
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle(
        f"Observed-frequency branch: {cohort_label}",
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


simpson_and_fate(FIG / "simpson_raw" / "01_simpson_metrics_and_fate.png")

fig, axes = plt.subplots(1, 2, figsize=(14.4, 5.9), constrained_layout=True)
z = master[
    primary_mask
    & pd.to_numeric(master["rarefied_log2_T2_over_T1"], errors="coerce").notna()
    & pd.to_numeric(master["raw_log2_T2_over_T1"], errors="coerce").notna()
].copy()
x = pd.to_numeric(z["rarefied_log2_T2_over_T1"], errors="coerce").to_numpy(float)
y = pd.to_numeric(z["raw_log2_T2_over_T1"], errors="coerce").to_numpy(float)
same = (
    z["rarefied_direction_0_10_percent"]
    .astype(str)
    .eq(z["raw_direction_0_10_percent"].astype(str))
    .to_numpy()
)
axes[0].scatter(
    x[same],
    y[same],
    s=25,
    color=COLORS["teal"],
    alpha=0.48,
    edgecolor="white",
    linewidth=0.3,
    label="Same direction class",
)
axes[0].scatter(
    x[~same],
    y[~same],
    s=34,
    color=COLORS["orange"],
    alpha=0.78,
    edgecolor="white",
    linewidth=0.3,
    label="Different direction class",
)
limit = (
    max(1.0, float(np.nanmax(np.abs(np.concatenate([x, y]))))) * 1.08 if len(x) else 1.0
)
axes[0].plot(
    [-limit, limit], [-limit, limit], linestyle="--", color=COLORS["gray"], linewidth=1
)
axes[0].axhline(0, color=COLORS["light_gray"], linewidth=0.8)
axes[0].axvline(0, color=COLORS["light_gray"], linewidth=0.8)
axes[0].set(
    xlim=(-limit, limit),
    ylim=(-limit, limit),
    xlabel="Repeated rarefaction: log2(T2/T1)",
    ylabel="Observed frequency: log2(T2/T1)",
    title="Clone-level agreement between branches\nTop 10 absolute disagreements are labelled",
)
finish_axis(axes[0], "both")
axes[0].legend(loc="lower right")
labels = np.array([label_for_clone(r) for r in z.itertuples(index=False)], dtype=object)
annotate_ranked_points(
    axes[0], x, y, labels, np.argsort(-np.abs(x - y))[:10], color=COLORS["orange"]
)
plot_threshold_lines(
    axes[1], threshold_final, ["pairwise_rarefied", "simpson_no_downsampling"]
)
axes[1].set_title(
    "Threshold sensitivity in both branches\nSolid: repeated rarefaction; dashed: observed frequency"
)
axes[1].legend(loc="best", fontsize=7.6)
panel_label(axes[0], "A")
panel_label(axes[1], "B")
fig.suptitle(
    "Primary and complementary branches: similarities and differences",
    fontsize=15,
    fontweight="bold",
    color=COLORS["ink"],
)
fig.savefig(FIG / "simpson_raw" / "02_shared_and_sensitivity.png", bbox_inches="tight")
plt.close(fig)

dumbbell_and_trbv(
    "raw",
    "T1_fraction",
    "T2_fraction",
    "raw_direction_0_10_percent",
    FIG / "simpson_raw" / "03_candidate_dumbbell_and_trbv.png",
)
enrichment_plot(
    "simpson_no_downsampling", FIG / "simpson_raw" / "04_enrichment_and_specificity.png"
)

old_val = FIG / "current" / "26_validation_identity_T3_batch.png"
if old_val.exists():
    shutil.copy2(old_val, FIG / "validation" / "01_identity_T3_batch.png")

direction_order = ["increased", "stable", "decreased", "not detected"]
z = master[primary_mask].copy()
tab = pd.crosstab(
    pd.Categorical(z["rarefied_direction_0_10_percent"], categories=direction_order),
    pd.Categorical(z["raw_direction_0_10_percent"], categories=direction_order),
    dropna=False,
)
row_pct = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0) * 100
fig, ax = plt.subplots(figsize=(8.2, 6.6), constrained_layout=True)
im = ax.imshow(row_pct.to_numpy(float), vmin=0, vmax=100, cmap="Blues", aspect="auto")
for i in range(tab.shape[0]):
    for j in range(tab.shape[1]):
        value, count = row_pct.iloc[i, j], int(tab.iloc[i, j])
        ax.text(
            j,
            i,
            f"{count}\n{value:.1f}%" if np.isfinite(value) else f"{count}\nNA",
            ha="center",
            va="center",
            color="white" if np.isfinite(value) and value >= 55 else COLORS["ink"],
            fontsize=9,
        )
ax.set_xticks(
    np.arange(len(direction_order)),
    [x.title() for x in direction_order],
    rotation=25,
    ha="right",
)
ax.set_yticks(np.arange(len(direction_order)), [x.title() for x in direction_order])
ax.set_xlabel(
    "Complementary branch: observed frequencies\n(Simpson metrics, no downsampling)"
)
ax.set_ylabel("Primary branch: repeated pair-specific rarefaction")
ax.set_title(
    "Clone-direction concordance at the 0.10% threshold\nCells show clone count and row percentage"
)
fig.colorbar(im, ax=ax, shrink=0.8, label="Percentage within rarefaction direction (%)")
fig.savefig(FIG / "validation" / "02_branch_concordance.png", bbox_inches="tight")
plt.close(fig)


def threshold_selection_plot(
    outfile, age_group=None, cohort_label="All paired patients"
):
    source = (
        threshold_final
        if age_group is None
        else threshold_final[threshold_final.age_group.eq(age_group)]
    )
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.7), constrained_layout=True)
    plot_threshold_lines(
        axes[0], source, ["pairwise_rarefied", "simpson_no_downsampling"], age_groups=()
    )
    axes[0].set_title(
        "T1 mass later not detected\nCircle/solid: rarefaction; hollow square/dashed: observed frequencies"
    )
    axes[0].legend(loc="best")
    for branch, label in [
        ("pairwise_rarefied", "Repeated rarefaction"),
        ("simpson_no_downsampling", "Observed frequency"),
    ]:
        q = source[source.branch.eq(branch)]
        col = (
            "n_T1_expanded_median"
            if "n_T1_expanded_median" in q.columns
            else "n_T1_expanded"
        )
        g = q.groupby("threshold_percent")[col].median().reindex(THRESHOLD_PCT)
        color = BRANCH_COLORS[branch]
        axes[1].plot(
            g.index,
            g.values,
            marker="o" if branch == "pairwise_rarefied" else "s",
            markerfacecolor=color if branch == "pairwise_rarefied" else "white",
            linestyle="-" if branch == "pairwise_rarefied" else "--",
            ms=4,
            color=color,
            label=label,
        )
    for ax in axes:
        ax.axvline(0.10, color=COLORS["ink"], linestyle=":", linewidth=1.2)
        ax.axvspan(0.075, 0.125, color=COLORS["gold"], alpha=0.14)
        ax.text(
            0.10,
            0.98,
            " selected 0.10%",
            transform=ax.get_xaxis_transform(),
            ha="left",
            va="top",
            fontsize=8,
        )
        ax.set_xticks(THRESHOLD_PCT)
        ax.tick_params(axis="x", rotation=45)
    axes[1].set_xlabel("T1 expansion threshold (%)")
    axes[1].set_ylabel("Median T1-expanded clone count per patient")
    axes[1].set_title(
        "Candidate coverage\nAll thresholds, including 0.90%, are computed"
    )
    axes[1].legend(loc="best")
    finish_axis(axes[1], "both")
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle(
        f"Threshold choice for {cohort_label}: 0.10% primary cutoff",
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.text(
        0.5,
        -0.015,
        "0.10% is the pre-specified balance: it avoids the sampling-sensitive 0.05% edge "
        "and retains more candidates than stricter cutoffs. Overlapping lines mean the two "
        "methods give nearly the same estimate. Line style identifies the analytical branch; "
        "missing eligibility is explicitly marked n=0.",
        ha="center",
        va="top",
        fontsize=8.4,
    )
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


if "age_group" not in convergence.columns and "patient" in convergence.columns:
    patient_age = master[["patient", "age_group"]].drop_duplicates("patient")
    convergence = convergence.merge(patient_age, on="patient", how="left")


def iteration_convergence_plot(
    outfile, age_group=None, cohort_label="All paired patients"
):
    source = (
        convergence
        if age_group is None
        else convergence[convergence.age_group.eq(age_group)]
    )
    iter_means = (
        source.groupby("iteration", as_index=True)
        .agg(
            not_detected_mass_percent=("not_detected_mass_percent", "mean"),
            expanded_jaccard=("expanded_jaccard", "mean"),
        )
        .sort_index()
    )
    jaccard_percent = 100 * iter_means["expanded_jaccard"]
    series = [
        (
            "Primary result stability\nT1-expanded mass not detected at T2",
            iter_means["not_detected_mass_percent"],
            COLORS["orange"],
        ),
        (
            "Set-overlap stability\nExpanded-clone Jaccard at the 0.10% threshold",
            jaccard_percent,
            COLORS["teal"],
        ),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(14.8, 7.1))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.80, bottom=0.33, wspace=0.14)
    for ax, (title, values, color) in zip(axes, series):
        running_mean = values.expanding(min_periods=1).mean()
        n = pd.Series(np.arange(1, len(values) + 1), index=values.index, dtype=float)
        mc_se = values.expanding(min_periods=2).std() / np.sqrt(n)
        ci = (1.96 * mc_se).fillna(0.0)
        ax.scatter(
            values.index,
            values.values,
            s=8,
            color=color,
            alpha=0.12,
            label="Single equal-depth resample",
        )
        ax.fill_between(
            values.index.to_numpy(float),
            (running_mean - ci).to_numpy(float),
            (running_mean + ci).to_numpy(float),
            color=color,
            alpha=0.16,
            label="Running 95% Monte Carlo CI",
        )
        ax.plot(
            running_mean.index,
            running_mean.values,
            color=color,
            linewidth=1.8,
            label="Running mean after k iterations",
        )
        final_value = float(running_mean.iloc[-1])
        final_half_width = float(ci.iloc[-1])
        ax.axhline(final_value, color=COLORS["ink"], linestyle=":", linewidth=1)
        ax.text(
            0.98,
            0.05,
            f"Final mean: {final_value:.2f}%\n95% Monte Carlo half-width: "
            f"±{final_half_width:.2f} percentage points",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
        )
        ax.set(
            xlabel="Number of rarefaction iterations accumulated",
            ylabel="Cohort estimate (%)",
            title=title,
        )
        finish_axis(ax, "both")
        ax.legend(loc="best", fontsize=7.6)
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle(
        f"Monte Carlo convergence of pair-specific rarefaction: {cohort_label}",
        y=0.965,
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    final_jaccard = float(jaccard_percent.expanding(min_periods=1).mean().iloc[-1])
    fig.text(
        0.07,
        0.235,
        "Rarefaction and panel A\n"
        "Each iteration independently samples T1 and T2 without replacement to the patient's\n"
        "common depth; only the selected reads change. Panel A tracks the abundance-weighted\n"
        "share of T1-expanded signal that is not detected at T2.",
        ha="left",
        va="top",
        fontsize=8.0,
        color=COLORS["gray"],
        linespacing=1.35,
        bbox={
            "boxstyle": "round,pad=0.45",
            "facecolor": "#F7F9FA",
            "edgecolor": COLORS["light_gray"],
            "linewidth": 0.7,
        },
    )
    fig.text(
        0.535,
        0.235,
        "Jaccard definition and interpretation\n"
        "J = |E₁ ∩ E₂| / |E₁ ∪ E₂|; E₁ and E₂ are the T1 and T2 clone sets at ≥0.10%;\n"
        "the plotted value is 100×J. 0%: no shared clones; 100%: identical sets; higher: greater overlap.\n"
        "J is calculated per patient per resample; each pale point is that resample's cohort mean.\n"
        f"Final cohort mean J = {final_jaccard:.2f}%. It was selected for symmetric set comparison:\n"
        "T1 losses and T2 gains both reduce the score; the 0.10% cutoff limits sensitivity to the rare-clone tail.",
        ha="left",
        va="top",
        fontsize=8.0,
        color=COLORS["gray"],
        linespacing=1.35,
        bbox={
            "boxstyle": "round,pad=0.45",
            "facecolor": "#F7F9FA",
            "edgecolor": COLORS["light_gray"],
            "linewidth": 0.7,
        },
    )
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


threshold_selection_plot(FIG / "validation" / "03_threshold_selection.png")
iteration_convergence_plot(FIG / "validation" / "04_downsampling_iterations.png")


def sampling_and_sensitivity_plot(
    outfile, age_group=None, cohort_label="All paired patients"
):
    keep_age = (
        pd.Series(True, index=master.index)
        if age_group is None
        else master.age_group.eq(age_group)
    )
    z = master[primary_mask & keep_age].copy()
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.8), constrained_layout=True)
    x = pd.to_numeric(
        z["P_detection_if_unchanged_empirical"], errors="coerce"
    ).to_numpy(float)
    y = pd.to_numeric(
        z["P_detection_in_observed_T2_rarefactions"], errors="coerce"
    ).to_numpy(float)
    robust = (
        truthy(z["sampling_robust_not_detected"])
        if "sampling_robust_not_detected" in z.columns
        else pd.Series(False, index=z.index)
    )
    axes[0].fill_between([0.8, 1.0], 0, 0.2, color=COLORS["orange"], alpha=0.10)
    axes[0].scatter(
        x[~robust.to_numpy()],
        y[~robust.to_numpy()],
        s=24,
        alpha=0.45,
        color=COLORS["blue"],
        edgecolor="white",
        linewidth=0.3,
        label="Other T1-expanded clones",
    )
    axes[0].scatter(
        x[robust.to_numpy()],
        y[robust.to_numpy()],
        marker="*",
        s=85,
        color=COLORS["orange"],
        edgecolor=COLORS["ink"],
        linewidth=0.35,
        label="Sampling-robust T2 absence:\nexpected detection ≥0.99, observed 0",
    )
    axes[0].plot([0, 1], [0, 1], linestyle="--", color=COLORS["gray"], linewidth=1)
    axes[0].set(
        xlim=(0, 1.01),
        ylim=(0, 1.01),
        xlabel="Expected detection if the T1 clone stayed equally abundant",
        ylabel="Fraction of T2 rarefactions where the clone was detected",
        title="Depth-only expectation versus real T2 observations\nBottom-right: disappearance exceeds sampling expectation",
    )
    finish_axis(axes[0], "both")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    labels = np.array(
        [label_for_clone(r) for r in z.itertuples(index=False)], dtype=object
    )
    annotate_ranked_points(
        axes[0],
        x,
        y,
        labels,
        np.argsort(-np.nan_to_num(x - y, nan=-np.inf))[:10],
        color=COLORS["orange"],
    )
    sens = threshold_final[threshold_final.branch.eq("pairwise_rarefied")]
    if age_group is not None:
        sens = sens[sens.age_group.eq(age_group)]
    plot_threshold_lines(axes[1], sens, ["pairwise_rarefied"], age_groups=())
    axes[1].set_title(
        "Threshold sensitivity\n× at zero means no eligible patient, not zero clone loss"
    )
    axes[1].legend(loc="best")
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle(
        f"Sampling robustness and threshold sensitivity: {cohort_label}",
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


def branch_comparison_plot(outfile, age_group=None, cohort_label="All paired patients"):
    keep_age = (
        pd.Series(True, index=master.index)
        if age_group is None
        else master.age_group.eq(age_group)
    )
    z = master[
        primary_mask
        & keep_age
        & pd.to_numeric(master["rarefied_log2_T2_over_T1"], errors="coerce").notna()
        & pd.to_numeric(master["raw_log2_T2_over_T1"], errors="coerce").notna()
    ].copy()
    x = pd.to_numeric(z["rarefied_log2_T2_over_T1"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(z["raw_log2_T2_over_T1"], errors="coerce").to_numpy(float)
    same = (
        z["rarefied_direction_0_10_percent"]
        .astype(str)
        .eq(z["raw_direction_0_10_percent"].astype(str))
        .to_numpy()
    )
    fig, axes = plt.subplots(1, 2, figsize=(14.4, 5.9), constrained_layout=True)
    axes[0].scatter(
        x[same],
        y[same],
        s=25,
        color=COLORS["teal"],
        alpha=0.48,
        edgecolor="white",
        linewidth=0.3,
        label="Same direction class",
    )
    axes[0].scatter(
        x[~same],
        y[~same],
        s=34,
        color=COLORS["orange"],
        alpha=0.78,
        edgecolor="white",
        linewidth=0.3,
        label="Different direction class",
    )
    limit = (
        max(1.0, float(np.nanmax(np.abs(np.concatenate([x, y]))))) * 1.08
        if len(x)
        else 1.0
    )
    axes[0].plot(
        [-limit, limit],
        [-limit, limit],
        linestyle="--",
        color=COLORS["gray"],
        linewidth=1,
    )
    axes[0].axhline(0, color=COLORS["light_gray"], linewidth=0.8)
    axes[0].axvline(0, color=COLORS["light_gray"], linewidth=0.8)
    axes[0].set(
        xlim=(-limit, limit),
        ylim=(-limit, limit),
        xlabel="Repeated rarefaction: log2(T2/T1)",
        ylabel="Observed frequency: log2(T2/T1)",
        title="Clone-level agreement\nDiagonal means the two estimates are equal",
    )
    finish_axis(axes[0], "both")
    axes[0].legend(loc="lower right")
    sens = (
        threshold_final
        if age_group is None
        else threshold_final[threshold_final.age_group.eq(age_group)]
    )
    plot_threshold_lines(
        axes[1], sens, ["pairwise_rarefied", "simpson_no_downsampling"], age_groups=()
    )
    axes[1].set_title(
        "Threshold sensitivity by analytical branch\nFilled circle/solid: rarefaction; hollow square/dashed: observed"
    )
    axes[1].legend(loc="best", fontsize=7.6)
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle(
        f"Rarefaction versus observed-frequency branch: {cohort_label}",
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.text(
        0.5,
        -0.01,
        "Direction classes use the same fold-change boundaries: increased = T2/T1 ≥ 2; stable = "
        "0.5 < T2/T1 < 2; decreased = T2/T1 ≤ 0.5. Not detected means median rarefied T2 "
        "fraction is zero in the primary branch and zero observed T2 reads in the complementary branch.",
        ha="center",
        va="top",
        fontsize=8.0,
        color=COLORS["gray"],
        wrap=True,
    )
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


def concordance_plot(outfile, age_group=None, cohort_label="All paired patients"):
    keep_age = (
        pd.Series(True, index=master.index)
        if age_group is None
        else master.age_group.eq(age_group)
    )
    z = master[primary_mask & keep_age].copy()
    tab = pd.crosstab(
        pd.Categorical(
            z["rarefied_direction_0_10_percent"], categories=direction_order
        ),
        pd.Categorical(z["raw_direction_0_10_percent"], categories=direction_order),
        dropna=False,
    )
    row_pct = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0) * 100
    fig, ax = plt.subplots(figsize=(8.2, 6.6), constrained_layout=True)
    im = ax.imshow(
        row_pct.to_numpy(float), vmin=0, vmax=100, cmap="Blues", aspect="auto"
    )
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            value, count = row_pct.iloc[i, j], int(tab.iloc[i, j])
            ax.text(
                j,
                i,
                f"{count}\n{value:.1f}%" if np.isfinite(value) else f"{count}\nn/a",
                ha="center",
                va="center",
                color="white" if np.isfinite(value) and value >= 55 else COLORS["ink"],
                fontsize=9,
            )
    ax.set_xticks(
        np.arange(len(direction_order)),
        [x.title() for x in direction_order],
        rotation=25,
        ha="right",
    )
    ax.set_yticks(np.arange(len(direction_order)), [x.title() for x in direction_order])
    ax.set_xlabel("Observed-frequency direction class")
    ax.set_ylabel("Repeated-rarefaction direction class")
    ax.set_title(
        f"Direction agreement at 0.10%: {cohort_label}\nRows: rarefaction; columns: observed frequencies"
    )
    fig.colorbar(im, ax=ax, shrink=0.8, label="Percentage within rarefaction row (%)")
    fig.text(
        0.5,
        -0.01,
        "Increased: T2/T1 ≥ 2. Stable: 0.5 < T2/T1 < 2. Decreased: T2/T1 ≤ 0.5 "
        "with the clone still detected. Not detected: median rarefied T2 fraction is zero in the "
        "primary branch; zero observed T2 reads in the observed-frequency branch.",
        ha="center",
        va="top",
        fontsize=8.0,
        color=COLORS["gray"],
        wrap=True,
    )
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


sampling_and_sensitivity_plot(FIG / "rarefied" / "02_sampling_and_sensitivity.png")
branch_comparison_plot(FIG / "simpson_raw" / "02_shared_and_sensitivity.png")
concordance_plot(FIG / "validation" / "02_branch_concordance.png")


# Fit and plot each age group separately.
for age_group, cohort_dir, cohort_label in [
    ("adult", "adults", "Adults only"),
    ("child", "children", "Children only"),
]:
    base = FIG / cohort_dir
    age_groups = (age_group,)
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.1), constrained_layout=True)
    shared_scatter(
        axes[0],
        "T1_fraction_rarefied_median",
        "T2_fraction_rarefied_median",
        f"Shared T1-expanded clones: {cohort_label}\nTop 10 by maximum T1/T2 fraction are labelled",
        "rarefied_direction_0_10_percent",
        age_group=age_group,
    )
    stacked_fate(
        axes[1],
        threshold_final[threshold_final.branch.eq("pairwise_rarefied")],
        "Mean patient-level fate after repeated rarefaction\nSelected T1 threshold: 0.10%",
        "_median",
        age_groups=age_groups,
    )
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle(
        f"Primary sampling-controlled analysis: {cohort_label}",
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.savefig(base / "rarefied" / "01_shared_and_fate.png", bbox_inches="tight")
    plt.close(fig)

    sampling_and_sensitivity_plot(
        base / "rarefied" / "02_sampling_and_sensitivity.png", age_group, cohort_label
    )
    dumbbell_and_trbv(
        "rarefied",
        "T1_fraction_rarefied_median",
        "T2_fraction_rarefied_median",
        "rarefied_direction_0_10_percent",
        base / "rarefied" / "03_candidate_dumbbell_and_trbv.png",
        age_group,
    )
    subgroup_master = master[master.age_group.eq(age_group)].copy()
    subgroup_enrichment = pd.concat(
        [
            enrichment_for_branch(subgroup_master, "pairwise_rarefied"),
            enrichment_for_branch(subgroup_master, "simpson_no_downsampling"),
        ],
        ignore_index=True,
        sort=False,
    )
    enrichment_plot(
        "pairwise_rarefied",
        base / "rarefied" / "04_enrichment_and_specificity.png",
        subgroup_enrichment,
        cohort_label,
    )

    simpson_and_fate(
        base / "simpson_raw" / "01_simpson_metrics_and_fate.png",
        age_groups,
        cohort_label,
    )
    branch_comparison_plot(
        base / "simpson_raw" / "02_shared_and_sensitivity.png", age_group, cohort_label
    )
    dumbbell_and_trbv(
        "raw",
        "T1_fraction",
        "T2_fraction",
        "raw_direction_0_10_percent",
        base / "simpson_raw" / "03_candidate_dumbbell_and_trbv.png",
        age_group,
    )
    enrichment_plot(
        "simpson_no_downsampling",
        base / "simpson_raw" / "04_enrichment_and_specificity.png",
        subgroup_enrichment,
        cohort_label,
    )
    concordance_plot(
        base / "validation" / "02_branch_concordance.png", age_group, cohort_label
    )
    threshold_selection_plot(
        base / "validation" / "03_threshold_selection.png", age_group, cohort_label
    )
    iteration_convergence_plot(
        base / "validation" / "04_downsampling_iterations.png", age_group, cohort_label
    )

(FIG / "FIGURE_INDEX.txt").write_text(
    "\n".join(
        [
            "FINAL FIGURE SET",
            "",
            "patients: combined index of one panel per accepted longitudinal patient comparison",
            "adults/: adult-only copies of every analysis and patient panel",
            "children/: child-only copies of every analysis and patient panel",
            "rarefied/01_shared_and_fate.png: labelled shared clones and mean patient-level fate by age",
            "rarefied/02_sampling_and_sensitivity.png: expected-versus-observed detection and complete threshold grid",
            "rarefied/03_candidate_dumbbell_and_trbv.png: paired top-5 increased, stable and decreased blocks",
            "rarefied/04_enrichment_and_specificity.png: patient-clustered GEE feature enrichment only",
            "simpson_raw/01_simpson_metrics_and_fate.png: non-redundant Simpson D and raw fate",
            "simpson_raw/02_shared_and_sensitivity.png: comparison with the primary rarefaction branch",
            "simpson_raw/03_candidate_dumbbell_and_trbv.png: observed-frequency candidate blocks and TRBV context",
            "simpson_raw/04_enrichment_and_specificity.png: observed-frequency GEE enrichment only",
            "validation/01_identity_T3_batch.png: identity, T3 and batch validation if available",
            "validation/02_branch_concordance.png: rows rarefaction, columns observed frequency",
            "validation/03_threshold_selection.png: evidence and rationale for the 0.10% threshold",
            "validation/04_downsampling_iterations.png: running mean and 95% Monte Carlo CI",
            "",
            "Primary clone key: CDR3 nucleotide + TRBV + TRBJ",
            "Primary expansion threshold: 0.10 percent",
            "Sensitivity grid: 0.05, 0.10, then 0.20 to 1.50 percent",
            "Repeated rarefaction: at least 500 iterations at each pair-specific common depth",
            "Final result table: results/final/candidate_clones.csv",
            "Supporting tables: results/supplementary_tables",
        ]
    )
    + "\n",
    encoding="utf-8",
)

if (FIG / "current").exists():
    shutil.rmtree(FIG / "current")
cache = RESULTS / "cache"
cache.mkdir(exist_ok=True)
for p in list(RESULTS.glob("*.csv")):
    dst = cache / p.name
    if dst.exists():
        dst.unlink()
    shutil.move(str(p), str(dst))
if CURRENT.exists():
    shutil.rmtree(CURRENT)

print("Curated specification reporting completed")
print("Primary result table:", FINAL / "candidate_clones.csv")
print("Supplementary tables:", SUPP)
print("Final figures:", FIG)
