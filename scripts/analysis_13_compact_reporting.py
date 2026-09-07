#!/usr/bin/env python3
"""Compact reporting for the nucleotide clone-fate analysis.

Only figures that answer a distinct biological or technical question are generated.
Old analyses remain available in the repository, but the final publisher exposes only
this compact set plus patient panels and a small supplementary set.
"""
from __future__ import annotations

from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
CURRENT = ROOT / "results" / "current"
FIG = ROOT / "figures" / "current"
FIG.mkdir(parents=True, exist_ok=True)

PRIMARY = 0.001
AGE_LABEL = {"child": "Children younger than 18 years", "adult": "Adults 18 years and older"}
FATE_ORDER = ["retained", "below", "not_detected"]
FATE_LABEL = {
    "retained": "Retained expanded",
    "below": "Detected below threshold",
    "not_detected": "Not detected at T2",
}

master = pd.read_csv(CURRENT / "paired_T1_clone_master_nt.csv")
rare = pd.read_csv(CURRENT / "pairwise_rarefaction_summary.csv")
sens = pd.read_csv(CURRENT / "pairwise_threshold_sensitivity.csv")

if master.patient.nunique() != 22 or rare.patient.nunique() != 22:
    raise RuntimeError("Compact reporting requires all 22 accepted T1-T2 pairs")

# A concise candidate table is the human-readable view of the full master table.
candidate = master[
    master["expanded_T1_0_10_percent"].astype(str).str.lower().isin(["true", "1", "yes"])
].copy()
candidate["candidate_priority"] = (
    100 * pd.to_numeric(candidate.T1_fraction, errors="coerce").fillna(0)
    * np.where(candidate.sampling_candidate_class.eq("sampling-robust not detected"), 3.0,
      np.where(candidate.sampling_candidate_class.eq("strong contraction"), 2.0, 1.0))
)
keep_cols = [
    "patient", "clone_id_nt", "clone_id_aa", "cdr3aa", "TRBV", "TRBJ",
    "T1_percent", "T2_percent", "fate_0_10_percent",
    "fold_contraction_T1_over_T2", "fold_contraction_lower_bound_if_T2_zero",
    "P_detection_if_unchanged", "P_detection_if_unchanged_analytic_raw_depth",
    "P_detection_in_observed_T2_rarefactions", "pair_rarefaction_depth",
    "sampling_robust_not_detected", "sampling_candidate_class",
    "public_private", "control_prevalence", "known_antigen", "HLA_compatible",
    "motif_cluster", "T2_sequence_neighbour_present", "T3_fate",
    "T1_batches", "T2_batches", "cross_batch_pair", "candidate_priority",
]
keep_cols = [c for c in keep_cols if c in candidate.columns]
candidate[keep_cols].sort_values("candidate_priority", ascending=False).to_csv(
    CURRENT / "candidate_clones_compact.csv", index=False
)


def patient_dots_with_mean(ax, data, value_col, category_x, age, offset):
    vals = pd.to_numeric(data.loc[data.age_group.eq(age), value_col], errors="coerce").dropna().to_numpy(float)
    if len(vals) == 0:
        return
    jitter = np.linspace(-0.035, 0.035, len(vals)) if len(vals) > 1 else np.array([0.0])
    ax.scatter(np.full(len(vals), category_x + offset) + jitter, vals, s=24, alpha=0.75)
    ax.hlines(np.mean(vals), category_x + offset - 0.11, category_x + offset + 0.11, linewidth=3)


# Figure 1. Primary fate after pair-specific repeated rarefaction.
primary = sens[np.isclose(sens.threshold_fraction, PRIMARY)].copy()
fig, axs = plt.subplots(1, 2, figsize=(14, 5.8), constrained_layout=True)
for ax, suffix, ylabel, title in [
    (axs[0], "count_percent_median", "Percent of T1 expanded clonotypes", "Fate by clone count"),
    (axs[1], "mass_percent_median", "Percent of T1 expanded repertoire mass", "Fate by T1 repertoire mass"),
]:
    for x, fate in enumerate(FATE_ORDER):
        col = f"{fate}_{suffix}"
        patient_dots_with_mean(ax, primary, col, x, "child", -0.13)
        patient_dots_with_mean(ax, primary, col, x, "adult", 0.13)
    ax.set_xticks(range(3), [FATE_LABEL[x] for x in FATE_ORDER], rotation=12, ha="right")
    ax.set_ylim(0, 100)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.2)
axs[0].scatter([], [], label="Children")
axs[0].scatter([], [], label="Adults")
axs[0].legend(fontsize=8)
fig.suptitle(
    "Fate of T1 expanded clonotypes after pair-specific repeated rarefaction\n"
    "Expanded threshold 0.10 percent, each patient is one biological replicate",
    fontsize=14, fontweight="bold",
)
fig.savefig(FIG / "21_primary_fate_pairwise_rarefied.png", dpi=220, bbox_inches="tight")
fig.savefig(FIG / "21_primary_fate_pairwise_rarefied.pdf", bbox_inches="tight")
plt.close(fig)


# Figure 2. Sampling robustness of apparent clone loss.
expanded = candidate.copy()
expanded["T1_percent"] = pd.to_numeric(expanded.T1_percent, errors="coerce")
expanded["T2_percent"] = pd.to_numeric(expanded.T2_percent, errors="coerce")
expanded["P_detection_if_unchanged"] = pd.to_numeric(expanded.P_detection_if_unchanged, errors="coerce")
expanded["T2_zero_upper_95_pair_depth"] = pd.to_numeric(expanded.get("T2_zero_upper_95_pair_depth", np.nan), errors="coerce")

fig, axs = plt.subplots(1, 2, figsize=(14, 5.8), constrained_layout=True)
ret = expanded[expanded.fate_0_10_percent.eq("retained expanded")]
below = expanded[expanded.fate_0_10_percent.eq("detected below 0.10 percent")]
lost = expanded[expanded.fate_0_10_percent.eq("not detected")]
robust = lost[lost.sampling_robust_not_detected.astype(str).str.lower().isin(["true", "1", "yes"])]
nonrobust = lost[~lost.index.isin(robust.index)]

for z, marker, label in [
    (ret, "o", "Retained expanded"),
    (below, "s", "Detected below 0.10 percent"),
    (nonrobust, "x", "Not detected, sampling uncertain"),
    (robust, "*", "Not detected, sampling robust"),
]:
    if z.empty:
        continue
    y = z.T2_percent.to_numpy(float).copy()
    zero = y <= 0
    if zero.any():
        upper = 100 * pd.to_numeric(z.loc[zero, "T2_zero_upper_95_pair_depth"], errors="coerce").to_numpy(float)
        y[zero] = np.where(np.isfinite(upper) & (upper > 0), upper, 1e-6)
    axs[0].scatter(z.T1_percent, y, marker=marker, s=34 if marker != "*" else 70, alpha=0.75, label=label)
axs[0].plot([0.05, 20], [0.05, 20], linestyle="--", linewidth=1)
axs[0].set_xscale("log")
axs[0].set_yscale("log")
axs[0].set_xlabel("T1 repertoire fraction, percent")
axs[0].set_ylabel("T2 repertoire fraction or 95 percent detection upper bound, percent")
axs[0].set_title("Observed contraction of T1 expanded clonotypes")
axs[0].grid(alpha=0.2)
axs[0].legend(fontsize=7)

if not lost.empty:
    axs[1].scatter(lost.T1_percent, lost.P_detection_if_unchanged, alpha=0.7, label="Not detected at T2")
if not robust.empty:
    axs[1].scatter(robust.T1_percent, robust.P_detection_if_unchanged, marker="*", s=80, label="Sampling robust")
axs[1].axhline(0.99, linestyle="--", linewidth=1)
axs[1].set_xscale("log")
axs[1].set_ylim(0, 1.02)
axs[1].set_xlabel("T1 repertoire fraction, percent")
axs[1].set_ylabel("P detection if unchanged from repeated internal null")
axs[1].set_title(f"Internal sampling null: {len(robust)} robust losses among {len(lost)} not detected clones")
axs[1].grid(alpha=0.2)
axs[1].legend(fontsize=8)
fig.suptitle("Apparent loss separated from finite-depth sampling", fontsize=14, fontweight="bold")
fig.savefig(FIG / "22_sampling_robust_clone_loss.png", dpi=220, bbox_inches="tight")
fig.savefig(FIG / "22_sampling_robust_clone_loss.pdf", bbox_inches="tight")
plt.close(fig)


# Figure 3. Threshold sensitivity, based on the same pair-specific rarefaction framework.
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True, sharey=True)
for ax, age in zip(axs, ["child", "adult"]):
    z = sens[sens.age_group.eq(age)].copy()
    g = z.groupby("threshold_percent", as_index=False).agg(
        retained=("retained_mass_percent_median", "mean"),
        below=("below_mass_percent_median", "mean"),
        not_detected=("not_detected_mass_percent_median", "mean"),
    )
    ax.plot(g.threshold_percent, g.retained, marker="o", label="Retained expanded")
    ax.plot(g.threshold_percent, g.below, marker="o", label="Detected below threshold")
    ax.plot(g.threshold_percent, g.not_detected, marker="o", label="Not detected")
    ax.set_xlabel("Expansion threshold, percent")
    ax.set_title(AGE_LABEL[age])
    ax.grid(alpha=0.2)
axs[0].set_ylabel("Mean percent of T1 expanded repertoire mass")
axs[0].legend(fontsize=8)
fig.suptitle(
    "Sensitivity of T1 clone fate to the expansion threshold\n"
    "Pair-specific repeated rarefaction at 0.05, 0.10, 0.20, 0.50, 0.75 and 1.00 percent",
    fontsize=14, fontweight="bold",
)
fig.savefig(FIG / "23_threshold_sensitivity_pairwise_rarefied.png", dpi=220, bbox_inches="tight")
fig.savefig(FIG / "23_threshold_sensitivity_pairwise_rarefied.pdf", bbox_inches="tight")
plt.close(fig)


# Figure 4 is generated only when the patient-aware model finds an FDR-supported feature.
enrichment_path = CURRENT / "lost_vs_retained_feature_enrichment.csv"
if enrichment_path.exists() and enrichment_path.stat().st_size > 0:
    e = pd.read_csv(enrichment_path)
    e["fdr"] = pd.to_numeric(e.get("fdr"), errors="coerce")
    e["coef"] = pd.to_numeric(e.get("log_odds_lost_vs_retained"), errors="coerce")
    e["se"] = pd.to_numeric(e.get("standard_error"), errors="coerce")
    sig = e[e.fdr.le(0.05) & e.coef.notna() & e.se.notna()].copy()
    if not sig.empty:
        sig = sig.sort_values(["fdr", "p_value"]).head(14).iloc[::-1]
        y = np.arange(len(sig))
        fig, ax = plt.subplots(figsize=(10.5, max(5.5, 0.38 * len(sig) + 2.2)))
        ax.errorbar(sig.coef, y, xerr=1.96 * sig.se, fmt="o", capsize=2)
        ax.axvline(0, linewidth=1)
        ax.set_yticks(y, [f"{r.feature_type}: {r.feature}" for r in sig.itertuples(index=False)], fontsize=8)
        ax.set_xlabel("Log odds of not detected at T2 versus retained expanded")
        ax.set_title("Features enriched among lost T1 expanded clonotypes\nPatient-aware GEE adjusted for T1 clone frequency, FDR at most 0.05")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        fig.savefig(FIG / "24_lost_vs_retained_enrichment.png", dpi=220, bbox_inches="tight")
        fig.savefig(FIG / "24_lost_vs_retained_enrichment.pdf", bbox_inches="tight")
        plt.close(fig)


# Figure 5 is generated only for a real cross-patient sequence-neighbour signal.
cluster_path = CURRENT / "sequence_neighbour_cluster_summary.csv"
if cluster_path.exists() and cluster_path.stat().st_size > 0:
    c = pd.read_csv(cluster_path)
    for col in ["n_patients", "n_clones", "n_robust_loss", "n_strong_contraction", "loss_or_strong_fraction"]:
        c[col] = pd.to_numeric(c[col], errors="coerce")
    c["n_loss_or_strong"] = c.n_robust_loss.fillna(0) + c.n_strong_contraction.fillna(0)
    signal = c[(c.n_patients >= 2) & (c.n_loss_or_strong >= 2)].copy()
    if not signal.empty:
        signal = signal.sort_values(["n_patients", "n_loss_or_strong", "n_clones"], ascending=False).head(12).iloc[::-1]
        fig, ax = plt.subplots(figsize=(10.5, max(5.5, 0.38 * len(signal) + 2.2)))
        ax.barh(np.arange(len(signal)), 100 * signal.loss_or_strong_fraction)
        ax.set_yticks(
            np.arange(len(signal)),
            [f"{r.motif_cluster}: {int(r.n_patients)} patients, {int(r.n_clones)} clones" for r in signal.itertuples(index=False)],
            fontsize=8,
        )
        ax.set_xlim(0, 100)
        ax.set_xlabel("Percent of cluster clones with robust loss or strong contraction")
        ax.set_title("Cross-patient sequence neighbourhoods affected after measles\nSame TRBV and nearby CDR3 amino-acid sequence")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        fig.savefig(FIG / "25_sequence_neighbourhoods.png", dpi=220, bbox_inches="tight")
        fig.savefig(FIG / "25_sequence_neighbourhoods.pdf", bbox_inches="tight")
        plt.close(fig)


# Figure 6. Validation combines sample identity, T3 follow-up and GA1_003 batch sensitivity.
identity = pd.read_csv(CURRENT / "repertoire_identity_qc.csv")
t3_path = CURRENT / "T3_clone_fate_validation.csv"
ga_path = CURRENT / "GA1_003_cross_batch_sensitivity.csv"
t3 = pd.read_csv(t3_path) if t3_path.exists() and t3_path.stat().st_size > 0 else pd.DataFrame()
ga = pd.read_csv(ga_path) if ga_path.exists() and ga_path.stat().st_size > 0 else pd.DataFrame()

fig, axs = plt.subplots(1, 3, figsize=(16, 5.4), constrained_layout=True)
identity = identity.sort_values("own_T1_T2_morisita")
x = np.arange(len(identity))
axs[0].scatter(x, pd.to_numeric(identity.own_T1_T2_morisita, errors="coerce"))
failed = ~identity.identity_pass_rank1.astype(str).str.lower().isin(["true", "1", "yes"])
if failed.any():
    axs[0].scatter(x[failed], pd.to_numeric(identity.loc[failed, "own_T1_T2_morisita"], errors="coerce"), marker="x", s=60, label="Own T2 not rank 1")
axs[0].set_xticks(x, identity.patient.astype(str), rotation=90, fontsize=6)
axs[0].set_ylim(0, 1)
axs[0].set_ylabel("Morisita-Horn similarity")
axs[0].set_title("Identity QC: own T1 to own T2")
axs[0].grid(axis="y", alpha=0.2)
if failed.any():
    axs[0].legend(fontsize=7)

if not t3.empty:
    y = np.arange(len(t3))
    absent = pd.to_numeric(t3.n_T2_lost_still_absent_T3, errors="coerce").fillna(0)
    reap = pd.to_numeric(t3.n_T2_lost_reappeared_T3, errors="coerce").fillna(0)
    axs[1].bar(y, absent, label="Still absent at T3")
    axs[1].bar(y, reap, bottom=absent, label="Reappeared at T3")
    axs[1].set_xticks(y, t3.patient.astype(str), rotation=30, ha="right")
    axs[1].set_ylabel("T1 expanded clones not detected at T2")
    axs[1].set_title("T3 validation of apparent loss")
    axs[1].legend(fontsize=7)
else:
    axs[1].axis("off")
    axs[1].text(0.5, 0.5, "No applicable T3 validation", ha="center", va="center")

if not ga.empty and {"age_group", "analysis_set", "mean_lost_T1_mass_percent"}.issubset(ga.columns):
    groups = []
    values = []
    for age in ["child", "adult"]:
        for analysis_set in ["all paired patients", "excluding GA1_003"]:
            q = ga[(ga.age_group.eq(age)) & (ga.analysis_set.eq(analysis_set))]
            if not q.empty:
                groups.append(f"{age}\n{analysis_set.replace(' paired patients', '')}")
                values.append(float(q.mean_lost_T1_mass_percent.iloc[0]))
    axs[2].bar(np.arange(len(values)), values)
    axs[2].set_xticks(np.arange(len(values)), groups, rotation=20, ha="right", fontsize=7)
    axs[2].set_ylabel("Mean T1 expanded mass not detected at T2, percent")
    axs[2].set_title("GA1_003 cross-batch sensitivity")
    axs[2].grid(axis="y", alpha=0.2)
else:
    axs[2].axis("off")
    axs[2].text(0.5, 0.5, "GA1_003 sensitivity table unavailable", ha="center", va="center")

fig.suptitle("Validation of clone-loss interpretation", fontsize=14, fontweight="bold")
fig.savefig(FIG / "26_validation_identity_T3_batch.png", dpi=220, bbox_inches="tight")
fig.savefig(FIG / "26_validation_identity_T3_batch.pdf", bbox_inches="tight")
plt.close(fig)

print("Compact reporting completed")
print("Primary fate figure generated")
print("Sampling robustness figure generated")
print("Threshold sensitivity figure generated")
print("Enrichment and sequence-neighbour figures are emitted only when supported by the data")
print("Validation figure generated")
