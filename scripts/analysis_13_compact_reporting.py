#!/usr/bin/env python3
"""Compact reporting for the specification-driven clone-fate analysis.

The detailed clone tables are the primary data product. Figures are deliberately few:
only plots that answer a distinct biological or technical question are generated.
The established patient panels remain available, while old exploratory plots are kept
only as supplementary context by the publisher.
"""
from __future__ import annotations

from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FIG = ROOT / "figures" / "current"
FIG.mkdir(parents=True, exist_ok=True)

PRIMARY = 0.001
KEY = ["cdr3nt", "cdr3aa", "v_gene", "j_gene"]
AGE_LABEL = {"child": "Children younger than 18 years", "adult": "Adults 18 years and older"}
FATE_ORDER = ["retained", "below", "not_detected"]
FATE_LABEL = {
    "retained": "Retained expanded",
    "below": "Detected below threshold",
    "not_detected": "Not detected at T2",
}

master_path = CURRENT / "paired_T1_clone_master_nt.csv"
master = pd.read_csv(master_path)
rare = pd.read_csv(CURRENT / "pairwise_rarefaction_summary.csv")
sens = pd.read_csv(CURRENT / "pairwise_threshold_sensitivity.csv")
meta = pd.read_csv(CURRENT / "sample_metadata_current.csv")
lib = pd.read_csv(RESULTS / "library_trb_clonotypes_nt.csv")

if master.patient.nunique() != 22 or rare.patient.nunique() != 22:
    raise RuntimeError("Compact reporting requires all 22 accepted T1-T2 pairs")
meta["analysis_include_final"] = meta.analysis_include_final.astype(str).str.lower().isin(["true", "1", "yes"])

# The empirical pair-specific sampling result is authoritative for the final candidate class.
primary_mask = master["expanded_T1_0_10_percent"].astype(str).str.lower().isin(["true", "1", "yes"])
if "candidate_class_analytic_raw_depth" not in master.columns:
    master["candidate_class_analytic_raw_depth"] = master.get("candidate_class", "")
master["candidate_class"] = np.where(
    primary_mask,
    master["sampling_candidate_class"].astype(str),
    "outside primary T1-expanded pool",
)
master.to_csv(master_path, index=False)

# Rebuild motif-cluster summaries after the empirical null has replaced the old analytic class.
candidate = master[primary_mask].copy()
if "motif_cluster" in candidate.columns:
    candidate["motif_cluster"] = candidate.motif_cluster.fillna("").astype(str)
cluster_rows = []
for cluster, z in candidate[candidate.motif_cluster.ne("")].groupby("motif_cluster"):
    antigens = sorted({x for x in z.get("known_antigen", pd.Series(dtype=str)).fillna("").astype(str) if x})
    cluster_rows.append({
        "motif_cluster": cluster,
        "n_clones": z.clone_id_nt.nunique(),
        "n_patients": z.patient.nunique(),
        "n_robust_loss": int(z.candidate_class.eq("sampling-robust not detected").sum()),
        "n_strong_contraction": int(z.candidate_class.eq("strong contraction").sum()),
        "n_retained_expanded": int(z.candidate_class.eq("retained expanded").sum()),
        "n_possible_replacement": int(z.get("T2_sequence_neighbour_present", False).astype(str).str.lower().isin(["true", "1", "yes"]).sum()) if "T2_sequence_neighbour_present" in z else 0,
        "known_antigens": "; ".join(antigens),
    })
clusters = pd.DataFrame(cluster_rows)
if not clusters.empty:
    clusters["loss_or_strong_fraction"] = (
        clusters.n_robust_loss + clusters.n_strong_contraction
    ) / clusters.n_clones
    clusters = clusters.sort_values(["n_patients", "n_clones"], ascending=False)
else:
    clusters = pd.DataFrame(columns=[
        "motif_cluster", "n_clones", "n_patients", "n_robust_loss",
        "n_strong_contraction", "n_retained_expanded", "n_possible_replacement",
        "known_antigens", "loss_or_strong_fraction",
    ])
clusters.to_csv(CURRENT / "sequence_neighbour_cluster_summary.csv", index=False)

# Recompute GA1_003 sensitivity from the pair-specific rarefaction, not from raw overlap.
primary_sens = sens[np.isclose(pd.to_numeric(sens.threshold_fraction, errors="coerce"), PRIMARY)].copy()
ga_rows = []
for age in ["child", "adult"]:
    za = primary_sens[primary_sens.age_group.eq(age)]
    for label, z in [
        ("all paired patients", za),
        ("excluding GA1_003", za[~za.patient.astype(str).eq("ГА1_003")]),
    ]:
        ga_rows.append({
            "age_group": age,
            "analysis_set": label,
            "n_patients": z.patient.nunique(),
            "mean_lost_T1_mass_percent": pd.to_numeric(z.not_detected_mass_percent_median, errors="coerce").mean(),
            "median_lost_T1_mass_percent": pd.to_numeric(z.not_detected_mass_percent_median, errors="coerce").median(),
            "mean_retained_T1_mass_percent": pd.to_numeric(z.retained_mass_percent_median, errors="coerce").mean(),
            "mean_expanded_jaccard": pd.to_numeric(z.expanded_jaccard_median, errors="coerce").mean(),
        })
ga = pd.DataFrame(ga_rows)
ga.to_csv(CURRENT / "GA1_003_cross_batch_sensitivity.csv", index=False)

# Extended T3 validation uses every valid longitudinal relation requested in the specification.
def point_table(patient, point):
    libs = meta.loc[
        meta.analysis_include_final
        & meta.patient_id.astype(str).eq(patient)
        & meta.clinical_timepoint.astype(str).eq(point),
        "library_id",
    ].astype(str).tolist()
    if not libs:
        return pd.DataFrame(columns=KEY + ["fraction"])
    z = lib[lib.library_id.astype(str).isin(libs)].copy()
    if z.empty:
        return pd.DataFrame(columns=KEY + ["fraction"])
    if "library_fraction" not in z.columns:
        totals = z.groupby("library_id")["count"].transform("sum")
        z["library_fraction"] = pd.to_numeric(z["count"], errors="coerce").fillna(0) / totals
    a = z.groupby(KEY, as_index=False, dropna=False).library_fraction.sum()
    a["fraction"] = a.library_fraction / len(libs)
    s = float(a.fraction.sum())
    if s > 0:
        a["fraction"] /= s
    return a[KEY + ["fraction"]]

validation_specs = [
    ("ЭГ1_042", "T1", "T3", "T1 to T3"),
    ("СА1_013", "T2", "T3", "T2 to T3"),
    ("ММ1_033", "T1", "T3", "T1 to T3"),
]
extended_rows = []
for patient, source_point, follow_point, comparison in validation_specs:
    source = point_table(patient, source_point)
    follow = point_table(patient, follow_point)
    if source.empty or follow.empty:
        continue
    fi = follow.set_index(KEY, drop=False)
    t2 = point_table(patient, "T2") if patient == "ЭГ1_042" and source_point == "T1" else pd.DataFrame()
    t2i = t2.set_index(KEY, drop=False) if not t2.empty else None
    for r in source[source.fraction.ge(PRIMARY)].itertuples(index=False):
        key = (str(r.cdr3nt), str(r.cdr3aa), str(r.v_gene), str(r.j_gene))
        f_follow = float(fi.loc[key].fraction) if key in fi.index else 0.0
        if isinstance(f_follow, pd.Series):
            f_follow = float(f_follow.iloc[0])
        fate = "retained expanded" if f_follow >= PRIMARY else ("detected below threshold" if f_follow > 0 else "not detected")
        if t2i is not None:
            f_t2 = float(t2i.loc[key].fraction) if key in t2i.index else 0.0
            if isinstance(f_t2, pd.Series):
                f_t2 = float(f_t2.iloc[0])
        else:
            f_t2 = np.nan
        extended_rows.append({
            "patient": patient,
            "comparison": comparison,
            "source_point": source_point,
            "followup_point": follow_point,
            "clone_id_nt": "|".join([str(r.cdr3nt), str(r.v_gene), str(r.j_gene)]),
            "cdr3aa": str(r.cdr3aa),
            "TRBV": str(r.v_gene),
            "TRBJ": str(r.j_gene),
            "source_fraction": float(r.fraction),
            "T2_fraction_for_T1_T2_T3_case": f_t2,
            "followup_fraction": f_follow,
            "followup_fate": fate,
            "lost_at_T2_then_reappeared_T3": bool(patient == "ЭГ1_042" and f_t2 == 0 and f_follow > 0),
        })
extended_t3 = pd.DataFrame(extended_rows)
extended_t3.to_csv(CURRENT / "extended_T3_clone_validation.csv", index=False)
if not extended_t3.empty:
    t3_summary = extended_t3.groupby(["patient", "comparison"], as_index=False).agg(
        n_source_expanded=("clone_id_nt", "size"),
        n_followup_retained_expanded=("followup_fate", lambda x: int((x == "retained expanded").sum())),
        n_followup_below_threshold=("followup_fate", lambda x: int((x == "detected below threshold").sum())),
        n_followup_not_detected=("followup_fate", lambda x: int((x == "not detected").sum())),
        n_lost_T2_reappeared_T3=("lost_at_T2_then_reappeared_T3", "sum"),
    )
else:
    t3_summary = pd.DataFrame(columns=[
        "patient", "comparison", "n_source_expanded", "n_followup_retained_expanded",
        "n_followup_below_threshold", "n_followup_not_detected", "n_lost_T2_reappeared_T3",
    ])
t3_summary.to_csv(CURRENT / "extended_T3_validation_summary.csv", index=False)

# Human-readable compact view of the full clone master table.
candidate = master[primary_mask].copy()
candidate["candidate_priority"] = (
    100 * pd.to_numeric(candidate.T1_fraction, errors="coerce").fillna(0)
    * np.where(candidate.candidate_class.eq("sampling-robust not detected"), 3.0,
      np.where(candidate.candidate_class.eq("strong contraction"), 2.0, 1.0))
)
keep_cols = [
    "patient", "clone_id_nt", "clone_id_aa", "cdr3aa", "TRBV", "TRBJ",
    "T1_raw_count_sum", "T2_raw_count_sum", "T1_percent", "T2_percent",
    "fate_0_10_percent", "absolute_delta_percentage_points", "log2_T2_over_T1",
    "T1_expanded_mass_share", "fold_contraction_T1_over_T2",
    "fold_contraction_lower_bound_if_T2_zero", "T2_zero_upper_95_pair_depth",
    "P_detection_if_unchanged", "P_detection_if_unchanged_analytic_raw_depth",
    "P_detection_in_observed_T2_rarefactions", "pair_rarefaction_depth",
    "rarefaction_iterations", "sampling_robust_not_detected", "candidate_class",
    "public_private", "control_prevalence", "known_antigen", "HLA_compatible",
    "motif_cluster", "T2_sequence_neighbour_present", "T3_fate",
    "T1_batches", "T2_batches", "cross_batch_pair", "candidate_priority",
]
keep_cols = [c for c in keep_cols if c in candidate.columns]
candidate[keep_cols].sort_values("candidate_priority", ascending=False).to_csv(
    CURRENT / "candidate_clones_compact.csv", index=False
)
master.to_csv(master_path, index=False)


def patient_dots_with_mean(ax, data, value_col, category_x, age, offset):
    vals = pd.to_numeric(data.loc[data.age_group.eq(age), value_col], errors="coerce").dropna().to_numpy(float)
    if len(vals) == 0:
        return
    jitter = np.linspace(-0.035, 0.035, len(vals)) if len(vals) > 1 else np.array([0.0])
    ax.scatter(np.full(len(vals), category_x + offset) + jitter, vals, s=24, alpha=0.75)
    ax.hlines(np.mean(vals), category_x + offset - 0.11, category_x + offset + 0.11, linewidth=3)


# Figure 21. Primary fate after pair-specific repeated rarefaction.
primary = primary_sens.copy()
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

# Figure 22. Sampling robustness of apparent clone loss.
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
axs[0].set_ylabel("T2 fraction or 95 percent detection upper bound, percent")
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
axs[1].set_ylabel("Empirical P detection if unchanged")
axs[1].set_title(f"Internal sampling null: {len(robust)} robust losses among {len(lost)} not detected clones")
axs[1].grid(alpha=0.2)
axs[1].legend(fontsize=8)
fig.suptitle("Apparent loss separated from finite-depth sampling", fontsize=14, fontweight="bold")
fig.savefig(FIG / "22_sampling_robust_clone_loss.png", dpi=220, bbox_inches="tight")
fig.savefig(FIG / "22_sampling_robust_clone_loss.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 23. Threshold sensitivity from the same pair-specific rarefaction framework.
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

# Figure 24 only when patient-aware enrichment survives FDR.
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
        ax.set_xlabel("Log odds of robust loss or strong contraction versus retained expanded")
        ax.set_title("Features enriched among affected T1 expanded clonotypes\nPatient-aware GEE adjusted for T1 clone frequency, FDR at most 0.05")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        fig.savefig(FIG / "24_lost_vs_retained_enrichment.png", dpi=220, bbox_inches="tight")
        fig.savefig(FIG / "24_lost_vs_retained_enrichment.pdf", bbox_inches="tight")
        plt.close(fig)

# Figure 25 only for a real cross-patient sequence-neighbour signal.
if not clusters.empty:
    signal = clusters[(clusters.n_patients >= 2) & ((clusters.n_robust_loss + clusters.n_strong_contraction) >= 2)].copy()
    if not signal.empty:
        signal = signal.head(12).iloc[::-1]
        fig, ax = plt.subplots(figsize=(10.5, max(5.5, 0.42 * len(signal) + 2.2)))
        ax.barh(np.arange(len(signal)), 100 * signal.loss_or_strong_fraction)
        labels = []
        for r in signal.itertuples(index=False):
            antigen = f", antigen {r.known_antigens}" if str(r.known_antigens) not in {"", "nan"} else ""
            labels.append(f"{r.motif_cluster}: {int(r.n_patients)} patients, {int(r.n_clones)} clones{antigen}")
        ax.set_yticks(np.arange(len(signal)), labels, fontsize=8)
        ax.set_xlim(0, 100)
        ax.set_xlabel("Percent of cluster clones with robust loss or strong contraction")
        ax.set_title("Cross-patient sequence neighbourhoods affected after measles\nSame TRBV and nearby CDR3 amino-acid sequence")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        fig.savefig(FIG / "25_sequence_neighbourhoods.png", dpi=220, bbox_inches="tight")
        fig.savefig(FIG / "25_sequence_neighbourhoods.pdf", bbox_inches="tight")
        plt.close(fig)

# Figure 26. One combined validation figure: identity, T3 and the GA1_003 cross-batch check.
identity = pd.read_csv(CURRENT / "repertoire_identity_qc.csv")
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

if not t3_summary.empty:
    y = np.arange(len(t3_summary))
    retained = pd.to_numeric(t3_summary.n_followup_retained_expanded, errors="coerce").fillna(0)
    below = pd.to_numeric(t3_summary.n_followup_below_threshold, errors="coerce").fillna(0)
    absent = pd.to_numeric(t3_summary.n_followup_not_detected, errors="coerce").fillna(0)
    axs[1].bar(y, retained, label="Retained expanded")
    axs[1].bar(y, below, bottom=retained, label="Detected below threshold")
    axs[1].bar(y, absent, bottom=retained + below, label="Not detected")
    axs[1].set_xticks(y, [f"{r.patient}\n{r.comparison}" for r in t3_summary.itertuples(index=False)], rotation=20, ha="right", fontsize=7)
    axs[1].set_ylabel("Expanded clones at source point")
    axs[1].set_title("Available T3 follow-up")
    axs[1].legend(fontsize=7)
else:
    axs[1].axis("off")
    axs[1].text(0.5, 0.5, "No applicable T3 validation", ha="center", va="center")

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
axs[2].set_title("GA1_003 cross-batch sensitivity after rarefaction")
axs[2].grid(axis="y", alpha=0.2)
fig.suptitle("Validation of clone-loss interpretation", fontsize=14, fontweight="bold")
fig.savefig(FIG / "26_validation_identity_T3_batch.png", dpi=220, bbox_inches="tight")
fig.savefig(FIG / "26_validation_identity_T3_batch.pdf", bbox_inches="tight")
plt.close(fig)

(CURRENT / "compact_reporting_README.txt").write_text(
    "\n".join([
        "Final reporting is intentionally compact",
        "Primary figure: fate of T1 expanded clonotypes after pair-specific repeated rarefaction",
        "Sampling figure: empirical P_detection_if_unchanged from the internal T1 null",
        "Sensitivity figure: 0.05, 0.10, 0.20, 0.50, 0.75 and 1.00 percent thresholds",
        "Enrichment figure appears only if at least one patient-aware feature passes FDR 0.05",
        "Sequence-neighbourhood figure appears only for a cross-patient affected cluster",
        "Validation combines repertoire identity, all usable T3 relations and GA1_003 cross-batch sensitivity",
        "Detailed results remain in paired_T1_clone_master_nt.csv and candidate_clones_compact.csv",
    ]) + "\n",
    encoding="utf-8",
)

print("Compact reporting completed")
print("Candidate clones", len(candidate))
print("Core figures: 21, 22, 23 and 26")
print("Figures 24 and 25 are conditional on data-supported results")
