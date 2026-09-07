#!/usr/bin/env python3
"""Curated final reporting for the clone-fate specification.

The final presentation has two parallel analytical branches:
1. pair-specific repeated rarefaction, the primary sampling-controlled branch;
2. no-downsampling observed frequencies with Simpson metrics implemented as in
   PMID 37051229 / PMCID PMC10084475.

Only figures and tables that answer an explicit specification question are published.
Old exploratory heatmaps and duplicate summaries are not retained in the final tree.
"""
from __future__ import annotations

from pathlib import Path
import math
import os
import shutil
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FINAL = RESULTS / "final"
FIG = ROOT / "figures"

PRIMARY = 0.001
THRESHOLD_PCT = [0.05, 0.10, 0.20, 0.50]

for p in [CURRENT / "paired_T1_clone_master_nt.csv", CURRENT / "pairwise_threshold_sensitivity.csv",
          CURRENT / "simpson_raw_patient_metrics.csv", CURRENT / "simpson_raw_threshold_sensitivity.csv",
          CURRENT / "simpson_raw_primary_candidates.csv"]:
    if not p.exists():
        raise FileNotFoundError(f"Required input is missing: {p}")

master = pd.read_csv(CURRENT / "paired_T1_clone_master_nt.csv")
rare_sens = pd.read_csv(CURRENT / "pairwise_threshold_sensitivity.csv")
rare_patient = pd.read_csv(CURRENT / "pairwise_rarefaction_summary.csv")
raw_metrics = pd.read_csv(CURRENT / "simpson_raw_patient_metrics.csv")
raw_sens = pd.read_csv(CURRENT / "simpson_raw_threshold_sensitivity.csv")
raw_candidates = pd.read_csv(CURRENT / "simpson_raw_primary_candidates.csv")

for c in ["T1_fraction", "T2_fraction", "T1_fraction_rarefied_median", "T2_fraction_rarefied_median",
          "P_detection_if_unchanged_empirical"]:
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

# Harmonize primary clone annotations.
known_col = find_col(master, ["known_antigen", "antigen", "matched_antigen", "specificity", "pathogen"])
if known_col is None:
    master["known_antigen_final"] = ""
else:
    master["known_antigen_final"] = master[known_col].fillna("").astype(str)


def specificity_class(x):
    s = str(x).strip()
    low = s.lower()
    if not s or low in {"nan", "none", "unknown", "no match", "not available"}:
        return "unknown"
    if "measles" in low or "rubeola" in low:
        return "measles-associated"
    return "known non-measles memory or other specificity"

master["specificity_class"] = master["known_antigen_final"].map(specificity_class)
master["literature_query"] = master.apply(
    lambda r: f'"{r.get("cdr3aa", "")}" "{r.get("TRBV", "")}" TCR viral infection'
              + (f' "{r.get("known_antigen_final", "")}"' if str(r.get("known_antigen_final", "")).strip() else ""),
    axis=1,
)
master["literature_status"] = np.where(
    master.specificity_class.eq("unknown"),
    "requires PubMed review after candidate selection",
    "exact database annotation available; confirm in literature",
)

# Rarefied descriptive direction among primary T1-expanded clones.
primary_mask = truthy(master["expanded_T1_0_10_percent"]) if "expanded_T1_0_10_percent" in master.columns else master.T1_fraction.ge(PRIMARY)
master["rarefied_log2_T2_over_T1"] = np.where(
    master.T1_fraction_rarefied_median.gt(0) & master.T2_fraction_rarefied_median.gt(0),
    np.log2(master.T2_fraction_rarefied_median / master.T1_fraction_rarefied_median),
    np.nan,
)

def rare_direction(r):
    if not bool(primary_mask.loc[r.name]):
        return "not primary T1-expanded"
    t2 = r.get("T2_fraction_rarefied_median", np.nan)
    x = r.get("rarefied_log2_T2_over_T1", np.nan)
    if not np.isfinite(t2) or t2 <= 0:
        return "not detected"
    if x >= 1:
        return "increased"
    if x <= -1:
        return "decreased"
    return "stable"
master["rarefied_direction_0_10_percent"] = master.apply(rare_direction, axis=1)

# Merge raw branch fields back into the one-row-per-T1-clonotype master.
raw_keep_cols = [c for c in ["patient", "clone_id_nt", "raw_log2_T2_over_T1", "raw_absolute_delta",
                              "raw_fold_contraction", "raw_T1_expanded_mass_share", "raw_direction_0_10_percent",
                              "raw_strong_contraction"] if c in raw_candidates.columns]
master = master.merge(raw_candidates[raw_keep_cols], on=["patient", "clone_id_nt"], how="left")

# Final candidate class, keeping the biological fork explicit.
def candidate_class(r, branch="rarefied"):
    if branch == "rarefied":
        d = str(r.get("rarefied_direction_0_10_percent", ""))
        robust = str(r.get("sampling_robust_not_detected", "")).lower() in {"true", "1", "yes"}
        strong = str(r.get("sampling_candidate_class", "")) == "strong contraction"
    else:
        d = str(r.get("raw_direction_0_10_percent", ""))
        robust = False
        strong = bool(r.get("raw_strong_contraction", False))
    spec = str(r.get("specificity_class", "unknown"))
    if d == "increased":
        return "candidate developing or persisting measles response"
    if d == "stable":
        return "stable circulating clone"
    if d in {"decreased", "not detected"} or robust or strong:
        if spec == "measles-associated":
            return "possible contraction of measles response"
        if spec.startswith("known non-measles"):
            return "candidate immune amnesia"
        return "contracted clone with unknown specificity"
    return "other"

master["candidate_class_rarefied"] = master.apply(candidate_class, axis=1, branch="rarefied")
master["candidate_class_simpson_raw"] = master.apply(candidate_class, axis=1, branch="raw")

# Patient-aware raw-branch enrichment, analogous to the rarefied analysis.
def enrichment_for_branch(df, branch):
    z = df[primary_mask.loc[df.index]].copy() if df.index.equals(master.index) else df.copy()
    if branch == "rarefied":
        affected = z.get("sampling_candidate_class", "").astype(str).isin(["sampling-robust not detected", "strong contraction"])
        retained = z.get("sampling_candidate_class", "").astype(str).eq("retained expanded")
    else:
        affected = z.get("raw_direction_0_10_percent", "").astype(str).eq("not detected") | z.get("raw_strong_contraction", False).fillna(False).astype(bool)
        retained = z.T2_fraction.ge(PRIMARY) & ~affected
    z = z[affected | retained].copy()
    if z.empty or affected.loc[z.index].nunique() < 2:
        return pd.DataFrame()
    z["affected"] = affected.loc[z.index].astype(int)
    z["log10_T1_fraction"] = np.log10(pd.to_numeric(z.T1_fraction, errors="coerce").clip(lower=1e-8))

    rows = []
    def fit(kind, feature, indicator):
        if int(indicator.sum()) < 5 or z.loc[indicator, "patient"].nunique() < 3 or indicator.nunique() < 2:
            return
        X = pd.DataFrame({"const": 1.0, "feature": indicator.astype(float).to_numpy(),
                          "log10_T1_fraction": z.log10_T1_fraction.to_numpy(float)})
        try:
            m = sm.GEE(z.affected.to_numpy(float), X, groups=z.patient.astype(str).to_numpy(),
                       family=sm.families.Binomial()).fit(maxiter=100)
            coef = float(m.params["feature"]); se = float(m.bse["feature"]); p = float(m.pvalues["feature"])
        except Exception:
            coef = se = p = np.nan
        rows.append({"branch": branch, "feature_type": kind, "feature": feature,
                     "n_candidate_clones": int(indicator.sum()),
                     "n_patients": int(z.loc[indicator, "patient"].nunique()),
                     "log_odds_affected_vs_retained": coef, "standard_error": se,
                     "odds_ratio_affected_vs_retained": math.exp(coef) if np.isfinite(coef) else np.nan,
                     "p_value": p})

    for col, kind in [("TRBV", "TRBV"), ("TRBJ", "TRBJ")]:
        if col not in z.columns:
            continue
        for value in sorted(v for v in z[col].dropna().astype(str).unique() if v):
            fit(kind, value, z[col].astype(str).eq(value))
    motifs = {}
    for idx, aa in z.cdr3aa.astype(str).items():
        for kmer in set(aa[i:i+3] for i in range(max(0, len(aa)-2))):
            motifs.setdefault(kmer, set()).add(idx)
    for motif, indices in motifs.items():
        fit("CDR3_3mer", motif, pd.Series(z.index.isin(indices), index=z.index))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr"] = np.nan
    for kind in out.feature_type.unique():
        idx = out.index[(out.feature_type.eq(kind)) & out.p_value.notna()]
        if len(idx):
            out.loc[idx, "fdr"] = multipletests(out.loc[idx, "p_value"], method="fdr_bh")[1]
    return out.sort_values(["fdr", "p_value"], na_position="last")

rare_enrich_path = CURRENT / "lost_vs_retained_feature_enrichment.csv"
rare_enrich = pd.read_csv(rare_enrich_path) if rare_enrich_path.exists() and rare_enrich_path.stat().st_size > 1 else pd.DataFrame()
if not rare_enrich.empty:
    rare_enrich = rare_enrich.rename(columns={"log_odds_lost_vs_retained": "log_odds_affected_vs_retained"})
    rare_enrich["branch"] = "pairwise_rarefied"
raw_enrich = enrichment_for_branch(master, "simpson_no_downsampling")
enrichment = pd.concat([rare_enrich, raw_enrich], ignore_index=True, sort=False)

# Curated final tables.
if FINAL.exists():
    shutil.rmtree(FINAL)
FINAL.mkdir(parents=True, exist_ok=True)

final_cols = [c for c in [
    "patient", "clone_id_nt", "cdr3nt", "cdr3aa", "TRBV", "TRBJ", "age", "age_group",
    "T1_count", "T2_count", "T1_fraction", "T2_fraction", "fate_0_10_percent",
    "absolute_delta", "log2_T2_over_T1", "fold_contraction", "T1_expanded_mass_share",
    "T2_zero_upper_95", "P_detection_if_unchanged_empirical", "sampling_robust_not_detected",
    "T1_fraction_rarefied_median", "T2_fraction_rarefied_median", "rarefied_log2_T2_over_T1",
    "rarefied_direction_0_10_percent", "raw_log2_T2_over_T1", "raw_absolute_delta",
    "raw_fold_contraction", "raw_direction_0_10_percent", "public_private", "control_prevalence",
    "known_antigen_final", "specificity_class", "HLA_compatible", "motif_cluster", "T3_fate",
    "sequencing_source_T1", "sequencing_source_T2", "candidate_class_rarefied",
    "candidate_class_simpson_raw", "literature_query", "literature_status"
] if c in master.columns]
master[final_cols].to_csv(FINAL / "01_clone_master.csv", index=False)

candidates = master[primary_mask].copy()
cand_cols = [c for c in final_cols if c in candidates.columns]
candidates[cand_cols].sort_values(["patient", "T1_fraction"], ascending=[True, False]).to_csv(FINAL / "02_candidates.csv", index=False)

# Patient summary: rarefied metrics plus article-style raw Simpson metrics.
raw_wide = raw_metrics.pivot(index=["patient", "age_group"], columns="timepoint",
                             values=["simpson_D", "simpson_clonality", "normalized_richness", "top_gt1_mass"])
raw_wide.columns = [f"{metric}_{point}" for metric, point in raw_wide.columns]
raw_wide = raw_wide.reset_index()
patient_summary = rare_patient.merge(raw_wide, on=["patient", "age_group"], how="left")
patient_summary.to_csv(FINAL / "03_patient_summary.csv", index=False)

rare_sens2 = rare_sens[rare_sens.threshold_percent.round(6).isin(THRESHOLD_PCT)].copy()
rare_sens2["branch"] = "pairwise_rarefied"
raw_sens2 = raw_sens[raw_sens.threshold_percent.round(6).isin(THRESHOLD_PCT)].copy()
threshold_final = pd.concat([rare_sens2, raw_sens2], ignore_index=True, sort=False)
threshold_final.to_csv(FINAL / "04_threshold_sensitivity.csv", index=False)
enrichment.to_csv(FINAL / "05_feature_enrichment.csv", index=False)

# Validation in one long purposeful table.
validation_parts = []
for label, filename in [
    ("identity", "repertoire_identity_qc.csv"),
    ("T3", "extended_T3_clone_validation.csv"),
    ("GA1_003_batch", "GA1_003_cross_batch_sensitivity.csv"),
]:
    p = CURRENT / filename
    if p.exists() and p.stat().st_size > 1:
        x = pd.read_csv(p); x.insert(0, "validation_type", label); validation_parts.append(x)
validation = pd.concat(validation_parts, ignore_index=True, sort=False) if validation_parts else pd.DataFrame()
validation.to_csv(FINAL / "06_validation.csv", index=False)

# Final figure tree. PNG only to avoid duplicate clutter.
for d in ["current", "group", "main", "downsampling_pairwise", "sensitivity_thresholds", "supplementary", "rarefied", "simpson_raw", "validation", "patients"]:
    p = FIG / d
    if d not in {"current"} and p.exists():
        shutil.rmtree(p)
for d in [FIG / "patients", FIG / "rarefied", FIG / "simpson_raw", FIG / "validation"]:
    d.mkdir(parents=True, exist_ok=True)

# Keep patient panels because they are the requested Venn layer before cohort aggregation.
patient_src = FIG / "current" / "patients"
if patient_src.exists():
    for p in sorted(patient_src.glob("*.png")):
        shutil.copy2(p, FIG / "patients" / p.name)


def mean_by_age(df, value):
    return [pd.to_numeric(df.loc[df.age_group.eq(a), value], errors="coerce").mean() for a in ["child", "adult"]]


def stacked_fate(ax, sens, title):
    z = sens[np.isclose(pd.to_numeric(sens.threshold_percent, errors="coerce"), 0.10)].copy()
    xs = np.arange(2)
    bottom = np.zeros(2)
    for col, label in [("retained_mass_percent_median", "Retained expanded"),
                       ("below_mass_percent_median", "Detected below threshold"),
                       ("not_detected_mass_percent_median", "Not detected")]:
        if col not in z.columns:
            col = col.replace("_median", "")
        vals = np.array(mean_by_age(z, col), dtype=float)
        ax.bar(xs, vals, bottom=bottom, label=label)
        bottom += np.nan_to_num(vals)
    ax.set_xticks(xs, ["Children", "Adults"])
    ax.set_ylim(0, 100)
    ax.set_ylabel("Percent of T1 expanded repertoire mass")
    ax.set_title(title)
    ax.legend(fontsize=8)


def shared_scatter(ax, df, xcol, ycol, title):
    z = df[primary_mask & pd.to_numeric(df[ycol], errors="coerce").gt(0)].copy()
    x = 100 * pd.to_numeric(z[xcol], errors="coerce")
    y = 100 * pd.to_numeric(z[ycol], errors="coerce")
    ax.scatter(x, y, s=14, alpha=0.35)
    lim = max(float(np.nanmax(x)) if len(x) else 1, float(np.nanmax(y)) if len(y) else 1)
    ax.plot([0, lim], [0, lim], linestyle="--", linewidth=1)
    ax.set_xlabel("T1 repertoire fraction, percent")
    ax.set_ylabel("T2 repertoire fraction, percent")
    ax.set_title(title)
    ax.grid(alpha=0.2)

# RAREFIED 01: shared segment plus complete fate.
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), constrained_layout=True)
shared_scatter(axes[0], master, "T1_fraction_rarefied_median", "T2_fraction_rarefied_median",
               "Shared T1-expanded clonotypes change in representation")
stacked_fate(axes[1], rare_sens2, "Fate of T1 expanded repertoire after pair-specific rarefaction")
fig.suptitle("Primary sampling-controlled T1 to T2 clone fate", fontsize=14, fontweight="bold")
fig.savefig(FIG / "rarefied" / "01_shared_and_fate.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# RAREFIED 02: sampling null and threshold sensitivity.
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), constrained_layout=True)
z = master[primary_mask].copy()
pdet = pd.to_numeric(z.P_detection_if_unchanged_empirical, errors="coerce")
axes[0].scatter(100*z.T1_fraction, pdet, s=14, alpha=0.35)
rob = truthy(z["sampling_robust_not_detected"]) if "sampling_robust_not_detected" in z.columns else pd.Series(False, index=z.index)
axes[0].scatter(100*z.loc[rob, "T1_fraction"], pdet.loc[rob], marker="*", s=60, label="Sampling-robust not detected")
axes[0].set_xlabel("T1 repertoire fraction, percent"); axes[0].set_ylabel("P detection if unchanged")
axes[0].set_ylim(0, 1.02); axes[0].set_title("Internal T1 sampling null"); axes[0].legend(fontsize=8); axes[0].grid(alpha=0.2)
for age, label in [("child", "Children"), ("adult", "Adults")]:
    q = rare_sens2[rare_sens2.age_group.eq(age)]
    g = q.groupby("threshold_percent")["not_detected_mass_percent_median"].mean()
    axes[1].plot(g.index, g.values, marker="o", label=label)
axes[1].set_xlabel("Expanded threshold, percent"); axes[1].set_ylabel("Mean T1 mass not detected at T2, percent")
axes[1].set_title("Threshold sensitivity of apparent loss"); axes[1].legend(); axes[1].grid(alpha=0.2)
fig.savefig(FIG / "rarefied" / "02_sampling_and_sensitivity.png", dpi=220, bbox_inches="tight"); plt.close(fig)


def dumbbell_and_trbv(branch, xcol, ycol, dircol, outfile):
    z = master[primary_mask].copy()
    z = z[pd.to_numeric(z[ycol], errors="coerce").gt(0)].copy()
    groups = []
    for cat in ["increased", "stable", "decreased"]:
        q = z[z[dircol].astype(str).eq(cat)].copy()
        if cat == "increased": q = q.sort_values("rarefied_log2_T2_over_T1" if branch=="rarefied" else "raw_log2_T2_over_T1", ascending=False).head(5)
        elif cat == "decreased": q = q.sort_values("rarefied_log2_T2_over_T1" if branch=="rarefied" else "raw_log2_T2_over_T1", ascending=True).head(5)
        else:
            col = "rarefied_log2_T2_over_T1" if branch=="rarefied" else "raw_log2_T2_over_T1"
            q = q.assign(absfc=pd.to_numeric(q[col], errors="coerce").abs()).sort_values("absfc").head(5)
        q["display_group"] = cat
        groups.append(q)
    top = pd.concat(groups, ignore_index=True) if groups else pd.DataFrame()
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True)
    if top.empty:
        axes[0].text(0.5,0.5,"No nonzero candidates",ha="center"); axes[0].axis("off")
    else:
        labels = [f'{r.patient} {short_seq(r.cdr3aa)} {r.TRBV}' for r in top.itertuples(index=False)]
        y = np.arange(len(top))
        x1 = 100*pd.to_numeric(top[xcol], errors="coerce").to_numpy(float)
        x2 = 100*pd.to_numeric(top[ycol], errors="coerce").to_numpy(float)
        for yi,a,b in zip(y,x1,x2): axes[0].plot([a,b],[yi,yi],linewidth=1)
        axes[0].scatter(x1,y,label="T1",s=32); axes[0].scatter(x2,y,label="T2",s=32)
        axes[0].set_yticks(y,labels); axes[0].invert_yaxis(); axes[0].set_xlabel("Repertoire fraction, percent")
        axes[0].set_title("Top 5 increased, stable and decreased nonzero clonotypes")
        axes[0].legend(); axes[0].grid(axis="x",alpha=0.2)
    dec = z[z[dircol].astype(str).eq("decreased")].copy()
    if dec.empty:
        axes[1].text(0.5,0.5,"No decreased nonzero clonotypes",ha="center"); axes[1].axis("off")
    else:
        g = dec.groupby("TRBV")[xcol].sum().sort_values(ascending=False).head(10).sort_values()
        axes[1].barh(g.index,100*g.values); axes[1].set_xlabel("T1 repertoire fraction carried by decreased clones, percent")
        axes[1].set_title("TRBV carried by the most decreased nonzero clone pool"); axes[1].grid(axis="x",alpha=0.2)
    fig.suptitle("Candidate clone trajectories and TRBV context", fontsize=14, fontweight="bold")
    fig.savefig(outfile,dpi=220,bbox_inches="tight"); plt.close(fig)

dumbbell_and_trbv("rarefied", "T1_fraction_rarefied_median", "T2_fraction_rarefied_median",
                  "rarefied_direction_0_10_percent", FIG/"rarefied"/"03_candidate_dumbbell_and_trbv.png")


def enrichment_specificity(branch, outfile, affected_mask):
    fig, axes = plt.subplots(1,2,figsize=(14,6),constrained_layout=True)
    e = enrichment[enrichment.branch.astype(str).eq(branch)].copy() if not enrichment.empty else pd.DataFrame()
    e = e[pd.to_numeric(e.get("fdr", np.nan), errors="coerce").le(0.05)].copy() if not e.empty else e
    if e.empty:
        axes[0].text(0.5,0.5,"No features pass FDR 0.05",ha="center",va="center"); axes[0].axis("off")
    else:
        e = e.assign(abscoef=pd.to_numeric(e.log_odds_affected_vs_retained, errors="coerce").abs()).sort_values("abscoef",ascending=False).head(12).sort_values("log_odds_affected_vs_retained")
        y=np.arange(len(e)); coef=pd.to_numeric(e.log_odds_affected_vs_retained,errors="coerce").to_numpy(float)
        se=pd.to_numeric(e.get("standard_error",np.nan),errors="coerce").to_numpy(float)
        axes[0].errorbar(coef,y,xerr=1.96*se,fmt="o"); axes[0].axvline(0,linewidth=1)
        axes[0].set_yticks(y,[f'{a}:{b}' for a,b in zip(e.feature_type,e.feature)])
        axes[0].set_xlabel("Log odds affected versus retained"); axes[0].set_title("Patient-aware enriched features, FDR at most 0.05")
    a = master[primary_mask & affected_mask].copy()
    if a.empty:
        axes[1].text(0.5,0.5,"No affected candidates",ha="center"); axes[1].axis("off")
    else:
        counts=a.specificity_class.value_counts().reindex(["measles-associated","known non-measles memory or other specificity","unknown"],fill_value=0)
        axes[1].barh(counts.index,counts.values); axes[1].set_xlabel("Affected T1-expanded clonotypes")
        axes[1].set_title("Specificity fork for contracted or lost candidates")
    fig.savefig(outfile,dpi=220,bbox_inches="tight"); plt.close(fig)

rare_aff = master.get("sampling_candidate_class", "").astype(str).isin(["sampling-robust not detected","strong contraction"])
enrichment_specificity("pairwise_rarefied", FIG/"rarefied"/"04_enrichment_and_specificity.png", rare_aff)

# SIMPSON RAW 01: article metrics plus raw fate.
fig, axes = plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
for ax, metric, title in [
    (axes[0,0],"simpson_clonality","Simpson clonality = sqrt(sum p squared)"),
    (axes[0,1],"simpson_D","Simpson D = sum p squared"),
    (axes[1,0],"normalized_richness","Normalized richness = unique productive clones divided by productive reads"),
]:
    for patient, z in raw_metrics.groupby("patient"):
        if {"T1","T2"}.issubset(set(z.timepoint)):
            a=float(z.loc[z.timepoint.eq("T1"),metric].iloc[0]); b=float(z.loc[z.timepoint.eq("T2"),metric].iloc[0])
            ax.plot([1,2],[a,b],alpha=0.25); ax.scatter([1,2],[a,b],s=12)
    ax.set_xticks([1,2],["T1","T2"]); ax.set_title(title); ax.grid(alpha=0.2)
stacked_fate(axes[1,1], raw_sens2, "Observed fate without downsampling")
fig.suptitle("No-downsampling Simpson branch following PMID 37051229",fontsize=14,fontweight="bold")
fig.savefig(FIG/"simpson_raw"/"01_simpson_metrics_and_fate.png",dpi=220,bbox_inches="tight"); plt.close(fig)

# SIMPSON RAW 02: shared changes plus threshold sensitivity.
fig, axes = plt.subplots(1,2,figsize=(13.5,5.5),constrained_layout=True)
shared_scatter(axes[0], master, "T1_fraction", "T2_fraction", "Observed shared T1-expanded clonotypes")
for age,label in [("child","Children"),("adult","Adults")]:
    q=raw_sens2[raw_sens2.age_group.eq(age)]; g=q.groupby("threshold_percent")["not_detected_mass_percent"].mean()
    axes[1].plot(g.index,g.values,marker="o",label=label)
axes[1].set_xlabel("Expanded threshold, percent"); axes[1].set_ylabel("Mean T1 mass not detected at T2, percent")
axes[1].set_title("Threshold sensitivity without downsampling"); axes[1].legend(); axes[1].grid(alpha=0.2)
fig.savefig(FIG/"simpson_raw"/"02_shared_and_sensitivity.png",dpi=220,bbox_inches="tight"); plt.close(fig)

dumbbell_and_trbv("raw", "T1_fraction", "T2_fraction", "raw_direction_0_10_percent",
                  FIG/"simpson_raw"/"03_candidate_dumbbell_and_trbv.png")
raw_aff = master.get("raw_direction_0_10_percent", "").astype(str).eq("not detected") | master.get("raw_strong_contraction", False).fillna(False).astype(bool)
enrichment_specificity("simpson_no_downsampling", FIG/"simpson_raw"/"04_enrichment_and_specificity.png", raw_aff)

# Validation: retain only the useful combined validation figure if generated upstream.
old_val = FIG/"current"/"26_validation_identity_T3_batch.png"
if old_val.exists(): shutil.copy2(old_val, FIG/"validation"/"01_identity_T3_batch.png")

# Branch concordance is useful because it shows whether conclusions survive the alternative normalization strategy.
fig, ax = plt.subplots(figsize=(8,5.5),constrained_layout=True)
z=master[primary_mask].copy(); tab=pd.crosstab(z.rarefied_direction_0_10_percent,z.raw_direction_0_10_percent,normalize="index")*100
bottom=np.zeros(len(tab)); xs=np.arange(len(tab))
for col in tab.columns:
    vals=tab[col].to_numpy(float); ax.bar(xs,vals,bottom=bottom,label=col); bottom+=vals
ax.set_xticks(xs,tab.index,rotation=15,ha="right"); ax.set_ylim(0,100); ax.set_ylabel("Percent of clones")
ax.set_title("Concordance of clone direction between rarefied and Simpson raw branches"); ax.legend(fontsize=8)
fig.savefig(FIG/"validation"/"02_branch_concordance.png",dpi=220,bbox_inches="tight"); plt.close(fig)

# Figure index.
index = FIG/"FIGURE_INDEX.txt"
index.write_text("\n".join([
    "FINAL FIGURE SET",
    "",
    "patients: one patient panel per accepted T1-T2 pair, retaining the complete and expanded Venn layer",
    "rarefied/01_shared_and_fate.png: shared-clone abundance changes plus complete T1-expanded fate",
    "rarefied/02_sampling_and_sensitivity.png: empirical internal null plus threshold sensitivity",
    "rarefied/03_candidate_dumbbell_and_trbv.png: top increased stable decreased nonzero clonotypes and TRBV context",
    "rarefied/04_enrichment_and_specificity.png: lost-versus-retained enrichment plus measles/non-measles/unknown fork",
    "simpson_raw/01_simpson_metrics_and_fate.png: no-downsampling article-style Simpson metrics and raw fate",
    "simpson_raw/02_shared_and_sensitivity.png: raw shared-clone changes and threshold sensitivity",
    "simpson_raw/03_candidate_dumbbell_and_trbv.png: analogous raw candidate trajectories and TRBV context",
    "simpson_raw/04_enrichment_and_specificity.png: analogous raw enrichment and specificity fork",
    "validation/01_identity_T3_batch.png: identity QC, available T3 validation and GA1_003 batch sensitivity",
    "validation/02_branch_concordance.png: robustness of clone direction across the two analytical branches",
    "",
    "Primary physical clone key: CDR3 nucleotide + TRBV + TRBJ",
    "Functional/public key: CDR3 amino acid + TRBV + TRBJ",
    "Primary expanded threshold: 0.10 percent",
    "Sensitivity thresholds: 0.05, 0.10, 0.20 and 0.50 percent",
    "Pair-specific rarefaction uses at least 500 iterations at the within-pair minimum depth",
    "No-downsampling branch uses Simpson D and sqrt(Simpson D) exactly as described in PMID 37051229",
    "Actual PubMed candidate review is intentionally performed after this run using the literature_query column in results/final/02_candidates.csv",
]) + "\n", encoding="utf-8")

# Remove the temporary figure staging tree.
if (FIG/"current").exists(): shutil.rmtree(FIG/"current")

# Keep results visually clean: all root CSV caches are moved to results/cache, and
# temporary results/current is removed after curated final tables have been written.
cache = RESULTS/"cache"; cache.mkdir(exist_ok=True)
for p in list(RESULTS.glob("*.csv")):
    dst=cache/p.name
    if dst.exists(): dst.unlink()
    shutil.move(str(p),str(dst))
if CURRENT.exists(): shutil.rmtree(CURRENT)

print("Curated specification reporting completed")
print("Final tables:", FINAL)
print("Final figures:", FIG)
