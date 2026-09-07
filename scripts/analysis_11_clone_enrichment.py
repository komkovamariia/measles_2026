#!/usr/bin/env python3
"""Patient-aware enrichment among lost versus retained T1-expanded nucleotide clonotypes.

The analysis uses the primary 0.10 percent T1-expanded candidate pool created by
analysis_10_clone_fate_v2.py. TRBV, TRBJ and CDR3 amino-acid 3-mer features are tested
with a binomial GEE. Patient ID defines the correlation group and log10 T1 fraction is
included as a covariate so clone abundance is not ignored.
"""
from pathlib import Path
import math
import os
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
CURRENT = ROOT / "results" / "current"
FIG = ROOT / "figures" / "current"
CURRENT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

master_path = CURRENT / "paired_T1_clone_master_nt.csv"
if not master_path.exists():
    subprocess.run([sys.executable, str(ROOT / "scripts" / "analysis_10_clone_fate_v2.py")], check=True)
master = pd.read_csv(master_path)

z = master[
    master.expanded_T1_0_10_percent.astype(str).str.lower().isin(["true","1","yes"])
    & master.fate_0_10_percent.isin(["retained expanded","not detected"])
].copy()
if z.empty:
    raise RuntimeError("No retained/lost primary T1-expanded clones are available")
z["lost"] = z.fate_0_10_percent.eq("not detected").astype(int)
z["log10_T1_fraction"] = np.log10(pd.to_numeric(z.T1_fraction, errors="coerce").clip(lower=1e-8))


def fit_feature(feature_type, feature, indicator):
    n_clones = int(indicator.sum())
    n_patients = int(z.loc[indicator, "patient"].nunique())
    if n_clones < 5 or n_patients < 3 or indicator.nunique() < 2:
        return None
    X = pd.DataFrame({
        "const": 1.0,
        "feature": indicator.astype(float).to_numpy(),
        "log10_T1_fraction": z.log10_T1_fraction.to_numpy(float),
    })
    try:
        model = sm.GEE(
            z.lost.to_numpy(float), X,
            groups=z.patient.astype(str).to_numpy(),
            family=sm.families.Binomial(),
        ).fit(maxiter=100)
        coef = float(model.params["feature"])
        se = float(model.bse["feature"])
        p = float(model.pvalues["feature"])
    except Exception:
        coef = se = p = np.nan
    return {
        "feature_type": feature_type,
        "feature": feature,
        "n_candidate_clones": n_clones,
        "n_patients": n_patients,
        "log_odds_lost_vs_retained": coef,
        "standard_error": se,
        "odds_ratio_lost_vs_retained": math.exp(coef) if np.isfinite(coef) else np.nan,
        "p_value": p,
    }

rows = []
for col, kind in [("TRBV","TRBV"),("TRBJ","TRBJ")]:
    for value in sorted(v for v in z[col].dropna().astype(str).unique() if v):
        rec = fit_feature(kind, value, z[col].astype(str).eq(value))
        if rec is not None:
            rows.append(rec)

# Short CDR3 amino-acid motifs. These are descriptive sequence features, not antigen
# specificity assignments.
motifs = {}
for idx, aa in z.cdr3aa.astype(str).items():
    for kmer in set(aa[i:i+3] for i in range(max(0, len(aa)-2))):
        motifs.setdefault(kmer, set()).add(idx)
for motif, indices in motifs.items():
    indicator = pd.Series(z.index.isin(indices), index=z.index)
    rec = fit_feature("CDR3_3mer", motif, indicator)
    if rec is not None:
        rows.append(rec)

enrichment = pd.DataFrame(rows)
if enrichment.empty:
    enrichment.to_csv(CURRENT / "lost_vs_retained_feature_enrichment.csv", index=False)
    print("No features passed the minimum clone and patient support filters")
    raise SystemExit(0)

enrichment["fdr"] = np.nan
for kind in enrichment.feature_type.unique():
    idx = enrichment.index[(enrichment.feature_type == kind) & enrichment.p_value.notna()]
    if len(idx):
        enrichment.loc[idx,"fdr"] = multipletests(enrichment.loc[idx,"p_value"], method="fdr_bh")[1]
enrichment = enrichment.sort_values(["fdr","p_value"], na_position="last")
enrichment.to_csv(CURRENT / "lost_vs_retained_feature_enrichment.csv", index=False)

# Compact figure: strongest data-supported effects, regardless of whether FDR reaches
# significance. FDR is shown in the table and the figure title avoids claiming significance.
e = enrichment[np.isfinite(pd.to_numeric(enrichment.log_odds_lost_vs_retained, errors="coerce"))].head(20).iloc[::-1]
if not e.empty:
    fig, ax = plt.subplots(figsize=(10.5, max(6.0, 0.34*len(e)+2.3)))
    y = np.arange(len(e))
    ax.barh(y, e.log_odds_lost_vs_retained.to_numpy(float))
    ax.axvline(0, color="black", linewidth=1)
    labels = [f"{r.feature_type}: {r.feature}  FDR {r.fdr:.3g}" if pd.notna(r.fdr) else f"{r.feature_type}: {r.feature}" for r in e.itertuples(index=False)]
    ax.set_yticks(y, labels, fontsize=8)
    ax.set_xlabel("Log odds of not detected at T2 versus retained expanded")
    ax.set_title("Features associated with loss of T1 expanded clonotypes\nPatient-aware GEE adjusted for baseline T1 clone frequency")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIG / "15_lost_vs_retained_feature_enrichment.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "15_lost_vs_retained_feature_enrichment.pdf", bbox_inches="tight")
    plt.close(fig)

# Existing sequence-neighbour clusters are summarized separately from the GEE analysis.
cluster_path = CURRENT / "sequence_neighbour_cluster_summary.csv"
if cluster_path.exists():
    c = pd.read_csv(cluster_path)
    c = c[c.n_clones >= 2].sort_values(["n_patients","n_clones"], ascending=False).head(20).iloc[::-1]
    if not c.empty:
        fig, ax = plt.subplots(figsize=(10, max(5.5,0.34*len(c)+2.3)))
        ax.barh(np.arange(len(c)), 100*c.loss_or_strong_fraction.to_numpy(float))
        ax.set_yticks(np.arange(len(c)), [f"{r.motif_cluster}  {int(r.n_patients)} patients  {int(r.n_clones)} clones" for r in c.itertuples(index=False)], fontsize=8)
        ax.set_xlim(0,100)
        ax.set_xlabel("Percent of cluster clones with robust loss or strong contraction")
        ax.set_title("Sequence-neighbour clusters among T1 expanded clonotypes\nSame TRBV and CDR3 amino-acid edit distance at most 2")
        ax.grid(axis="x",alpha=0.2)
        fig.tight_layout()
        fig.savefig(FIG / "16_sequence_neighbour_clusters.png",dpi=220,bbox_inches="tight")
        fig.savefig(FIG / "16_sequence_neighbour_clusters.pdf",bbox_inches="tight")
        plt.close(fig)

print("Patient-aware lost versus retained enrichment completed")
print("Candidate clones", len(z))
print("Tested features", len(enrichment))
