#!/usr/bin/env python3
"""Patient-aware enrichment among affected versus retained T1-expanded clonotypes.

The primary candidate pool is T1 expanded at at least 0.10 percent. Affected clones are
those classified after pair-specific repeated rarefaction as either sampling-robust not
detected or strong contraction. Retained expanded clones are the reference group.
TRBV, TRBJ and CDR3 amino-acid 3-mer features are tested with binomial GEE. Patient ID
defines the correlation group and log10 T1 fraction is included as a covariate.
"""
from pathlib import Path
import math
import os
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
CURRENT = ROOT / "results" / "current"
CURRENT.mkdir(parents=True, exist_ok=True)

master_path = CURRENT / "paired_T1_clone_master_nt.csv"
if not master_path.exists():
    subprocess.run([sys.executable, str(ROOT / "scripts" / "analysis_10_clone_fate_v2.py")], check=True)
if "P_detection_if_unchanged_empirical" not in pd.read_csv(master_path, nrows=1).columns:
    subprocess.run([sys.executable, str(ROOT / "scripts" / "analysis_12_exact_rarefaction.py")], check=True)
master = pd.read_csv(master_path)

primary = master[
    master.expanded_T1_0_10_percent.astype(str).str.lower().isin(["true", "1", "yes"])
].copy()
affected_classes = {"sampling-robust not detected", "strong contraction"}
primary["analysis_class"] = np.where(
    primary.sampling_candidate_class.astype(str).isin(affected_classes),
    "affected",
    np.where(primary.sampling_candidate_class.astype(str).eq("retained expanded"), "retained", "other"),
)
z = primary[primary.analysis_class.isin(["affected", "retained"])].copy()
if z.empty or z.analysis_class.nunique() < 2:
    pd.DataFrame().to_csv(CURRENT / "lost_vs_retained_feature_enrichment.csv", index=False)
    print("No affected versus retained comparison is available")
    raise SystemExit(0)

z["affected"] = z.analysis_class.eq("affected").astype(int)
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
            z.affected.to_numpy(float),
            X,
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
        "affected_definition": "sampling-robust not detected or strong contraction",
        "reference_definition": "retained expanded",
    }

rows = []
for col, kind in [("TRBV", "TRBV"), ("TRBJ", "TRBJ")]:
    for value in sorted(v for v in z[col].dropna().astype(str).unique() if v):
        rec = fit_feature(kind, value, z[col].astype(str).eq(value))
        if rec is not None:
            rows.append(rec)

motifs = {}
for idx, aa in z.cdr3aa.astype(str).items():
    for kmer in set(aa[i:i + 3] for i in range(max(0, len(aa) - 2))):
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
        enrichment.loc[idx, "fdr"] = multipletests(enrichment.loc[idx, "p_value"], method="fdr_bh")[1]
enrichment = enrichment.sort_values(["fdr", "p_value"], na_position="last")
enrichment.to_csv(CURRENT / "lost_vs_retained_feature_enrichment.csv", index=False)

print("Patient-aware affected versus retained enrichment completed")
print("Primary expanded clones", len(primary))
print("Affected clones used", int(z.affected.sum()))
print("Retained clones used", int((z.affected == 0).sum()))
print("Tested features", len(enrichment))
