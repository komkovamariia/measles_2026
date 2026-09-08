#!/usr/bin/env python3
"""No-downsampling branch using Simpson repertoire metrics from PMCID PMC10084475.

This branch is complementary to the pair-specific rarefaction analysis. It never
subsamples reads. Simpson's D is sum(p_i^2). Simpson clonality is sqrt(Simpson's D),
matching the Methods definition in PMID 37051229. Normalized richness is the number
of unique productive nucleotide clonotypes divided by total productive TRB reads.

The same T1 clone-fate specification is applied to observed normalized frequencies
at 0.05, 0.10, then 0.20 to 1.50 percent in 0.10 percentage-point steps.
The primary threshold remains 0.10 percent.
"""
from __future__ import annotations

from pathlib import Path
import os
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
CURRENT.mkdir(parents=True, exist_ok=True)

PRIMARY = 0.001
THRESHOLD_PCT = [0.05, 0.10] + [round(x / 10, 2) for x in range(2, 16)]
THRESHOLDS = [x / 100.0 for x in THRESHOLD_PCT]
KEY = ["cdr3nt", "cdr3aa", "v_gene", "j_gene"]

meta = pd.read_csv(CURRENT / "sample_metadata_current.csv")
master = pd.read_csv(CURRENT / "paired_T1_clone_master_nt.csv")
lib = pd.read_csv(RESULTS / "library_trb_clonotypes_nt.csv")

meta["analysis_include_final"] = meta["analysis_include_final"].astype(str).str.lower().isin(["true", "1", "yes"])
lib["count"] = pd.to_numeric(lib["count"], errors="coerce").fillna(0.0)

case = meta[meta.analysis_include_final & meta.disease_status.eq("measles")].copy()
paired = []
for patient, z in case.groupby("patient_id"):
    if {"T1", "T2"}.issubset(set(z.clinical_timepoint.astype(str))):
        paired.append(str(patient))
paired = sorted(paired)
if len(paired) != 22:
    raise RuntimeError(f"Expected 22 paired T1-T2 patients, found {len(paired)}")


def lib_metrics(z: pd.DataFrame) -> dict:
    counts = pd.to_numeric(z["count"], errors="coerce").fillna(0.0).to_numpy(float)
    counts = counts[counts > 0]
    total = float(counts.sum())
    if total <= 0 or len(counts) == 0:
        return {"productive_reads": 0.0, "richness": 0.0, "simpson_D": np.nan,
                "simpson_clonality": np.nan, "normalized_richness": np.nan,
                "top_gt1_mass": np.nan}
    p = counts / total
    simpson_d = float(np.sum(p * p))
    return {
        "productive_reads": total,
        "richness": float(len(counts)),
        "simpson_D": simpson_d,
        "simpson_clonality": float(np.sqrt(simpson_d)),
        "normalized_richness": float(len(counts) / total),
        "top_gt1_mass": float(p[p > 0.01].sum()),
    }


metric_rows = []
for patient in paired:
    pmeta = case[case.patient_id.astype(str).eq(patient) & case.clinical_timepoint.isin(["T1", "T2"])].copy()
    for point in ["T1", "T2"]:
        lids = pmeta.loc[pmeta.clinical_timepoint.eq(point), "library_id"].astype(str).tolist()
        vals = []
        for lid in lids:
            z = lib[lib.library_id.astype(str).eq(lid)].copy()
            if z.empty:
                raise RuntimeError(f"Missing nucleotide clonotypes for {lid}")
            vals.append(lib_metrics(z))
        row = {
            "patient": patient,
            "age_group": str(pmeta.age_group.iloc[0]),
            "timepoint": point,
            "libraries": ";".join(lids),
            "n_libraries": len(lids),
        }
        for key in vals[0]:
            row[key] = float(np.nanmean([v[key] for v in vals]))
        metric_rows.append(row)

metrics = pd.DataFrame(metric_rows)
metrics.to_csv(CURRENT / "simpson_raw_patient_metrics.csv", index=False)

# The nucleotide master created by analysis_10 stores counts as T1_raw_count_sum and
# T2_raw_count_sum. Standardize those aliases here for the final presentation table.
# Fractions are mandatory inputs; missing fractions indicate a broken upstream master
# and must not be silently replaced by a scalar NaN.
for c in ["T1_fraction", "T2_fraction"]:
    if c not in master.columns:
        raise RuntimeError(f"Required master-table column is missing: {c}")
    master[c] = pd.to_numeric(master[c], errors="coerce").fillna(0.0)

count_aliases = {
    "T1_count": "T1_raw_count_sum",
    "T2_count": "T2_raw_count_sum",
}
for target, source in count_aliases.items():
    if target in master.columns:
        master[target] = pd.to_numeric(master[target], errors="coerce").fillna(0.0)
    elif source in master.columns:
        master[target] = pd.to_numeric(master[source], errors="coerce").fillna(0.0)
    else:
        raise RuntimeError(f"Required count information is missing: neither {target} nor {source} exists")

raw = master.copy()
raw["raw_log2_T2_over_T1"] = np.where(
    (raw.T1_fraction > 0) & (raw.T2_fraction > 0),
    np.log2(raw.T2_fraction / raw.T1_fraction),
    np.nan,
)
raw["raw_absolute_delta"] = raw.T2_fraction - raw.T1_fraction
raw["raw_fold_contraction"] = np.where(raw.T2_fraction > 0, raw.T1_fraction / raw.T2_fraction, np.nan)
raw["raw_T1_expanded_mass_share"] = np.nan

for patient, idx in raw.groupby("patient").groups.items():
    idx = list(idx)
    p = raw.loc[idx]
    mask = p.T1_fraction.ge(PRIMARY)
    total = float(p.loc[mask, "T1_fraction"].sum())
    if total > 0:
        raw.loc[p.index[mask], "raw_T1_expanded_mass_share"] = p.loc[mask, "T1_fraction"] / total


def direction(row):
    if row.T1_fraction < PRIMARY:
        return "not primary T1-expanded"
    if row.T2_fraction <= 0:
        return "not detected"
    x = row.raw_log2_T2_over_T1
    if x >= 1:
        return "increased"
    if x <= -1:
        return "decreased"
    return "stable"

raw["raw_direction_0_10_percent"] = raw.apply(direction, axis=1)
raw["raw_strong_contraction"] = (
    raw.T1_fraction.ge(PRIMARY)
    & raw.T2_fraction.gt(0)
    & raw.T2_fraction.le(0.2 * raw.T1_fraction)
)

threshold_rows = []
for patient in paired:
    p = raw[raw.patient.astype(str).eq(patient)].copy()
    age = str(p.age_group.iloc[0])
    for pct, th in zip(THRESHOLD_PCT, THRESHOLDS):
        e = p[p.T1_fraction.ge(th)].copy()
        base = {
            "branch": "simpson_no_downsampling",
            "patient": patient,
            "age_group": age,
            "threshold_fraction": th,
            "threshold_percent": pct,
        }
        if e.empty:
            threshold_rows.append({
                **base,
                "n_T1_expanded": 0,
                "n_T2_expanded": int(p.T2_fraction.ge(th).sum()),
                "retained_count_percent": np.nan,
                "below_count_percent": np.nan,
                "not_detected_count_percent": np.nan,
                "retained_mass_percent": np.nan,
                "below_mass_percent": np.nan,
                "not_detected_mass_percent": np.nan,
            })
            continue
        retained = e.T2_fraction.ge(th)
        detected = e.T2_fraction.gt(0)
        below = detected & ~retained
        lost = ~detected
        mass_total = float(e.T1_fraction.sum())
        threshold_rows.append({
            **base,
            "n_T1_expanded": int(len(e)),
            "n_T2_expanded": int(p.T2_fraction.ge(th).sum()),
            "retained_count_percent": 100 * float(retained.mean()),
            "below_count_percent": 100 * float(below.mean()),
            "not_detected_count_percent": 100 * float(lost.mean()),
            "retained_mass_percent": 100 * float(e.loc[retained, "T1_fraction"].sum()) / mass_total,
            "below_mass_percent": 100 * float(e.loc[below, "T1_fraction"].sum()) / mass_total,
            "not_detected_mass_percent": 100 * float(e.loc[lost, "T1_fraction"].sum()) / mass_total,
        })

threshold = pd.DataFrame(threshold_rows)
observed = {round(float(x), 8) for x in threshold.threshold_percent.dropna().unique()}
if observed != set(THRESHOLD_PCT):
    raise RuntimeError(f"Unexpected Simpson threshold grid: {sorted(observed)}")
counts_per_threshold = threshold.groupby("threshold_percent").patient.nunique().to_dict()
if any(v != 22 for v in counts_per_threshold.values()) or len(counts_per_threshold) != len(THRESHOLD_PCT):
    raise RuntimeError(f"Incomplete Simpson threshold coverage: {counts_per_threshold}")
threshold.to_csv(CURRENT / "simpson_raw_threshold_sensitivity.csv", index=False)

keep = raw.T1_fraction.ge(PRIMARY)
raw.loc[keep].to_csv(CURRENT / "simpson_raw_primary_candidates.csv", index=False)

print("No-downsampling Simpson branch completed")
print("Paired patients:", len(paired))
print("Primary T1-expanded candidates:", int(keep.sum()))
print("Thresholds percent:", THRESHOLD_PCT)
