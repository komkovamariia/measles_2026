#!/usr/bin/env python3
"""Extended expansion-threshold sensitivity for both analytical branches.

Threshold grid in percent:
0.05, 0.10, then 0.20 to 1.50 in 0.10 percentage-point steps.

The rarefied branch repeats the pair-specific, without-replacement sampling rule used
by analysis_12. The raw branch uses observed normalized frequencies without downsampling.
The output is intentionally one compact long table used by the final publisher.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import os

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
N_ITER = int(os.environ.get("MEASLES_RAREFACTION_ITERATIONS", "500"))
SEED = int(os.environ.get("MEASLES_RAREFACTION_SEED", "1031"))
if N_ITER < 500:
    raise RuntimeError("Extended threshold analysis requires at least 500 rarefaction iterations")

THRESHOLD_PCT = [0.05, 0.10] + [round(x / 10, 2) for x in range(2, 16)]
THRESHOLDS = [x / 100.0 for x in THRESHOLD_PCT]
KEY = ["cdr3nt", "cdr3aa", "v_gene", "j_gene"]

meta = pd.read_csv(CURRENT / "sample_metadata_current.csv")
master = pd.read_csv(CURRENT / "paired_T1_clone_master_nt.csv")
lib = pd.read_csv(RESULTS / "library_trb_clonotypes_nt.csv")
meta["analysis_include_final"] = meta["analysis_include_final"].astype(str).str.lower().isin(["true", "1", "yes"])
lib["count"] = pd.to_numeric(lib["count"], errors="coerce").fillna(0).round().astype(np.int64)

case = meta[meta.analysis_include_final & meta.disease_status.eq("measles")].copy()
paired = sorted(
    str(patient)
    for patient, z in case.groupby("patient_id")
    if {"T1", "T2"}.issubset(set(z.clinical_timepoint.astype(str)))
)
if len(paired) != 22:
    raise RuntimeError(f"Expected 22 paired patients, found {len(paired)}")

rng = np.random.default_rng(SEED)
rows = []


def summarize(values):
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if len(a) == 0:
        return np.nan, np.nan, np.nan
    return float(np.median(a)), float(np.quantile(a, 0.025)), float(np.quantile(a, 0.975))


def metrics_for_pair(f1, f2, th):
    e1 = np.flatnonzero(f1 >= th)
    e2 = set(np.flatnonzero(f2 >= th))
    if len(e1) == 0:
        return {k: np.nan for k in [
            "expanded_jaccard", "retained_count_percent", "below_count_percent",
            "not_detected_count_percent", "retained_mass_percent", "below_mass_percent",
            "not_detected_mass_percent",
        ]}
    retained = np.array([i in e2 for i in e1], dtype=bool)
    detected = f2[e1] > 0
    below = detected & ~retained
    lost = ~detected
    union = set(e1) | e2
    mass = f1[e1]
    total = float(mass.sum())
    return {
        "expanded_jaccard": len(set(e1) & e2) / len(union) if union else np.nan,
        "retained_count_percent": 100 * float(retained.mean()),
        "below_count_percent": 100 * float(below.mean()),
        "not_detected_count_percent": 100 * float(lost.mean()),
        "retained_mass_percent": 100 * float(mass[retained].sum()) / total if total else np.nan,
        "below_mass_percent": 100 * float(mass[below].sum()) / total if total else np.nan,
        "not_detected_mass_percent": 100 * float(mass[lost].sum()) / total if total else np.nan,
    }


# Branch 1: pair-specific repeated rarefaction.
for patient in paired:
    pmeta = case[case.patient_id.astype(str).eq(patient) & case.clinical_timepoint.isin(["T1", "T2"])].copy()
    point_libraries = {
        point: pmeta.loc[pmeta.clinical_timepoint.eq(point), "library_id"].astype(str).tolist()
        for point in ["T1", "T2"]
    }
    all_libraries = point_libraries["T1"] + point_libraries["T2"]
    z = lib[lib.library_id.astype(str).isin(all_libraries)].copy()
    keys = sorted(set(map(tuple, z[KEY].astype(str).to_numpy())))
    key_to_i = {k: i for i, k in enumerate(keys)}
    vectors, depths = {}, {}
    for lid, zz in z.groupby("library_id"):
        vec = np.zeros(len(keys), dtype=np.int64)
        for r in zz.itertuples(index=False):
            k = (str(r.cdr3nt), str(r.cdr3aa), str(r.v_gene), str(r.j_gene))
            vec[key_to_i[k]] += int(r.count)
        vectors[str(lid)] = vec
        depths[str(lid)] = int(vec.sum())
    missing = [x for x in all_libraries if x not in vectors]
    if missing:
        raise RuntimeError(f"{patient}: missing clonotype vectors for {missing}")
    pair_depth = min(depths[x] for x in all_libraries)

    def point_frequency(point):
        parts = []
        for lid in point_libraries[point]:
            vec = vectors[lid]
            sampled = vec.copy() if int(vec.sum()) == pair_depth else rng.multivariate_hypergeometric(vec, pair_depth)
            parts.append(sampled / float(pair_depth))
        f = np.mean(np.stack(parts), axis=0)
        s = float(f.sum())
        return f / s if s > 0 else f

    store = {th: defaultdict(list) for th in THRESHOLDS}
    for _ in range(N_ITER):
        f1 = point_frequency("T1")
        f2 = point_frequency("T2")
        for th in THRESHOLDS:
            for name, value in metrics_for_pair(f1, f2, th).items():
                store[th][name].append(value)

    for pct, th in zip(THRESHOLD_PCT, THRESHOLDS):
        row = {
            "branch": "pairwise_rarefied",
            "patient": patient,
            "age_group": str(pmeta.age_group.iloc[0]),
            "threshold_fraction": th,
            "threshold_percent": pct,
            "pair_rarefaction_depth": pair_depth,
            "iterations": N_ITER,
        }
        for name, values in store[th].items():
            med, lo, hi = summarize(values)
            row[f"{name}_median"] = med
            row[f"{name}_q025"] = lo
            row[f"{name}_q975"] = hi
        rows.append(row)

# Branch 2: no downsampling, same threshold grid.
for patient in paired:
    p = master[master.patient.astype(str).eq(patient)].copy()
    if p.empty:
        continue
    f1 = pd.to_numeric(p["T1_fraction"], errors="coerce").fillna(0).to_numpy(float)
    f2 = pd.to_numeric(p["T2_fraction"], errors="coerce").fillna(0).to_numpy(float)
    age = str(p.age_group.iloc[0])
    for pct, th in zip(THRESHOLD_PCT, THRESHOLDS):
        vals = metrics_for_pair(f1, f2, th)
        row = {
            "branch": "simpson_no_downsampling",
            "patient": patient,
            "age_group": age,
            "threshold_fraction": th,
            "threshold_percent": pct,
            "pair_rarefaction_depth": np.nan,
            "iterations": 0,
        }
        for name, value in vals.items():
            row[f"{name}_median"] = value
            row[f"{name}_q025"] = value
            row[f"{name}_q975"] = value
        rows.append(row)

out = pd.DataFrame(rows)
expected = set(THRESHOLD_PCT)
observed = {round(float(x), 8) for x in out.threshold_percent.dropna().unique()}
if observed != expected:
    raise RuntimeError(f"Threshold grid mismatch: {sorted(observed)}")
if out.groupby("branch").patient.nunique().to_dict() != {
    "pairwise_rarefied": 22,
    "simpson_no_downsampling": 22,
}:
    raise RuntimeError("Extended threshold table does not contain both complete 22-patient branches")

out.to_csv(RESULTS / "expanded_threshold_sensitivity_extended.csv", index=False)
print("Extended threshold sensitivity completed")
print("Thresholds percent:", THRESHOLD_PCT)
print("Rows:", len(out))
