#!/usr/bin/env python3
"""Estimate clone fate and sampling uncertainty for 22 paired T1 and T2 repertoires.

Each library is sampled without replacement to the minimum productive-read
depth within its patient pair. Replicate frequencies are averaged after
sampling. At least 500 iterations are required; the default seed is 1031.

The unchanged-frequency null averages exact conditional detection probabilities
across T1 draws. Expanded-set Jaccard measures intersection divided by union.
Read downsampling controls depth, while PCR bias remains uncorrected without UMI."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import os
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
CURRENT.mkdir(parents=True, exist_ok=True)

N_ITER = int(os.environ.get("MEASLES_RAREFACTION_ITERATIONS", "500"))
SEED = int(os.environ.get("MEASLES_RAREFACTION_SEED", "1031"))
if N_ITER < 500:
    raise RuntimeError("Pair-specific rarefaction requires at least 500 iterations")
THRESHOLD_PCT = [0.05, 0.10] + [round(x / 10, 2) for x in range(2, 16)]
THRESHOLDS = [x / 100.0 for x in THRESHOLD_PCT]
PRIMARY = 0.001
ROBUST_P = 0.99
KEY = ["cdr3nt", "cdr3aa", "v_gene", "j_gene"]

nt_path = RESULTS / "library_trb_clonotypes_nt.csv"
meta_path = CURRENT / "sample_metadata_current.csv"
master_path = CURRENT / "paired_T1_clone_master_nt.csv"
if not nt_path.exists() or not master_path.exists():
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "analysis_10_clone_fate_v2.py")],
        check=True,
    )

lib = pd.read_csv(nt_path)
meta = pd.read_csv(meta_path)
master = pd.read_csv(master_path)
meta["analysis_include_final"] = (
    meta["analysis_include_final"].astype(str).str.lower().isin(["true", "1", "yes"])
)
lib["count"] = (
    pd.to_numeric(lib["count"], errors="coerce").fillna(0).round().astype(np.int64)
)

case = meta[meta.analysis_include_final & meta.disease_status.eq("measles")].copy()
paired = []
for patient, z in case.groupby("patient_id"):
    if {"T1", "T2"}.issubset(set(z.clinical_timepoint.astype(str))):
        paired.append(str(patient))
paired = sorted(paired)
if len(paired) != 22:
    raise RuntimeError(f"Expected 22 T1-T2 patients, found {len(paired)}")

rng = np.random.default_rng(SEED)


def summarize(values):
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if len(finite_values) == 0:
        return np.nan, np.nan, np.nan
    return (
        float(np.median(finite_values)),
        float(np.quantile(finite_values, 0.025)),
        float(np.quantile(finite_values, 0.975)),
    )


def repertoire_metrics(f):
    positive_frequencies = f[f > 0]
    richness = len(positive_frequencies)
    shannon = (
        float(-(positive_frequencies * np.log(positive_frequencies)).sum())
        if richness
        else np.nan
    )
    clonality = float(1 - shannon / np.log(richness)) if richness > 1 else np.nan
    simpson_d = (
        float(np.sum(positive_frequencies * positive_frequencies))
        if richness
        else np.nan
    )
    return richness, shannon, clonality, simpson_d


if "P_detection_if_unchanged_analytic_raw_depth" not in master.columns:
    master["P_detection_if_unchanged_analytic_raw_depth"] = pd.to_numeric(
        master.get("P_detection_if_unchanged", np.nan), errors="coerce"
    )

pair_summaries = []
threshold_rows = []
clone_null_rows = []
convergence_rows = []

for patient in paired:
    pmeta = case[
        case.patient_id.astype(str).eq(patient)
        & case.clinical_timepoint.isin(["T1", "T2"])
    ].copy()
    point_libraries = {
        point: pmeta.loc[pmeta.clinical_timepoint.eq(point), "library_id"]
        .astype(str)
        .tolist()
        for point in ["T1", "T2"]
    }
    all_libraries = point_libraries["T1"] + point_libraries["T2"]
    z = lib[lib.library_id.astype(str).isin(all_libraries)].copy()

    keys = sorted(set(map(tuple, z[KEY].astype(str).to_numpy())))
    key_to_i = {k: i for i, k in enumerate(keys)}
    vectors = {}
    depths = {}
    for library_id, zz in z.groupby("library_id"):
        vec = np.zeros(len(keys), dtype=np.int64)
        for r in zz.itertuples(index=False):
            key = (str(r.cdr3nt), str(r.cdr3aa), str(r.v_gene), str(r.j_gene))
            vec[key_to_i[key]] += int(r.count)
        vectors[str(library_id)] = vec
        depths[str(library_id)] = int(vec.sum())

    missing = [x for x in all_libraries if x not in vectors]
    if missing:
        raise RuntimeError(f"{patient}: missing nucleotide count vectors for {missing}")
    pair_depth = min(depths[x] for x in all_libraries)
    if pair_depth <= 0:
        raise RuntimeError(f"{patient}: non-positive pair rarefaction depth")

    def rarefy_library(library_id):
        vec = vectors[library_id]
        total = int(vec.sum())
        if total == pair_depth:
            return vec.copy()
        return rng.multivariate_hypergeometric(vec, pair_depth).astype(np.int64)

    def point_frequency(point):
        fractions = []
        for library_id in point_libraries[point]:
            sampled = rarefy_library(library_id)
            fractions.append(sampled / float(pair_depth))
        f = np.mean(np.stack(fractions, axis=0), axis=0)
        total = float(f.sum())
        return f / total if total > 0 else f

    pmaster = master[master.patient.astype(str).eq(patient)].copy()
    master_indices = []
    key_indices = []
    for idx, r in pmaster.iterrows():
        key = (str(r.cdr3nt), str(r.cdr3aa), str(r.TRBV), str(r.TRBJ))
        if key not in key_to_i:
            raise RuntimeError(
                f"{patient}: master clonotype not found in pair key space: {r.clone_id_nt}"
            )
        master_indices.append(idx)
        key_indices.append(key_to_i[key])
    key_indices = np.asarray(key_indices, dtype=int)

    null_detection_probability_sum = np.zeros(len(master_indices), dtype=float)
    observed_t2_detect_hits = np.zeros(len(master_indices), dtype=np.int64)
    rare_t1_store = np.zeros((N_ITER, len(master_indices)), dtype=np.float32)
    rare_t2_store = np.zeros((N_ITER, len(master_indices)), dtype=np.float32)
    metric_store = defaultdict(list)
    threshold_store = {th: defaultdict(list) for th in THRESHOLDS}

    n_t2_libraries = len(point_libraries["T2"])
    null_total_reads = pair_depth * n_t2_libraries

    for iteration in range(N_ITER):
        f1 = point_frequency("T1")
        f2 = point_frequency("T2")
        rare_t1_store[iteration, :] = f1[key_indices]
        rare_t2_store[iteration, :] = f2[key_indices]
        observed_t2_detect_hits += (f2[key_indices] > 0).astype(np.int64)

        p_null = 1.0 - np.power(
            1.0 - np.clip(f1[key_indices], 0.0, 1.0), null_total_reads
        )
        # Averaging exact probabilities avoids an extra Bernoulli draw;
        # variation in this null comes only from rarefied T1 frequencies.
        null_detection_probability_sum += p_null

        r1, h1, c1, s1 = repertoire_metrics(f1)
        r2, h2, c2, s2 = repertoire_metrics(f2)
        for name, value in [
            ("T1_richness", r1),
            ("T2_richness", r2),
            ("T1_shannon", h1),
            ("T2_shannon", h2),
            ("T1_clonality", c1),
            ("T2_clonality", c2),
            ("T1_simpson_D", s1),
            ("T2_simpson_D", s2),
        ]:
            metric_store[name].append(value)

        all1 = set(np.flatnonzero(f1 > 0))
        all2 = set(np.flatnonzero(f2 > 0))
        union = all1 | all2
        complete_jaccard = len(all1 & all2) / len(union) if union else np.nan
        metric_store["complete_jaccard"].append(complete_jaccard)

        for th in THRESHOLDS:
            e1 = np.flatnonzero(f1 >= th)
            e2 = set(np.flatnonzero(f2 >= th))
            threshold_store[th]["n_T1_expanded"].append(float(len(e1)))
            threshold_store[th]["n_T2_expanded"].append(float(len(e2)))
            eu = set(e1) | e2
            threshold_store[th]["expanded_jaccard"].append(
                len(set(e1) & e2) / len(eu) if eu else np.nan
            )
            if len(e1) == 0:
                for name in [
                    "retained_count_percent",
                    "below_count_percent",
                    "not_detected_count_percent",
                    "retained_mass_percent",
                    "below_mass_percent",
                    "not_detected_mass_percent",
                ]:
                    threshold_store[th][name].append(np.nan)
                if np.isclose(th, PRIMARY):
                    convergence_rows.append(
                        {
                            "patient": patient,
                            "age_group": str(pmeta.age_group.iloc[0]),
                            "iteration": iteration + 1,
                            "pair_rarefaction_depth": pair_depth,
                            "not_detected_mass_percent": np.nan,
                            "expanded_jaccard": (
                                len(set(e1) & e2) / len(eu) if eu else np.nan
                            ),
                            "complete_jaccard": complete_jaccard,
                            "T1_simpson_D": s1,
                            "T2_simpson_D": s2,
                        }
                    )
                continue

            retained_mask = np.array([i in e2 for i in e1], dtype=bool)
            detected_mask = f2[e1] > 0
            below_mask = detected_mask & ~retained_mask
            lost_mask = ~detected_mask
            threshold_store[th]["retained_count_percent"].append(
                100 * retained_mask.mean()
            )
            threshold_store[th]["below_count_percent"].append(100 * below_mask.mean())
            threshold_store[th]["not_detected_count_percent"].append(
                100 * lost_mask.mean()
            )

            mass = f1[e1]
            mass_total = float(mass.sum())
            if mass_total > 0:
                threshold_store[th]["retained_mass_percent"].append(
                    100 * float(mass[retained_mask].sum()) / mass_total
                )
                threshold_store[th]["below_mass_percent"].append(
                    100 * float(mass[below_mask].sum()) / mass_total
                )
                not_detected_mass = 100 * float(mass[lost_mask].sum()) / mass_total
                threshold_store[th]["not_detected_mass_percent"].append(
                    not_detected_mass
                )
            else:
                not_detected_mass = np.nan
                for name in [
                    "retained_mass_percent",
                    "below_mass_percent",
                    "not_detected_mass_percent",
                ]:
                    threshold_store[th][name].append(np.nan)
            if np.isclose(th, PRIMARY):
                convergence_rows.append(
                    {
                        "patient": patient,
                        "age_group": str(pmeta.age_group.iloc[0]),
                        "iteration": iteration + 1,
                        "pair_rarefaction_depth": pair_depth,
                        "not_detected_mass_percent": not_detected_mass,
                        "expanded_jaccard": (
                            len(set(e1) & e2) / len(eu) if eu else np.nan
                        ),
                        "complete_jaccard": complete_jaccard,
                        "T1_simpson_D": s1,
                        "T2_simpson_D": s2,
                    }
                )

    summary = {
        "patient": patient,
        "age_group": str(pmeta.age_group.iloc[0]),
        "pair_rarefaction_depth": pair_depth,
        "iterations": N_ITER,
        "n_T1_libraries": len(point_libraries["T1"]),
        "n_T2_libraries": len(point_libraries["T2"]),
        "T1_libraries": ";".join(point_libraries["T1"]),
        "T2_libraries": ";".join(point_libraries["T2"]),
    }
    for name, values in metric_store.items():
        med, lo, hi = summarize(values)
        summary[f"{name}_median"] = med
        summary[f"{name}_q025"] = lo
        summary[f"{name}_q975"] = hi
    pair_summaries.append(summary)

    for th in THRESHOLDS:
        row = {
            "patient": patient,
            "age_group": str(pmeta.age_group.iloc[0]),
            "threshold_fraction": th,
            "threshold_percent": 100 * th,
            "pair_rarefaction_depth": pair_depth,
            "iterations": N_ITER,
        }
        for name, values in threshold_store[th].items():
            med, lo, hi = summarize(values)
            row[f"{name}_median"] = med
            row[f"{name}_q025"] = lo
            row[f"{name}_q975"] = hi
        threshold_rows.append(row)

    empirical_p = null_detection_probability_sum / float(N_ITER)
    observed_detect_p = observed_t2_detect_hits / float(N_ITER)
    t1_median = np.median(rare_t1_store, axis=0)
    t2_median = np.median(rare_t2_store, axis=0)
    zero_upper_pair = 1.0 - 0.05 ** (1.0 / max(null_total_reads, 1))

    for pos, master_idx in enumerate(master_indices):
        clone_null_rows.append(
            {
                "patient": patient,
                "clone_id_nt": str(master.at[master_idx, "clone_id_nt"]),
                "pair_rarefaction_depth": pair_depth,
                "iterations": N_ITER,
                "n_T2_libraries": n_t2_libraries,
                "null_total_reads": null_total_reads,
                "P_detection_if_unchanged_empirical": empirical_p[pos],
                "P_detection_in_observed_T2_rarefactions": observed_detect_p[pos],
                "T1_fraction_rarefied_median": float(t1_median[pos]),
                "T2_fraction_rarefied_median": float(t2_median[pos]),
                "T2_zero_upper_95_pair_depth": zero_upper_pair,
            }
        )
        master.at[master_idx, "pair_rarefaction_depth"] = pair_depth
        master.at[master_idx, "rarefaction_iterations"] = N_ITER
        master.at[master_idx, "P_detection_if_unchanged_empirical"] = empirical_p[pos]
        master.at[master_idx, "P_detection_if_unchanged"] = empirical_p[pos]
        master.at[master_idx, "P_detection_in_observed_T2_rarefactions"] = (
            observed_detect_p[pos]
        )
        master.at[master_idx, "T1_fraction_rarefied_median"] = float(t1_median[pos])
        master.at[master_idx, "T2_fraction_rarefied_median"] = float(t2_median[pos])
        master.at[master_idx, "T2_zero_upper_95_pair_depth"] = zero_upper_pair

pair_summary = pd.DataFrame(pair_summaries)
threshold_summary = pd.DataFrame(threshold_rows)
clone_null = pd.DataFrame(clone_null_rows)
convergence = pd.DataFrame(convergence_rows)
pair_summary.to_csv(CURRENT / "pairwise_rarefaction_summary.csv", index=False)
threshold_summary.to_csv(CURRENT / "pairwise_threshold_sensitivity.csv", index=False)
clone_null.to_csv(CURRENT / "clone_sampling_null.csv", index=False)
convergence.to_csv(CURRENT / "downsampling_iteration_convergence.csv", index=False)

expanded_primary = (
    master["expanded_T1_0_10_percent"]
    .astype(str)
    .str.lower()
    .isin(["true", "1", "yes"])
)
t2_zero = pd.to_numeric(master["T2_fraction"], errors="coerce").fillna(0).eq(0)
master["sampling_robust_not_detected"] = (
    expanded_primary
    & t2_zero
    & pd.to_numeric(master["P_detection_if_unchanged_empirical"], errors="coerce").ge(
        ROBUST_P
    )
)
strong_contraction = (
    expanded_primary
    & pd.to_numeric(master["T2_fraction"], errors="coerce").gt(0)
    & (
        pd.to_numeric(master["T2_fraction"], errors="coerce")
        <= 0.2 * pd.to_numeric(master["T1_fraction"], errors="coerce")
    )
)
master["sampling_candidate_class"] = np.select(
    [
        master["sampling_robust_not_detected"],
        strong_contraction,
        master["fate_0_10_percent"].eq("detected below 0.10 percent"),
        master["fate_0_10_percent"].eq("retained expanded"),
    ],
    [
        "sampling-robust not detected",
        "strong contraction",
        "detected below primary threshold",
        "retained expanded",
    ],
    default="other",
)
master.to_csv(master_path, index=False)

(CURRENT / "pairwise_rarefaction_README.txt").write_text(
    "\n".join(
        [
            "Pair-specific repeated rarefaction for the 22 accepted T1-T2 patient pairs",
            "Depth rule: each patient has its own pair depth equal to the minimum productive-read depth among all accepted T1 and T2 libraries for that patient",
            f"Iterations: {N_ITER}",
            f"Seed: {SEED}",
            "Observed libraries: NumPy multivariate hypergeometric sampling without replacement",
            "Biological replicate libraries: rarefied separately, converted to frequencies, then averaged within a time point",
            "Internal null: after each T1 rarefaction, exact clone-detection probability at the same T2 sampling effort is calculated from the rarefied T1 frequency",
            "P_detection_if_unchanged is the mean exact conditional detection probability across T1 rarefactions; no extra Bernoulli noise is added",
            "Analytical raw-depth probability is retained separately in P_detection_if_unchanged_analytic_raw_depth",
            "Threshold sensitivity: 0.05, 0.10, then 0.20 to 1.50 percent in 0.10 percentage-point steps",
            "This controls sequencing read depth only and cannot remove PCR amplification bias because UMI are absent",
        ]
    )
    + "\n",
    encoding="utf-8",
)

print("Pair-specific rarefaction completed")
print("Paired patients", len(pair_summary))
print("Iterations", N_ITER)
print("Thresholds percent", THRESHOLD_PCT)
print(
    "Pair depths",
    int(pair_summary.pair_rarefaction_depth.min()),
    "to",
    int(pair_summary.pair_rarefaction_depth.max()),
)
