#!/usr/bin/env python3
"""Prepare patient TRB repertoires from config/analysis_metadata.csv.

Apply library QC before aggregation. Average normalized clonotype frequencies
across accepted libraries at each patient time point to give replicates equal weight."""

from pathlib import Path
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
CONFIG = ROOT / "config"
RESULTS.mkdir(exist_ok=True)
CURRENT.mkdir(parents=True, exist_ok=True)

QC_MIN_PRODUCTIVE_CLONOTYPES = 4400
EXPANDED_FRACTION = 0.001

manifest = pd.read_csv(CONFIG / "sample_manifest.csv")
meta = pd.read_csv(CONFIG / "analysis_metadata.csv")

for col in ["include", "is_control"]:
    manifest[col] = manifest[col].astype(str).str.lower().isin(["true", "1", "yes"])
meta["passes_qc_4400"] = (
    meta["passes_qc_4400"].astype(str).str.lower().isin(["true", "1", "yes"])
)
for col in [
    "productive_trb_reads",
    "productive_trb_clonotypes",
    "largest_clone_percent",
]:
    meta[col] = pd.to_numeric(meta[col], errors="coerce")

joined = manifest.merge(
    meta,
    on=[
        "library_id",
        "patient_id",
        "clinical_timepoint",
        "disease_status",
        "age_group",
        "batch",
    ],
    how="outer",
    suffixes=("_manifest", ""),
    indicator=True,
)
if not joined["_merge"].eq("both").all():
    bad = joined.loc[~joined["_merge"].eq("both"), ["library_id", "_merge"]]
    raise RuntimeError(
        "Manifest and analysis metadata are inconsistent\n" + bad.to_string(index=False)
    )

joined["analysis_include_final"] = (
    joined["include"]
    & joined["passes_qc_4400"]
    & joined["clinical_timepoint"].isin(["CONTROL", "T1", "T2", "T3"])
    & joined["disease_status"].isin(["healthy", "measles"])
)
joined["final_exclusion_reason"] = ""
joined.loc[~joined["include"], "final_exclusion_reason"] = joined.loc[
    ~joined["include"], "exclusion_reason"
].fillna("")
low = joined["include"] & ~joined["passes_qc_4400"]
joined.loc[low, "final_exclusion_reason"] = (
    "Productive TRB clonotypes below 4400 or library level quality control failure"
)

# Controls have no longitudinal study point
bad_controls = joined[
    joined["is_control"] & ~joined["clinical_timepoint"].eq("CONTROL")
]
if not bad_controls.empty:
    raise RuntimeError("Healthy controls must use CONTROL only")
if joined["clinical_timepoint"].astype(str).eq("T0").any():
    raise RuntimeError("T0 does not exist in this study")

# Explicit resolved biological repeat
rep = joined[joined["replicate_group"].fillna("").eq("БА1_071_T2")]
if (
    set(rep["library_id"]) != {"new_106", "new_111"}
    or not rep["analysis_include_final"].all()
):
    raise RuntimeError("БА1_071 T2 replicate definition is incomplete")

joined.to_csv(CURRENT / "sample_metadata_current.csv", index=False)


def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def clean_gene(value):
    if pd.isna(value):
        return ""
    text = str(value).split(",")[0].split("(")[0].split("*")[0].strip()
    return text


library_rows = []
for r in joined.itertuples(index=False):
    path = RESULTS / "clonotypes" / f"{r.library_id}.tsv"
    if not path.exists():
        if r.analysis_include_final:
            raise FileNotFoundError(
                f"Missing clonotype table for included library {r.library_id}"
            )
        continue
    df = pd.read_csv(path, sep="\t")
    aa_col = pick_col(df, ["aaSeqCDR3", "cdr3aa", "cdr3_aa"])
    count_col = pick_col(df, ["cloneCount", "readCount", "count"])
    v_col = pick_col(df, ["allVHitsWithScore", "bestVHit", "vGene", "v_call"])
    j_col = pick_col(df, ["allJHitsWithScore", "bestJHit", "jGene", "j_call"])
    if aa_col is None or count_col is None:
        raise ValueError(f"Cannot identify CDR3 or count columns in {path.name}")
    small = pd.DataFrame(
        {
            "library_id": r.library_id,
            "patient_id": r.patient_id,
            "clinical_timepoint": r.clinical_timepoint,
            "disease_status": r.disease_status,
            "age_group": r.age_group,
            "batch": r.batch,
            "analysis_include_final": bool(r.analysis_include_final),
            "cdr3aa": df[aa_col].astype(str),
            "count": pd.to_numeric(df[count_col], errors="coerce").fillna(0.0),
            "v_gene": df[v_col].map(clean_gene) if v_col else "",
            "j_gene": df[j_col].map(clean_gene) if j_col else "",
        }
    )
    small = small[
        (small["cdr3aa"] != "nan") & (small["cdr3aa"] != "") & (small["count"] > 0)
    ].copy()
    library_rows.append(small)

if not library_rows:
    raise RuntimeError("No TRB clonotype tables were found")

library_clones = pd.concat(library_rows, ignore_index=True)
library_clones.to_csv(RESULTS / "library_trb_clonotypes.csv", index=False)

# Collapse duplicated aaV rows within each library before normalization
lib = (
    library_clones[library_clones["analysis_include_final"]]
    .groupby(
        [
            "library_id",
            "patient_id",
            "clinical_timepoint",
            "disease_status",
            "age_group",
            "batch",
            "cdr3aa",
            "v_gene",
        ],
        as_index=False,
    )["count"]
    .sum()
)
lib_totals = (
    lib.groupby("library_id", as_index=False)["count"]
    .sum()
    .rename(columns={"count": "library_total"})
)
lib = lib.merge(lib_totals, on="library_id", how="left")
lib["library_fraction"] = lib["count"] / lib["library_total"]

# For every patient and clinical point average normalized frequencies over accepted libraries
point_keys = ["patient_id", "clinical_timepoint", "disease_status", "age_group"]
point_libraries = (
    lib[point_keys + ["library_id", "library_total"]]
    .drop_duplicates()
    .groupby(point_keys, as_index=False)
    .agg(
        n_replicate_libraries=("library_id", "nunique"),
        mean_library_reads=("library_total", "mean"),
    )
)

freq = (
    lib.groupby(point_keys + ["cdr3aa", "v_gene"], as_index=False)
    .agg(sum_library_fraction=("library_fraction", "sum"))
    .merge(point_libraries, on=point_keys, how="left")
)
freq["fraction"] = freq["sum_library_fraction"] / freq["n_replicate_libraries"]
# Numerical guard, although averaging normalized libraries already sums to one
norm = freq.groupby(point_keys)["fraction"].transform("sum")
freq["fraction"] = np.where(norm > 0, freq["fraction"] / norm, 0.0)
freq["count"] = freq["fraction"] * freq["mean_library_reads"]
freq["expanded"] = freq["fraction"] >= EXPANDED_FRACTION
freq["sample_id"] = np.where(
    freq["clinical_timepoint"].eq("CONTROL"),
    freq["patient_id"] + "__CONTROL",
    freq["patient_id"] + "__" + freq["clinical_timepoint"],
)
combined = freq[
    [
        "sample_id",
        "patient_id",
        "clinical_timepoint",
        "disease_status",
        "age_group",
        "cdr3aa",
        "v_gene",
        "count",
        "fraction",
        "expanded",
        "n_replicate_libraries",
        "mean_library_reads",
    ]
].copy()
combined.to_csv(RESULTS / "combined_trb_clonotypes.csv", index=False)


def repertoire_metrics(z):
    p = z["fraction"].to_numpy(float)
    p = p[p > 0]
    richness = int(len(p))
    shannon = float(-(p * np.log(p)).sum()) if richness else np.nan
    clonality = float(1.0 - shannon / np.log(richness)) if richness > 1 else np.nan
    top10 = float(np.sort(p)[::-1][:10].sum()) if richness else np.nan
    exp_mask = z["fraction"] >= EXPANDED_FRACTION
    expanded_fraction = float(z.loc[exp_mask, "fraction"].sum())
    return pd.Series(
        {
            "richness": richness,
            "shannon": shannon,
            "clonality": clonality,
            "top10_fraction": top10,
            "n_expanded_clonotypes": int(exp_mask.sum()),
            "expanded_read_fraction": expanded_fraction,
            "n_replicate_libraries": int(z["n_replicate_libraries"].iloc[0]),
            "mean_library_reads": float(z["mean_library_reads"].iloc[0]),
        }
    )


metric_rows = []
for keys, z in combined.groupby(
    ["sample_id", "patient_id", "clinical_timepoint", "disease_status", "age_group"],
    sort=False,
):
    rec = dict(
        zip(
            [
                "sample_id",
                "patient_id",
                "clinical_timepoint",
                "disease_status",
                "age_group",
            ],
            keys,
        )
    )
    rec.update(repertoire_metrics(z).to_dict())
    metric_rows.append(rec)
metrics = pd.DataFrame(metric_rows)
metrics["expanded_threshold_fraction"] = EXPANDED_FRACTION
metrics.to_csv(RESULTS / "diversity_metrics.csv", index=False)

# A compact audit of the biological observation counts
cohort = metrics.groupby(["age_group", "clinical_timepoint"], as_index=False).agg(
    n_patient_observations=("patient_id", "nunique")
)
cohort.to_csv(CURRENT / "cohort_counts.csv", index=False)

print(
    f"Final library QC threshold: {QC_MIN_PRODUCTIVE_CLONOTYPES} productive TRB clonotypes"
)
print(
    f"Expanded clone threshold: {100 * EXPANDED_FRACTION:.2f} percent of productive TRB reads"
)
print(f"Included sequencing libraries: {int(joined.analysis_include_final.sum())}")
print(
    f"Biological patients: {joined.loc[joined.analysis_include_final, 'patient_id'].nunique()}"
)
print(f"Patient timepoint observations: {metrics.shape[0]}")
print("БА1_071 T2 libraries are averaged after within library normalization")
