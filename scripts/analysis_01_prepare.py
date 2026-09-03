from pathlib import Path
import os, math, warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests
import requests

warnings.filterwarnings("ignore", category=FutureWarning)
sns.set_context("notebook")

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)
SEED = 1031
PILOT = os.environ.get("MEASLES_PILOT", "0") == "1"

manifest_path = RESULTS / "sample_manifest_resolved.csv"
if not manifest_path.exists():
    raise FileNotFoundError("Run notebook 00 first: results/sample_manifest_resolved.csv is missing.")

manifest = pd.read_csv(manifest_path)
for col in ["include", "is_control", "fastq_ready"]:
    if col in manifest.columns:
        manifest[col] = manifest[col].astype(str).str.lower().isin(["true", "1", "yes"])

manifest["analysis_timepoint"] = pd.to_numeric(manifest["analysis_timepoint"], errors="coerce")

def stage_from_row(row):
    if bool(row["is_control"]):
        return "healthy"
    tp = row["analysis_timepoint"]
    if tp == 1:
        return "acute"
    if tp == 2:
        return "recovery"
    if tp in (3, 4):
        return "follow_up"
    return "unknown"

manifest["stage"] = manifest.apply(stage_from_row, axis=1)
manifest[["library_id","patient_id","age_group","stage","batch","include"]].head()

def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def clean_gene(value):
    if pd.isna(value):
        return np.nan
    text = str(value).split(",")[0].strip()
    text = text.split("(")[0].strip()
    text = text.split("*")[0].strip()
    return text if text else np.nan

rows = []
selected_manifest = manifest[(manifest["include"]) & (manifest["fastq_ready"])].copy()
if PILOT:
    pilot_ids = ["new_10", "new_11", "new_12", "legacy_UDP0107", "legacy_UDP0108"]
    selected_manifest = selected_manifest[selected_manifest.library_id.isin(pilot_ids)]

for r in selected_manifest.itertuples(index=False):
    path = RESULTS / "clonotypes" / f"{r.library_id}.tsv"
    if not path.exists():
        continue
    df = pd.read_csv(path, sep="\t")
    aa_col = pick_col(df, ["aaSeqCDR3", "cdr3aa", "cdr3_aa"])
    count_col = pick_col(df, ["cloneCount", "readCount", "count"])
    v_col = pick_col(df, ["allVHitsWithScore", "bestVHit", "vGene", "v_call"])
    j_col = pick_col(df, ["allJHitsWithScore", "bestJHit", "jGene", "j_call"])
    if aa_col is None or count_col is None:
        raise ValueError(f"Cannot identify CDR3/count columns in {path.name}: {list(df.columns)}")

    small = pd.DataFrame({
        "library_id": r.library_id,
        "patient_id": r.patient_id,
        "age_group": r.age_group,
        "stage": r.stage,
        "analysis_timepoint": r.analysis_timepoint,
        "batch": r.batch,
        "cdr3aa": df[aa_col].astype(str),
        "count": pd.to_numeric(df[count_col], errors="coerce").fillna(0),
        "v_gene": df[v_col].map(clean_gene) if v_col else np.nan,
        "j_gene": df[j_col].map(clean_gene) if j_col else np.nan,
    })
    small = small[(small["cdr3aa"].notna()) & (small["cdr3aa"] != "nan") & (small["count"] > 0)]
    rows.append(small)

if not rows:
    raise RuntimeError("No MiXCR TRB export tables were found. Run notebook 00 first.")

library_clones = pd.concat(rows, ignore_index=True)

# Sum technical libraries that resolve to one biological patient/timepoint.
# V/J are taken from the highest-count observation of a CDR3 when technical calls disagree.
best_annotation = (
    library_clones.sort_values("count", ascending=False)
    .drop_duplicates(["patient_id", "stage", "cdr3aa"])
    [["patient_id", "stage", "cdr3aa", "v_gene", "j_gene", "age_group", "analysis_timepoint"]]
)
counts = (
    library_clones.groupby(["patient_id", "stage", "cdr3aa"], as_index=False)["count"].sum()
)
clones = counts.merge(best_annotation, on=["patient_id", "stage", "cdr3aa"], how="left")
clones["sample_id"] = clones["patient_id"].astype(str) + "__" + clones["stage"].astype(str)

clones.to_csv(RESULTS / "combined_trb_clonotypes.csv", index=False)
print("Biological samples:", clones["sample_id"].nunique())
print("Unique CDR3β rows:", len(clones))

def shannon_from_counts(counts):
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts > 0]
    if counts.size == 0:
        return np.nan
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum())

def rarefy_counts(counts, depth, rng):
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total < depth or depth <= 0:
        return None
    p = counts / total
    return rng.multinomial(int(depth), p)

sample_totals = clones.groupby("sample_id")["count"].sum().sort_values()
if len(sample_totals) >= 3:
    rarefaction_depth = int(max(1, np.floor(np.quantile(sample_totals, 0.10))))
else:
    rarefaction_depth = int(max(1, sample_totals.min()))

rng = np.random.default_rng(SEED)
diversity_rows = []
for sample_id, sdf in clones.groupby("sample_id"):
    raw_counts = sdf["count"].to_numpy(dtype=float)
    total = int(raw_counts.sum())
    rare = rarefy_counts(raw_counts, rarefaction_depth, rng)
    if rare is None:
        richness = np.nan
        shannon = np.nan
        hill_q1 = np.nan
        rarefied_top10_fraction = np.nan
    else:
        richness = int((rare > 0).sum())
        shannon = shannon_from_counts(rare)
        hill_q1 = float(np.exp(shannon)) if np.isfinite(shannon) else np.nan
        rarefied_top10_fraction = float(np.sort(rare)[::-1][:10].sum() / rare.sum()) if rare.sum() > 0 else np.nan

    sorted_counts = np.sort(raw_counts)[::-1]
    top10_fraction = float(sorted_counts[:10].sum() / total) if total > 0 else np.nan
    meta = sdf.iloc[0]
    diversity_rows.append({
        "sample_id": sample_id,
        "patient_id": meta.patient_id,
        "age_group": meta.age_group,
        "stage": meta.stage,
        "analysis_timepoint": meta.analysis_timepoint,
        "productive_trb_reads": total,
        "raw_richness": int(len(sdf)),
        "rarefaction_depth": rarefaction_depth,
        "included_in_rarefaction": bool(total >= rarefaction_depth),
        "rarefied_richness": richness,
        "shannon": shannon,
        "hill_q1": hill_q1,
        "rarefied_top10_fraction": rarefied_top10_fraction,
        "top10_fraction": top10_fraction,
    })

diversity = pd.DataFrame(diversity_rows)
diversity.to_csv(RESULTS / "diversity_metrics.csv", index=False)
print("Common rarefaction depth:", rarefaction_depth)
print("Samples retained for normalized diversity:", int(diversity.included_in_rarefaction.sum()), "/", len(diversity))
