from pathlib import Path
import os

import numpy as np
import pandas as pd
import requests

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"

# The notebook runs analysis_01_prepare.py first. Loading the persisted table here
# also makes this script robust to namespace changes between notebook sections.
clones_path = RESULTS / "combined_trb_clonotypes.csv"
if not clones_path.exists():
    raise FileNotFoundError("Run analysis_01_prepare.py first: combined_trb_clonotypes.csv is missing.")
clones = pd.read_csv(clones_path)
clones["count"] = pd.to_numeric(clones["count"], errors="coerce").fillna(0)
clones["analysis_timepoint"] = pd.to_numeric(clones["analysis_timepoint"], errors="coerce")

IEDB_URL = "https://raw.githubusercontent.com/IEDB/TCRMatch/master/data/IEDB_data.tsv"
iedb_override = os.environ.get("MEASLES_IEDB_FILE", "").strip()
iedb_path = Path(iedb_override).expanduser() if iedb_override else ROOT / "resources" / "IEDB_data.tsv"
iedb_path.parent.mkdir(parents=True, exist_ok=True)

if not iedb_path.exists():
    response = requests.get(IEDB_URL, timeout=120)
    response.raise_for_status()
    iedb_path.write_bytes(response.content)

iedb = pd.read_csv(iedb_path, sep="\t")
expected = ["trimmed_seq", "original_sequence", "receptor_group", "epitope", "organism", "antigen"]
if not set(expected).issubset(iedb.columns):
    if iedb.shape[1] >= 6:
        iedb = iedb.iloc[:, :6].copy()
        iedb.columns = expected
    else:
        raise ValueError(f"Unexpected IEDB TCRMatch table columns: {list(iedb.columns)}")


def trim_cdr3(seq):
    seq = str(seq).strip()
    if len(seq) >= 2 and seq.startswith("C") and seq[-1] in {"F", "W"}:
        return seq[1:-1]
    return seq


def pathogen_label(row):
    text = (str(row.get("organism", "")) + " " + str(row.get("antigen", ""))).lower()
    if "measles" in text:
        return "Measles"
    if "hepatitis b" in text or "hepatitis b virus" in text or "hbv" in text:
        return "Hepatitis B"
    if "poliovirus" in text or "polio virus" in text:
        return "Poliovirus"
    if "clostridium tetani" in text or "tetanus" in text:
        return "Tetanus"
    return np.nan


iedb["pathogen"] = iedb.apply(pathogen_label, axis=1)
reference = iedb[iedb["pathogen"].notna()].copy()
reference["cdr3_trimmed"] = reference["trimmed_seq"].astype(str).str.strip()
reference = reference.drop_duplicates(["cdr3_trimmed", "pathogen", "epitope", "organism", "antigen"])
reference.to_csv(RESULTS / "antigen_reference_used.csv", index=False)

query = clones.copy()
query["cdr3_trimmed"] = query["cdr3aa"].map(trim_cdr3)
matches = query.merge(
    reference[["cdr3_trimmed", "pathogen", "epitope", "organism", "antigen"]],
    on="cdr3_trimmed",
    how="inner",
)
matched_unique = (
    matches.sort_values(["sample_id", "pathogen", "cdr3aa"])
    .drop_duplicates(["sample_id", "pathogen", "cdr3aa"])
)

sample_meta = clones.groupby("sample_id", as_index=False).agg(
    patient_id=("patient_id", "first"),
    age_group=("age_group", "first"),
    stage=("stage", "first"),
    analysis_timepoint=("analysis_timepoint", "first"),
    total_clonotypes=("cdr3aa", "nunique"),
    total_reads=("count", "sum"),
)

pathogens = ["Measles", "Hepatitis B", "Poliovirus", "Tetanus"]
summary_rows = []
for s in sample_meta.itertuples(index=False):
    for pathogen in pathogens:
        sub = matched_unique[(matched_unique.sample_id == s.sample_id) & (matched_unique.pathogen == pathogen)]
        matched_n = int(sub["cdr3aa"].nunique())
        matched_reads = float(sub.drop_duplicates("cdr3aa")["count"].sum())
        summary_rows.append({
            "sample_id": s.sample_id,
            "patient_id": s.patient_id,
            "age_group": s.age_group,
            "stage": s.stage,
            "analysis_timepoint": s.analysis_timepoint,
            "pathogen": pathogen,
            "matched_clonotypes": matched_n,
            "matched_clonotype_fraction": matched_n / s.total_clonotypes if s.total_clonotypes else np.nan,
            "matched_read_fraction": matched_reads / s.total_reads if s.total_reads else np.nan,
        })

antigen = pd.DataFrame(summary_rows)
antigen.to_csv(RESULTS / "antigen_specificity_by_sample.csv", index=False)
matches.to_csv(RESULTS / "antigen_specific_clonotypes.csv", index=False)

print("IEDB records used:", len(reference))
print("Exact sample-clonotype-pathogen matches:", len(matched_unique))
