#!/usr/bin/env python3
"""Exact CDR3 beta matching against the local IEDB TCRMatch reference"""
from pathlib import Path
import os

import numpy as np
import pandas as pd
import requests

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
clones = pd.read_csv(RESULTS / "combined_trb_clonotypes.csv")
clones["fraction"] = pd.to_numeric(clones["fraction"], errors="coerce").fillna(0.0)

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
    if iedb.shape[1] < 6:
        raise ValueError("Unexpected IEDB TCRMatch table format")
    iedb = iedb.iloc[:, :6].copy()
    iedb.columns = expected


def trim_cdr3(seq):
    seq = str(seq).strip()
    if len(seq) >= 2 and seq.startswith("C") and seq[-1] in {"F", "W"}:
        return seq[1:-1]
    return seq


def pathogen_label(row):
    text = (str(row.get("organism", "")) + " " + str(row.get("antigen", ""))).lower()
    if "measles" in text:
        return "Measles"
    if "hepatitis b" in text or "hbv" in text:
        return "Hepatitis B"
    if "poliovirus" in text or "polio virus" in text:
        return "Poliovirus"
    if "clostridium tetani" in text or "tetanus" in text:
        return "Tetanus"
    return np.nan

iedb["pathogen"] = iedb.apply(pathogen_label, axis=1)
reference = iedb[iedb.pathogen.notna()].copy()
reference["cdr3_trimmed"] = reference.trimmed_seq.astype(str).str.strip()
reference = reference.drop_duplicates(["cdr3_trimmed", "pathogen", "epitope", "organism", "antigen"])
reference.to_csv(RESULTS / "antigen_reference_used.csv", index=False)

query = clones.copy()
query["cdr3_trimmed"] = query.cdr3aa.map(trim_cdr3)
matches = query.merge(
    reference[["cdr3_trimmed", "pathogen", "epitope", "organism", "antigen"]],
    on="cdr3_trimmed",
    how="inner",
)
matches.to_csv(RESULTS / "antigen_specific_clonotypes.csv", index=False)

pathogens = ["Measles", "Hepatitis B", "Poliovirus", "Tetanus"]
rows = []
for sample_id, s in clones.groupby("sample_id"):
    for pathogen in pathogens:
        sub = matches[(matches.sample_id == sample_id) & matches.pathogen.eq(pathogen)].drop_duplicates(["cdr3aa", "v_gene"])
        rows.append({
            "sample_id": sample_id,
            "patient_id": s.patient_id.iloc[0],
            "clinical_timepoint": s.clinical_timepoint.iloc[0],
            "age_group": s.age_group.iloc[0],
            "pathogen": pathogen,
            "matched_clonotypes": int(len(sub)),
            "matched_read_fraction": float(sub.fraction.sum()) if len(sub) else 0.0,
        })
pd.DataFrame(rows).to_csv(RESULTS / "antigen_specificity_by_sample.csv", index=False)
print(f"IEDB records used: {len(reference)}")
print(f"Exact sample clonotype reference matches: {len(matches)}")
