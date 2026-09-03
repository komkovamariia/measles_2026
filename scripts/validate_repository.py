#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def fail(msg):
    raise SystemExit(f"VALIDATION FAILED: {msg}")


required = [
    "README.md", "RUN_GUIDE.md", "environment.yml", "config/sample_manifest.csv",
    "config/metadata_known_issues.csv", "notebooks/00_prepare_and_mixcr.ipynb",
    "notebooks/01_measles_analysis.ipynb", "scripts/preprocess_human_trb.py",
    "scripts/analysis_01_prepare.py", "scripts/analysis_02_antigens.py",
    "scripts/analysis_03_figures.py", "scripts/analysis_04_stats.py",
    "scripts/analysis_05_longitudinal_venn.py", "scripts/run_pipeline.py",
    "scripts/run_slurm.sh", "run_all_from_scratch.sh",
]
missing = [p for p in required if not (ROOT / p).exists()]
if missing:
    fail(f"missing required files: {missing}")

# Source notebooks must be clean and valid JSON.
for rel in ["notebooks/00_prepare_and_mixcr.ipynb", "notebooks/01_measles_analysis.ipynb"]:
    nb = json.loads((ROOT / rel).read_text(encoding="utf-8"))
    ids = [c.get("id") for c in nb.get("cells", [])]
    if any(x is None for x in ids) or len(ids) != len(set(ids)):
        fail(f"invalid/duplicate notebook cell IDs: {rel}")
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code":
            if cell.get("outputs") or cell.get("execution_count") is not None:
                fail(f"source notebook contains outputs/execution count: {rel}")

# All active Python must parse.
for path in (ROOT / "scripts").glob("*.py"):
    try:
        ast.parse(path.read_text(encoding="utf-8"))
    except SyntaxError as exc:
        fail(f"Python syntax error in {path.name}: {exc}")

manifest = pd.read_csv(ROOT / "config/sample_manifest.csv")
new = manifest[manifest.batch.eq("new")].copy()
legacy = manifest[manifest.batch.eq("legacy")].copy()
if len(new) != 96:
    fail(f"expected 96 new study-table libraries, found {len(new)}")
if len(legacy) != 16:
    fail(f"expected 16 legacy UDP libraries, found {len(legacy)}")
if new.library_id.duplicated().any() or new.sequencing_index.duplicated().any():
    fail("new sequencing library IDs/indices must be unique")

# Patient identity must retain the full code including the numeric part after letters.
if set(new.loc[new.patient_id.str.startswith("СА1_", na=False), "patient_id"]) < {"СА1_013", "СА1_017"}:
    fail("patient IDs СА1_013 and СА1_017 are not kept as distinct people")

# The new sample index and final timepoint must agree with the study-table record.
for r in new.itertuples(index=False):
    if not str(r.raw_sample_id).startswith(str(r.sequencing_index) + "_"):
        fail(f"sequencing prefix mismatch: {r.library_id} {r.raw_sample_id}")
    m = re.search(r"\.([1-4])(?:[A-Za-zА-Яа-я])?$", str(r.raw_sample_id))
    if not m or int(m.group(1)) != int(float(r.analysis_timepoint)):
        fail(f"timepoint mismatch: {r.library_id} {r.raw_sample_id}")

expected_new_controls = {
    "ТУ1_010", "СМ1_021", "СА1_024", "КА1_025", "ПЛ1_032",
    "АЕ1_035", "СА1_046", "АЖ1_050", "АА1_051",
}
observed_new_controls = set(new.loc[new.is_control.astype(str).str.lower().isin(["true", "1"]), "patient_id"])
if observed_new_controls != expected_new_controls:
    fail(f"new control mapping changed: {sorted(observed_new_controls)}")

# Explicit legacy rules from the colleague comment.
legacy = legacy.set_index("sequencing_index")
for udp in [f"UDP{i:04d}" for i in range(97, 106)]:
    if str(legacy.loc[udp, "disease_status"]) != "measles" or int(float(legacy.loc[udp, "analysis_timepoint"])) != 1:
        fail(f"legacy case rule broken for {udp}")
for udp in [f"UDP{i:04d}" for i in range(106, 113)]:
    if str(legacy.loc[udp, "disease_status"]) != "healthy":
        fail(f"legacy CGIMU control rule broken for {udp}")
if str(legacy.loc["UDP0108", "include"]).lower() not in {"false", "0"}:
    fail("UDP0108 must remain excluded from biological inference")

env = (ROOT / "environment.yml").read_text(encoding="utf-8")
for token in ["milaboratories", "mixcr=4.7.0", "matplotlib-venn=1.1", "python=3.12.13"]:
    if token not in env:
        fail(f"environment missing reproducibility token: {token}")

prep = (ROOT / "scripts/preprocess_human_trb.py").read_text(encoding="utf-8")
for token in ["EXPECTED_NEW_LANES", "--force-overwrite", "fastq_audit.csv", "UDP0107"]:
    if token not in prep:
        fail(f"preprocessing contract missing token: {token}")

analysis = (ROOT / "scripts/analysis_01_prepare.py").read_text(encoding="utf-8")
for token in ["analysis_timepoint", "multivariate_hypergeometric", "source_patient_id"]:
    if token not in analysis:
        fail(f"analysis preparation contract missing token: {token}")

longitudinal = (ROOT / "scripts/analysis_05_longitudinal_venn.py").read_text(encoding="utf-8")
for token in ["style_venn", "aav_sets", "trbv_group_statistics.csv", "adult vs child"]:
    if token not in longitudinal:
        fail(f"longitudinal contract missing token: {token}")

print("Repository validation passed.")
