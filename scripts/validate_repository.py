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
    "README.md", "RUN_GUIDE.md", "environment.yml",
    "config/sample_manifest.csv", "config/clinical_metadata_sequenced.csv",
    "config/metadata_known_issues.csv", "config/SOURCE_NOTES.md",
    "notebooks/00_prepare_and_mixcr.ipynb", "notebooks/01_measles_analysis.ipynb",
    "scripts/preprocess_human_trb.py", "scripts/mixcr_smoke_test.py",
    "scripts/validate_clonotype_exports.py",
    "scripts/analysis_01_prepare.py", "scripts/analysis_02_antigens.py",
    "scripts/analysis_03_figures.py", "scripts/analysis_04_stats.py",
    "scripts/analysis_05_longitudinal_venn.py", "scripts/analysis_06_clinical_summary.py",
    "scripts/run_pipeline.py", "scripts/run_slurm.sh", "run_all_from_scratch.sh",
]
missing = [p for p in required if not (ROOT / p).exists()]
if missing:
    fail(f"missing required files: {missing}")

for rel in ["notebooks/00_prepare_and_mixcr.ipynb", "notebooks/01_measles_analysis.ipynb"]:
    nb = json.loads((ROOT / rel).read_text(encoding="utf-8"))
    ids = [c.get("id") for c in nb.get("cells", [])]
    if any(x is None for x in ids) or len(ids) != len(set(ids)):
        fail(f"invalid/duplicate notebook cell IDs: {rel}")
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code" and (cell.get("outputs") or cell.get("execution_count") is not None):
            fail(f"source notebook contains outputs/execution count: {rel}")

for path in (ROOT / "scripts").glob("*.py"):
    try:
        ast.parse(path.read_text(encoding="utf-8"))
    except SyntaxError as exc:
        fail(f"Python syntax error in {path.name}: {exc}")

manifest = pd.read_csv(ROOT / "config/sample_manifest.csv")
required_cols = {
    "library_id", "sequencing_index", "batch", "raw_sample_id", "patient_id",
    "analysis_timepoint", "clinical_timepoint", "disease_status", "age_group",
    "is_control", "include", "fastq_key",
}
if not required_cols.issubset(manifest.columns):
    fail(f"manifest columns missing: {sorted(required_cols - set(manifest.columns))}")

new = manifest[manifest.batch.eq("new")].copy()
legacy = manifest[manifest.batch.eq("legacy")].copy()
if len(new) != 96:
    fail(f"expected 96 new sequencing libraries, found {len(new)}")
if len(legacy) != 16:
    fail(f"expected 16 legacy libraries, found {len(legacy)}")
if new.library_id.duplicated().any() or new.sequencing_index.duplicated().any():
    fail("new library IDs and sequencing indices must be unique")

# Raw labels retain their original suffix for traceability; clinical stage is separate.
for r in new.itertuples(index=False):
    if not str(r.raw_sample_id).startswith(str(r.sequencing_index) + "_"):
        fail(f"sequencing prefix mismatch: {r.library_id} {r.raw_sample_id}")
    m = re.search(r"\.([1-4])(?:[A-Za-zА-Яа-яЁё])?$", str(r.raw_sample_id))
    if not m or int(m.group(1)) != int(float(r.analysis_timepoint)):
        fail(f"raw suffix / analysis_timepoint mismatch: {r.library_id} {r.raw_sample_id}")

new_by_lib = new.set_index("library_id")
expected_patient = {
    "new_21": "ОЛ1_018",
    "new_22": "СА1_019",
    "new_23": "МС1_020",
}
for lib, patient in expected_patient.items():
    if str(new_by_lib.loc[lib, "patient_id"]) != patient:
        fail(f"clinical patient reconciliation changed for {lib}: expected {patient}")

expected_phase = {
    "new_15": "T1", "new_19": "T2", "new_72": "T3",   # СА1_013
    "new_38": "T1", "new_39": "T2", "new_43": "T2",   # ММ1_033
    "new_51": "T0", "new_56": "T1", "new_62": "T2",   # ЭГ1_042
    "new_10": "T2", "new_12": "T3",                       # ПМ1_009
    "new_46": "T3", "new_75": "T3", "new_57": "T2",
}
for lib, phase in expected_phase.items():
    if str(new_by_lib.loc[lib, "clinical_timepoint"]) != phase:
        fail(f"clinical timepoint changed for {lib}: expected {phase}")

included_new = new[new.include.astype(str).str.lower().isin(["true", "1", "yes"])].copy()
if included_new.age_group.astype(str).eq("unknown").any():
    bad = sorted(included_new.loc[included_new.age_group.astype(str).eq("unknown"), "library_id"])
    fail(f"included new libraries must have adult/child metadata: {bad}")
if included_new.clinical_timepoint.fillna("").astype(str).eq("").any():
    bad = sorted(included_new.loc[included_new.clinical_timepoint.fillna("").astype(str).eq(""), "library_id"])
    fail(f"included new libraries must have clinical_timepoint: {bad}")
if not set(included_new.clinical_timepoint.astype(str)).issubset({"T0", "T1", "T2", "T3"}):
    fail("included new clinical_timepoint values must be T0/T1/T2/T3")

expected_excluded_new = {"new_24", "new_39", "new_88", "new_111"}
observed_excluded_new = set(new.loc[~new.include.astype(str).str.lower().isin(["true", "1", "yes"]), "library_id"])
if observed_excluded_new != expected_excluded_new:
    fail(f"new-library exclusion set changed: {sorted(observed_excluded_new)}")

# The only unambiguous included new T3 libraries in the current workbook.
observed_t3 = set(included_new.loc[included_new.clinical_timepoint.eq("T3"), "library_id"])
expected_t3 = {"new_12", "new_46", "new_72", "new_75"}
if observed_t3 != expected_t3:
    fail(f"clinical T3 set changed: {sorted(observed_t3)}")

expected_new_controls = {
    "ТУ1_010", "СМ1_021", "СА1_024", "КА1_025", "ПЛ1_032",
    "АЕ1_035", "СА1_046", "АЖ1_050", "АА1_051",
    "ВА1_060", "ЛВ1_061", "КП1_064",
}
observed_new_controls = set(
    included_new.loc[included_new.is_control.astype(str).str.lower().isin(["true", "1", "yes"]), "patient_id"]
)
if observed_new_controls != expected_new_controls:
    fail(f"new control mapping changed: {sorted(observed_new_controls)}")

clinical = pd.read_csv(ROOT / "config/clinical_metadata_sequenced.csv").set_index("patient_id")
for r in included_new.itertuples(index=False):
    if r.patient_id not in clinical.index:
        fail(f"included patient missing from clinical metadata table: {r.patient_id}")
    c = clinical.loc[r.patient_id]
    if str(c["disease_status"]) != str(r.disease_status):
        fail(f"disease-status mismatch for {r.patient_id}")
    if str(c["age_group"]) != str(r.age_group):
        fail(f"age-group mismatch for {r.patient_id}")

legacy_idx = legacy.set_index("sequencing_index")
if set(legacy_idx.age_group.astype(str)) != {"adult"}:
    fail("all legacy libraries must be adults")
if legacy_idx.patient_id.duplicated().any():
    fail("legacy libraries must represent distinct people")
for udp in [f"UDP{i:04d}" for i in range(97, 113)]:
    if str(legacy_idx.loc[udp, "patient_id"]) != udp:
        fail(f"legacy patient {udp} must remain independent")
for udp in [f"UDP{i:04d}" for i in range(97, 106)]:
    row = legacy_idx.loc[udp]
    if str(row.disease_status) != "measles" or str(row.clinical_timepoint) != "legacy_acute":
        fail(f"legacy case rule broken for {udp}")
for udp in [f"UDP{i:04d}" for i in range(106, 113)]:
    row = legacy_idx.loc[udp]
    if str(row.disease_status) != "healthy" or str(row.clinical_timepoint) != "control":
        fail(f"legacy control rule broken for {udp}")
if str(legacy_idx.loc["UDP0108", "include"]).lower() not in {"false", "0"}:
    fail("UDP0108 must remain excluded from biological inference")

env = (ROOT / "environment.yml").read_text(encoding="utf-8")
for token in ["python=3.12.13", "matplotlib-venn=1.1", "milaboratories", "mixcr=4.7.0"]:
    if token not in env:
        fail(f"environment missing reproducibility token: {token}")

print("Repository validation passed.")
