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
    "config/metadata_known_issues.csv", "config/SOURCE_NOTES.md",
    "notebooks/00_prepare_and_mixcr.ipynb", "notebooks/01_measles_analysis.ipynb",
    "scripts/preprocess_human_trb.py", "scripts/mixcr_smoke_test.py",
    "scripts/analysis_01_prepare.py", "scripts/analysis_02_antigens.py",
    "scripts/analysis_03_figures.py", "scripts/analysis_04_stats.py",
    "scripts/analysis_05_longitudinal_venn.py", "scripts/run_pipeline.py",
    "scripts/run_slurm.sh", "run_all_from_scratch.sh",
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
new = manifest[manifest.batch.eq("new")].copy()
legacy = manifest[manifest.batch.eq("legacy")].copy()
if len(new) != 96:
    fail(f"expected 96 new study-table libraries, found {len(new)}")
if len(legacy) != 16:
    fail(f"expected 16 legacy UDP libraries, found {len(legacy)}")
if new.library_id.duplicated().any() or new.sequencing_index.duplicated().any():
    fail("new sequencing library IDs/indices must be unique")

if set(new.loc[new.patient_id.str.startswith("СА1_", na=False), "patient_id"]) < {"СА1_013", "СА1_017"}:
    fail("patient IDs СА1_013 and СА1_017 are not kept as distinct people")

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

legacy = legacy.set_index("sequencing_index")
if set(legacy["age_group"].astype(str)) != {"adult"}:
    fail("all legacy UDP0097-UDP0112 samples must be adults")
if legacy["patient_id"].duplicated().any():
    fail("legacy libraries must represent distinct people")
for udp in [f"UDP{i:04d}" for i in range(97, 113)]:
    if str(legacy.loc[udp, "patient_id"]) != udp:
        fail(f"legacy patient {udp} must remain independent from the new cohort")
for udp in [f"UDP{i:04d}" for i in range(97, 106)]:
    row = legacy.loc[udp]
    if str(row["disease_status"]) != "measles" or str(row["age_group"]) != "adult":
        fail(f"legacy adult case rule broken for {udp}")
    if int(float(row["analysis_timepoint"])) != 1 or str(row["is_control"]).lower() not in {"false", "0"}:
        fail(f"legacy adult case timepoint/control rule broken for {udp}")
for udp in [f"UDP{i:04d}" for i in range(106, 113)]:
    row = legacy.loc[udp]
    if str(row["disease_status"]) != "healthy" or str(row["age_group"]) != "adult":
        fail(f"legacy adult control rule broken for {udp}")
    if str(row["is_control"]).lower() not in {"true", "1"}:
        fail(f"legacy control flag broken for {udp}")
if str(legacy.loc["UDP0108", "include"]).lower() not in {"false", "0"}:
    fail("UDP0108 must remain excluded from biological inference")

env = (ROOT / "environment.yml").read_text(encoding="utf-8")
for token in ["matplotlib-venn=1.1", "python=3.12.13"]:
    if token not in env:
        fail(f"environment missing reproducibility token: {token}")
if "mixcr=" in env.lower() or "milaboratories" in env.lower():
    fail("MiXCR must stay external to conda, matching the working mouse-project setup")

prep = (ROOT / "scripts/preprocess_human_trb.py").read_text(encoding="utf-8")
for token in [
    "EXPECTED_NEW_LANES", "~/soft/mixcr/mixcr", "-Xmx", "\"-f\"",
    "--threads", "fastq_audit.csv", "_1.fq.gz", "UDP0107", "L{{{{n}}}}",
]:
    if token not in prep:
        fail(f"preprocessing contract missing token: {token}")

smoke = (ROOT / "scripts/mixcr_smoke_test.py").read_text(encoding="utf-8")
for token in ["--limit-input", "~/soft/mixcr/mixcr", "-Xmx4g", "MiXCR smoke test passed"]:
    if token not in smoke:
        fail(f"MiXCR smoke-test contract missing token: {token}")

analysis = (ROOT / "scripts/analysis_01_prepare.py").read_text(encoding="utf-8")
for token in ["analysis_timepoint", "multivariate_hypergeometric", "source_patient_id"]:
    if token not in analysis:
        fail(f"analysis preparation contract missing token: {token}")

longitudinal = (ROOT / "scripts/analysis_05_longitudinal_venn.py").read_text(encoding="utf-8")
for token in ["style_venn", "aav_sets", "trbv_group_statistics.csv", "adult vs child"]:
    if token not in longitudinal:
        fail(f"longitudinal contract missing token: {token}")

print("Repository validation passed.")
