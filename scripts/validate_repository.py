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
    "README.md", "RUN_GUIDE.md", "VALIDATION.md", "TZ_COVERAGE.md", "environment.yml",
    "config/sample_manifest.csv", "config/analysis_metadata.csv",
    "config/clinical_metadata_sequenced.csv", "config/metadata_known_issues.csv", "config/SOURCE_NOTES.md",
    "notebooks/00_prepare_and_mixcr.ipynb", "notebooks/01_measles_analysis.ipynb",
    "scripts/analysis_01_prepare.py", "scripts/analysis_02_antigens.py",
    "scripts/analysis_06_clinical_summary.py", "scripts/analysis_07_longitudinal_expanded.py",
    "scripts/analysis_08_final_reporting.py", "scripts/analysis_09_patient_panels.py",
    "scripts/analysis_10_clone_fate_v2.py", "scripts/analysis_11_clone_enrichment.py",
    "scripts/analysis_12_exact_rarefaction.py", "scripts/analysis_14_simpson_no_downsampling.py",
    "scripts/analysis_15_tz_final_reporting.py", "scripts/analysis_16_extended_thresholds.py",
    "scripts/analysis_17_apply_extended_thresholds.py", "scripts/plot_style.py",
    "scripts/run_analysis_refresh_slurm.sh",
]
missing = [p for p in required if not (ROOT / p).exists()]
if missing:
    fail(f"missing required files: {missing}")
if (ROOT / "scripts/analysis_07_age6_reporting.py").exists():
    fail("obsolete age 6 reporting script must be removed")
if (ROOT / "scripts/analysis_08_final_reporting_safe.py").exists():
    fail("obsolete compatibility wrapper must be removed")

for path in (ROOT / "scripts").glob("*.py"):
    try:
        ast.parse(path.read_text(encoding="utf-8"))
    except SyntaxError as exc:
        fail(f"Python syntax error in {path.name}: {exc}")

for rel in ["notebooks/00_prepare_and_mixcr.ipynb", "notebooks/01_measles_analysis.ipynb"]:
    nb = json.loads((ROOT / rel).read_text(encoding="utf-8"))
    ids = [c.get("id") for c in nb.get("cells", [])]
    if any(x is None for x in ids) or len(ids) != len(set(ids)):
        fail(f"invalid notebook cell IDs in {rel}")
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code" and (cell.get("outputs") or cell.get("execution_count") is not None):
            fail(f"source notebook contains outputs in {rel}")

manifest = pd.read_csv(ROOT / "config/sample_manifest.csv")
meta = pd.read_csv(ROOT / "config/analysis_metadata.csv")
for col in ["include", "is_control"]:
    manifest[col] = manifest[col].astype(str).str.lower().isin(["true", "1", "yes"])
meta["passes_qc_4400"] = meta["passes_qc_4400"].astype(str).str.lower().isin(["true", "1", "yes"])
meta["productive_trb_clonotypes"] = pd.to_numeric(meta.productive_trb_clonotypes, errors="coerce")

if len(manifest) != 112 or len(meta) != 112:
    fail(f"expected 112 libraries, found manifest {len(manifest)} metadata {len(meta)}")
if len(manifest[manifest.batch.eq("new")]) != 96 or len(manifest[manifest.batch.eq("legacy")]) != 16:
    fail("expected 96 new and 16 legacy libraries")
if manifest.library_id.duplicated().any() or meta.library_id.duplicated().any():
    fail("library IDs must be unique")

joined = manifest.merge(meta, on="library_id", suffixes=("_manifest", ""), how="outer", indicator=True)
if not joined._merge.eq("both").all():
    fail("sample manifest and analysis metadata have different library sets")

if manifest.clinical_timepoint.astype(str).eq("T0").any():
    fail("T0 does not exist in this study")
new = manifest[manifest.batch.eq("new")]
for r in new.itertuples(index=False):
    m = re.search(r"\.([1-3])(?:[A-Za-zА-Яа-яЁё])?$", str(r.raw_sample_id))
    if not m:
        fail(f"cannot parse new library suffix for {r.library_id}")
    suffix = int(m.group(1))
    if int(float(r.analysis_timepoint)) != suffix:
        fail(f"analysis time point mismatch for {r.library_id}")
    expected = "CONTROL" if bool(r.is_control) else f"T{suffix}"
    if str(r.clinical_timepoint) != expected:
        fail(f"clinical time point mismatch for {r.library_id}: expected {expected}")

if not manifest.loc[manifest.is_control, "clinical_timepoint"].eq("CONTROL").all():
    fail("all healthy controls must use CONTROL")
case_points = set(manifest.loc[manifest.disease_status.eq("measles"), "clinical_timepoint"].astype(str))
if not case_points.issubset({"T1", "T2", "T3"}):
    fail(f"unexpected measles time points: {sorted(case_points)}")

if set(meta.age_group) - {"child", "adult", "unknown"}:
    fail("age_group must be child adult or unknown")
for r in meta.itertuples(index=False):
    text = str(r.age_text).lower()
    m = re.search(r"(\d+)\s*(?:год|года|лет)", text)
    if m:
        years = int(m.group(1))
        expected = "child" if years < 18 else "adult"
        if r.age_group != expected:
            fail(f"age rule mismatch for {r.library_id}: {r.age_text} expected {expected}")

if not (meta.loc[meta.passes_qc_4400, "productive_trb_clonotypes"] >= 4400).all():
    fail("a library marked as passing QC has fewer than 4400 productive TRB clonotypes")
manifest_include = manifest.set_index("library_id")["include"]
final = meta[meta.passes_qc_4400 & meta.disease_status.isin(["healthy", "measles"])].copy()
final = final[manifest_include.loc[final.library_id].to_numpy()]
if final.patient_id.nunique() != 78:
    fail(f"expected 78 unique final patients, found {final.patient_id.nunique()}")
patients = final.drop_duplicates("patient_id")
counts = patients.age_group.value_counts().to_dict()
if counts.get("child", 0) != 33 or counts.get("adult", 0) != 45:
    fail(f"expected 33 children and 45 adults, found {counts}")

obs = final.drop_duplicates(["patient_id", "clinical_timepoint"])
expected_counts = {
    ("child", "CONTROL"): 4, ("child", "T1"): 27, ("child", "T2"): 8,
    ("adult", "CONTROL"): 13, ("adult", "T1"): 29, ("adult", "T2"): 19,
}
for key, expected in expected_counts.items():
    got = len(obs[(obs.age_group == key[0]) & (obs.clinical_timepoint == key[1])])
    if got != expected:
        fail(f"unexpected cohort count for {key}: expected {expected}, found {got}")

paired = []
for patient, z in obs[obs.disease_status.eq("measles")].groupby("patient_id"):
    if {"T1", "T2"}.issubset(set(z.clinical_timepoint.astype(str))):
        paired.append(str(patient))
if len(paired) != 22:
    fail(f"expected 22 paired T1-T2 patients, found {len(paired)}")
paired_meta = obs[obs.patient_id.astype(str).isin(paired)]
paired_counts = paired_meta.drop_duplicates("patient_id").age_group.value_counts().to_dict()
if paired_counts.get("child", 0) != 6 or paired_counts.get("adult", 0) != 16:
    fail(f"expected paired cohort 6 children and 16 adults, found {paired_counts}")

rep = meta[meta.replicate_group.fillna("").eq("БА1_071_T2")]
if set(rep.library_id) != {"new_106", "new_111"} or not rep.passes_qc_4400.all():
    fail("БА1_071 T2 replicate group must contain accepted new_106 and new_111")
mr = manifest.set_index("library_id")
if not bool(mr.loc["new_111", "include"]) or str(mr.loc["new_111", "clinical_timepoint"]) != "T2":
    fail("new_111 must be included as T2")

excluded = set(manifest.loc[~manifest.include, "library_id"])
if excluded != {"new_24", "new_39", "legacy_UDP0108"}:
    fail(f"unexpected config exclusion set: {sorted(excluded)}")

legacy_expected = {
    "legacy_UDP0097": ("МР1_001", "T1", "measles", "child"),
    "legacy_UDP0098": ("БА1_002", "T1", "measles", "child"),
    "legacy_UDP0099": ("ГА1_003", "T1", "measles", "child"),
    "legacy_UDP0100": ("ХИ1_004", "T1", "measles", "child"),
    "legacy_UDP0101": ("ИМ1_005", "T1", "measles", "child"),
    "legacy_UDP0102": ("ПИ1_006", "T1", "measles", "adult"),
    "legacy_UDP0103": ("СС1_007", "T1", "measles", "child"),
    "legacy_UDP0104": ("ХИ1_004", "T2", "measles", "child"),
    "legacy_UDP0105": ("ШЕ1_008", "T1", "measles", "child"),
}
for lib, expected in legacy_expected.items():
    r = mr.loc[lib]
    got = (str(r.patient_id), str(r.clinical_timepoint), str(r.disease_status), str(r.age_group))
    if got != expected:
        fail(f"legacy mapping mismatch for {lib}: {got}")
for i in range(106, 113):
    lib = f"legacy_UDP{i:04d}"
    r = mr.loc[lib]
    if str(r.patient_id) != f"CGIMU_N{i-105}" or str(r.clinical_timepoint) != "CONTROL" or str(r.disease_status) != "healthy":
        fail(f"legacy control mapping mismatch for {lib}")

threshold_grid_token = "THRESHOLD_PCT = [0.05, 0.10] + [round(x / 10, 2) for x in range(2, 16)]"
threshold_scripts = [
    "scripts/analysis_12_exact_rarefaction.py",
    "scripts/analysis_14_simpson_no_downsampling.py",
    "scripts/analysis_16_extended_thresholds.py",
    "scripts/analysis_17_apply_extended_thresholds.py",
]
for rel in threshold_scripts:
    text = (ROOT / rel).read_text(encoding="utf-8")
    if threshold_grid_token not in text:
        fail(f"extended threshold grid missing in {rel}")

contracts = {
    "scripts/analysis_01_prepare.py": ["QC_MIN_PRODUCTIVE_CLONOTYPES = 4400", "EXPANDED_FRACTION = 0.001", "БА1_071_T2"],
    "scripts/analysis_07_longitudinal_expanded.py": ["expanded_jaccard", "new_T2_expansions", "expected = 6 if age == \"child\" else 16"],
    "scripts/analysis_09_patient_panels.py": ["Expanded clonotype overlap", "EXPANDED_FRACTION = 0.001"],
    "scripts/analysis_10_clone_fate_v2.py": [
        "PRIMARY = 0.001", "clone_id_nt", "control_prevalence", "repertoire_identity_morisita_matrix.csv",
    ],
    "scripts/analysis_11_clone_enrichment.py": [
        "sampling_candidate_class", "sampling-robust not detected", "strong contraction",
        "sm.GEE", "log10_T1_fraction", "TRBV", "TRBJ", "CDR3_3mer",
    ],
    "scripts/analysis_12_exact_rarefaction.py": [
        "pair_depth = min", "multivariate_hypergeometric", "N_ITER < 500", "PRIMARY = 0.001",
        "P_detection_if_unchanged_empirical", "clone_sampling_null.csv", "pairwise_threshold_sensitivity.csv",
        "downsampling_iteration_convergence.csv",
    ],
    "scripts/analysis_14_simpson_no_downsampling.py": [
        "PRIMARY = 0.001", "simpson_D", "simpson_clonality", "normalized_richness",
        "simpson_raw_threshold_sensitivity.csv",
    ],
    "scripts/analysis_15_tz_final_reporting.py": [
        "FINAL = RESULTS / \"final\"", "candidate_class_rarefied", "candidate_class_simpson_raw",
        "01_shared_and_fate.png", "03_candidate_dumbbell_and_trbv.png",
        "03_threshold_selection.png", "04_downsampling_iterations.png",
    ],
    "scripts/analysis_16_extended_thresholds.py": [
        "pairwise_threshold_sensitivity.csv", "simpson_raw_threshold_sensitivity.csv",
        "pairwise_rarefied", "simpson_no_downsampling", "expanded_threshold_sensitivity_extended.csv",
        "does not repeat rarefaction",
    ],
    "scripts/analysis_17_apply_extended_thresholds.py": [
        "03_threshold_sensitivity.csv", "03_threshold_selection.png",
        "Published threshold and iteration outputs validated",
    ],
    "scripts/run_analysis_refresh_slurm.sh": [
        "PAIR_ITERATIONS", "MEASLES_RAREFACTION_ITERATIONS=1", "analysis_12_exact_rarefaction.py",
        "analysis_11_clone_enrichment.py", "analysis_14_simpson_no_downsampling.py",
        "analysis_16_extended_thresholds.py", "analysis_15_tz_final_reporting.py",
        "analysis_17_apply_extended_thresholds.py",
        "Final sensitivity thresholds: 0.05 0.10 and 0.20 to 1.50 percent in 0.10 percentage-point steps",
    ],
}
for rel, tokens in contracts.items():
    text = (ROOT / rel).read_text(encoding="utf-8")
    for token in tokens:
        if token not in text:
            fail(f"analysis contract token missing in {rel}: {token}")

print("Repository validation passed")
print("Authoritative metadata: config/analysis_metadata.csv")
print("QC threshold: 4400 productive TRB clonotypes")
print("Age rule: child younger than 18 years, adult 18 years and older")
print("Controls use CONTROL only")
print("new_111 is included as a T2 biological replicate of new_106")
print("Physical clone tracking: CDR3 nucleotide plus TRBV plus TRBJ")
print("Primary expanded threshold: 0.10 percent")
print("Sensitivity thresholds: 0.05, 0.10, then 0.20 to 1.50 percent in 0.10 percentage-point steps")
print("Pair-specific repeated rarefaction: at least 500 iterations at each pair's own common depth")
print("Same-depth internal null provides mean exact P_detection_if_unchanged")
print("No-downsampling Simpson branch uses the same expansion-threshold grid")
