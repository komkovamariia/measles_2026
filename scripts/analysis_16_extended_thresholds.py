#!/usr/bin/env python3
"""Consolidate the full expansion-threshold grid from both completed branches.

The expensive calculations are already performed by:
- analysis_12_exact_rarefaction.py for pair-specific repeated rarefaction
- analysis_14_simpson_no_downsampling.py for observed frequencies without downsampling

This script does not repeat rarefaction. It only validates and concatenates the two
full sensitivity tables so the final publisher can keep a single authoritative table.

Threshold grid in percent:
0.05, 0.10, then 0.20 to 1.50 in 0.10 percentage-point steps.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
THRESHOLD_PCT = [0.05, 0.10] + [round(x / 10, 2) for x in range(2, 16)]
EXPECTED = set(THRESHOLD_PCT)

rare_path = CURRENT / "pairwise_threshold_sensitivity.csv"
raw_path = CURRENT / "simpson_raw_threshold_sensitivity.csv"
for p in [rare_path, raw_path]:
    if not p.exists() or p.stat().st_size == 0:
        raise FileNotFoundError(p)

rare = pd.read_csv(rare_path)
raw = pd.read_csv(raw_path)
rare["branch"] = "pairwise_rarefied"
raw["branch"] = "simpson_no_downsampling"

for name, df in [("pairwise_rarefied", rare), ("simpson_no_downsampling", raw)]:
    observed = {round(float(x), 8) for x in df.threshold_percent.dropna().unique()}
    if observed != EXPECTED:
        raise RuntimeError(f"{name} threshold grid mismatch: {sorted(observed)}")
    if df.patient.nunique() != 22:
        raise RuntimeError(f"{name} has {df.patient.nunique()} patients instead of 22")

# Harmonize schemas. Rarefied summaries have median and interval columns; the raw
# branch has direct point estimates. Duplicate raw values into the median columns so
# downstream plotting can use one schema without altering the underlying numbers.
for base in [
    "expanded_jaccard", "retained_count_percent", "below_count_percent",
    "not_detected_count_percent", "retained_mass_percent", "below_mass_percent",
    "not_detected_mass_percent",
]:
    if base in raw.columns and f"{base}_median" not in raw.columns:
        raw[f"{base}_median"] = raw[base]
        raw[f"{base}_q025"] = raw[base]
        raw[f"{base}_q975"] = raw[base]

out = pd.concat([rare, raw], ignore_index=True, sort=False)
out = out.sort_values(["branch", "patient", "threshold_percent"]).reset_index(drop=True)
observed = {round(float(x), 8) for x in out.threshold_percent.dropna().unique()}
if observed != EXPECTED:
    raise RuntimeError(f"Combined threshold grid mismatch: {sorted(observed)}")
if out.groupby("branch").patient.nunique().to_dict() != {
    "pairwise_rarefied": 22,
    "simpson_no_downsampling": 22,
}:
    raise RuntimeError("Combined threshold table does not contain both complete branches")

out.to_csv(RESULTS / "expanded_threshold_sensitivity_extended.csv", index=False)
print("Full threshold grid consolidated without repeating rarefaction")
print("Thresholds percent:", THRESHOLD_PCT)
print("Rows:", len(out))
