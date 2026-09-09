#!/usr/bin/env python3
"""Validate the published full threshold grid without redrawing figures."""
from pathlib import Path
import os
import pandas as pd

ROOT = Path(os.environ.get("MEASLES_REPO", Path(__file__).resolve().parents[1])).resolve()
RESULTS = ROOT / "results"
FINAL = RESULTS / "final"
SUPP = RESULTS / "supplementary_tables"
FIG = ROOT / "figures"

THRESHOLD_PCT = [0.05, 0.10] + [round(x / 10, 2) for x in range(2, 16)]
source = SUPP / "03_threshold_sensitivity.csv"
candidate_table = FINAL / "candidate_clones.csv"
for path in [source, candidate_table]:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)

sens = pd.read_csv(source)
observed = {round(float(x), 8) for x in sens.threshold_percent.dropna().unique()}
if observed != set(THRESHOLD_PCT):
    raise RuntimeError(f"Unexpected extended threshold set: {sorted(observed)}")
coverage = sens.groupby(["branch", "threshold_percent"]).patient.nunique()
if (coverage != 22).any():
    raise RuntimeError(f"Incomplete branch-threshold coverage: {coverage[coverage != 22].to_dict()}")

for path in [
    FIG / "rarefied" / "02_sampling_and_sensitivity.png",
    FIG / "simpson_raw" / "02_shared_and_sensitivity.png",
    FIG / "validation" / "03_threshold_selection.png",
    FIG / "validation" / "04_downsampling_iterations.png",
]:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)

for cohort in ["adults", "children"]:
    for relative in [
        "rarefied/02_sampling_and_sensitivity.png",
        "simpson_raw/02_shared_and_sensitivity.png",
        "validation/03_threshold_selection.png",
        "validation/04_downsampling_iterations.png",
    ]:
        path = FIG / cohort / relative
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

print("Published threshold and iteration outputs validated")
print("Thresholds percent:", THRESHOLD_PCT)
