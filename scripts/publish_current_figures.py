#!/usr/bin/env python3
"""Publish authoritative QC-filtered figures into the top-level figures/ directory.

The analysis scripts intentionally write intermediate/history outputs into
figures/corrected* and authoritative final outputs into figures/current.
This publisher makes the top-level figures/ directory user-facing and removes
stale top-level PNG/PDF files from older iterations.
"""
from __future__ import annotations

from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
CURRENT = FIG / "current"
CURRENT_PATIENTS = CURRENT / "patients"
ROOT_PATIENTS = FIG / "patients"

marker = CURRENT / "GENERATED_CURRENT.txt"
if not marker.exists() or marker.stat().st_size == 0:
    raise RuntimeError(
        "Authoritative figures/current outputs are missing. "
        "Run scripts/analysis_08_final_reporting.py successfully before publishing."
    )

current_main = sorted([p for p in CURRENT.iterdir() if p.suffix.lower() in {".png", ".pdf"}])
current_patient = sorted([p for p in CURRENT_PATIENTS.glob("*") if p.suffix.lower() in {".png", ".pdf"}])
if not current_main:
    raise RuntimeError("No authoritative PNG/PDF files found in figures/current.")
if not current_patient:
    raise RuntimeError("No authoritative patient PNG/PDF files found in figures/current/patients.")

# Remove only published top-level figure files from earlier iterations.
for p in FIG.iterdir():
    if p.is_file() and p.suffix.lower() in {".png", ".pdf"}:
        p.unlink()

if ROOT_PATIENTS.exists():
    shutil.rmtree(ROOT_PATIENTS)
ROOT_PATIENTS.mkdir(parents=True, exist_ok=True)

for src in current_main:
    shutil.copy2(src, FIG / src.name)
for src in current_patient:
    shutil.copy2(src, ROOT_PATIENTS / src.name)

publish_note = FIG / "CURRENT_OUTPUTS.txt"
publish_note.write_text(
    marker.read_text(encoding="utf-8")
    + "\n\nTop-level figures/ was refreshed from figures/current by publish_current_figures.py.\n"
    + f"Main figure files published: {len(current_main)}\n"
    + f"Patient figure files published: {len(current_patient)}\n",
    encoding="utf-8",
)

print(f"Published {len(current_main)} main files to {FIG}")
print(f"Published {len(current_patient)} patient files to {ROOT_PATIENTS}")
print(f"Marker: {publish_note}")
