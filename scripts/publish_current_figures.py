#!/usr/bin/env python3
"""Publish only final group and patient figures"""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
CURRENT = FIG / "current"
CURRENT_PATIENTS = CURRENT / "patients"
GROUP = FIG / "group"
PATIENTS = FIG / "patients"

marker = CURRENT / "GENERATED_CURRENT.txt"
if not marker.exists():
    raise RuntimeError("Final figure marker is missing")

for d in [GROUP, PATIENTS]:
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)

main_files = sorted(p for p in CURRENT.iterdir() if p.suffix.lower() in {".png", ".pdf"})
patient_files = sorted(p for p in CURRENT_PATIENTS.iterdir() if p.suffix.lower() in {".png", ".pdf"})
if not main_files:
    raise RuntimeError("No final group figures were generated")
if not patient_files:
    raise RuntimeError("No final patient figures were generated")

for src in main_files:
    shutil.copy2(src, GROUP / src.name)
for src in patient_files:
    shutil.copy2(src, PATIENTS / src.name)

index = FIG / "FIGURE_INDEX.txt"
lines = [
    "FINAL FIGURES",
    "",
    marker.read_text(encoding="utf-8").strip(),
    "",
    "GROUP ANALYSIS",
]
lines.extend(f"group/{p.name}" for p in main_files)
lines.extend(["", "PATIENT ANALYSIS"])
lines.extend(f"patients/{p.name}" for p in patient_files)
lines.extend([
    "",
    "Analysis logic",
    "Cross sectional figures use all accepted T1 observations, all accepted T2 observations and healthy controls",
    "Paired figures use only patients with both accepted T1 and T2 observations",
    "Patient overlap is calculated only for expanded aaV clonotypes",
    "Expanded means at least 0.10 percent of productive TRB reads",
    "БА1_071 T2 is one patient observation formed by averaging normalized frequencies from new_106 and new_111",
])
index.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"Published {len(main_files)} group figure files")
print(f"Published {len(patient_files)} patient figure files")
