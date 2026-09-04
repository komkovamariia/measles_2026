#!/usr/bin/env python3
"""Publish only the final user-facing figures.

Final layout:
figures/
  group/      final cohort-level figures
  patients/   final patient-level panels
  FIGURE_INDEX.txt

All stale top-level PNG/PDF files and intermediate corrected/current figure
folders are removed during publication.
"""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
CURRENT = FIG / "current"
CURRENT_PATIENTS = CURRENT / "patients"
GROUP = FIG / "group"
PATIENTS = FIG / "patients"

marker = CURRENT / "GENERATED_CURRENT.txt"
if not marker.exists() or marker.stat().st_size == 0:
    raise RuntimeError(
        "Authoritative figures/current outputs are missing. "
        "Run final reporting successfully before publishing."
    )

main_files = sorted(p for p in CURRENT.iterdir() if p.suffix.lower() in {".png", ".pdf"})
patient_files = sorted(p for p in CURRENT_PATIENTS.iterdir() if p.suffix.lower() in {".png", ".pdf"})
if not main_files:
    raise RuntimeError("No final group PNG/PDF files found in figures/current.")
if not patient_files:
    raise RuntimeError("No final patient PNG/PDF files found in figures/current/patients.")

# Clean stale published outputs.
for d in [GROUP, PATIENTS]:
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)

for p in FIG.iterdir():
    if p.is_file() and (p.suffix.lower() in {".png", ".pdf"} or p.name in {"CURRENT_OUTPUTS.txt", "FIGURE_INDEX.txt"}):
        p.unlink()

for src in main_files:
    shutil.copy2(src, GROUP / src.name)
for src in patient_files:
    shutil.copy2(src, PATIENTS / src.name)

index = FIG / "FIGURE_INDEX.txt"
lines = [
    "FINAL USER-FACING FIGURES",
    "",
    marker.read_text(encoding="utf-8").strip(),
    "",
    f"Group figures: {len(main_files)}",
]
lines.extend(f"  group/{p.name}" for p in main_files)
lines.append("")
lines.append(f"Patient figures: {len(patient_files)}")
lines.extend(f"  patients/{p.name}" for p in patient_files)
lines.extend([
    "",
    "Patient panels include:",
    "- Venn counts with readable white background boxes;",
    "- biological labels for T1/T2/T3;",
    "- top shared aaV clonotypes by productive-read share;",
    "- top 10 TRBV changes by summed clone read share for the primary transition.",
])
index.write_text("\n".join(lines) + "\n", encoding="utf-8")

print(f"Published {len(main_files)} group figure files to {GROUP}")
print(f"Published {len(patient_files)} patient figure files to {PATIENTS}")
print(f"Index: {index}")
