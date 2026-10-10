#!/usr/bin/env python3
"""Publish a compact final figure set with a shallow directory structure."""

from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
CURRENT = FIG / "current"
CURRENT_PATIENTS = CURRENT / "patients"

OUT = {
    "main": FIG / "main",
    "patients": FIG / "patients",
    "downsampling_pairwise": FIG / "downsampling_pairwise",
    "sensitivity_thresholds": FIG / "sensitivity_thresholds",
    "validation": FIG / "validation",
    "supplementary": FIG / "supplementary",
}

marker = CURRENT / "GENERATED_CURRENT.txt"
if not marker.exists():
    raise RuntimeError("Final figure marker is missing")

# Remove stale published trees, including the former catch-all group directory.
for name in ["group", *OUT]:
    d = FIG / name
    if d.exists():
        shutil.rmtree(d)
for d in OUT.values():
    d.mkdir(parents=True, exist_ok=True)


def copy_stem(stem, target, required=False):
    copied = []
    for ext in [".png", ".pdf"]:
        src = CURRENT / f"{stem}{ext}"
        if src.exists() and src.stat().st_size > 0:
            dst = target / src.name
            shutil.copy2(src, dst)
            copied.append(dst)
    if required and not copied:
        raise RuntimeError(f"Required final figure is missing: {stem}")
    return copied


published = {k: [] for k in OUT}

# Main biological conclusions from the clone-level specification.
published["main"] += copy_stem(
    "21_primary_fate_pairwise_rarefied", OUT["main"], required=True
)
published["main"] += copy_stem(
    "24_lost_vs_retained_enrichment", OUT["main"], required=False
)
published["main"] += copy_stem(
    "25_sequence_neighbourhoods", OUT["main"], required=False
)

# Sampling and threshold sensitivity each have one dedicated conclusion figure.
published["downsampling_pairwise"] += copy_stem(
    "22_sampling_robust_clone_loss", OUT["downsampling_pairwise"], required=True
)
published["sensitivity_thresholds"] += copy_stem(
    "23_threshold_sensitivity_pairwise_rarefied",
    OUT["sensitivity_thresholds"],
    required=True,
)
published["validation"] += copy_stem(
    "26_validation_identity_T3_batch", OUT["validation"], required=True
)

# Patient panels remain an explicit, separate biological layer.
patient_files = (
    sorted(
        p for p in CURRENT_PATIENTS.iterdir() if p.suffix.lower() in {".png", ".pdf"}
    )
    if CURRENT_PATIENTS.exists()
    else []
)
if not patient_files:
    raise RuntimeError("No final patient figures were generated")
for src in patient_files:
    dst = OUT["patients"] / src.name
    shutil.copy2(src, dst)
    published["patients"].append(dst)

# Preserve only useful old context. Old TRBV heatmaps and duplicate exploratory plots are not published.
for stem in [
    "09_specificity_reference_coverage",
    "10_complete_clonotype_overlap_summary",
    "11_expanded_clonotype_overlap_summary",
    "12_expanded_T1_fate_summary",
    "13_patient_overlap_similarity_summary",
]:
    published["supplementary"] += copy_stem(stem, OUT["supplementary"], required=False)

index = FIG / "FIGURE_INDEX.txt"
lines = [
    "FINAL FIGURE SET",
    "",
    marker.read_text(encoding="utf-8").strip(),
    "",
    "MAIN",
    "21 primary fate: retained expanded, detected below threshold and not detected after pair-specific repeated rarefaction",
    "24 enrichment is published only when at least one patient-aware feature passes FDR 0.05",
    "25 sequence neighbourhoods are published only when a cross-patient affected cluster exists",
    "",
    "DOWNSAMPLING PAIRWISE",
    "22 sampling robustness: observed contraction and empirical P_detection_if_unchanged from the internal T1 null",
    "",
    "SENSITIVITY THRESHOLDS",
    "23 one combined sensitivity figure for 0.05, 0.10, 0.20, 0.50, 0.75 and 1.00 percent",
    "",
    "VALIDATION",
    "26 combined identity QC, available T3 follow-up and GA1_003 cross-batch sensitivity",
    "",
    "PATIENTS",
    "One panel per accepted longitudinal patient, retaining complete and expanded repertoire views",
    "",
    "SUPPLEMENTARY",
    "Only selected contextual figures from the established analysis are retained. Old TRBV heatmaps are intentionally excluded from the final figure tree.",
    "",
    "Published files",
]
for section in [
    "main",
    "downsampling_pairwise",
    "sensitivity_thresholds",
    "validation",
    "patients",
    "supplementary",
]:
    lines.append(f"[{section}]")
    lines.extend(str(p.relative_to(FIG)) for p in published[section])

lines.extend(
    [
        "",
        "Analysis logic",
        "Physical tracking uses CDR3 nucleotide plus TRBV plus TRBJ",
        "Functional and public analysis uses CDR3 amino acid plus TRBV plus TRBJ",
        "Primary T1 expanded threshold is 0.10 percent",
        "Sensitivity thresholds are 0.05 0.10 0.20 0.50 0.75 and 1.00 percent",
        "Each T1-T2 pair is repeatedly rarefied to its own within-pair minimum productive-read depth",
        "At least 500 rarefaction iterations are required",
        "P_detection_if_unchanged is empirical from the repeated internal T1 sampling null",
        "T2 zero values use a detection upper bound and no arbitrary pseudocount",
        "Patient-aware enrichment adjusts for baseline T1 clone frequency and within-patient dependence",
        "Sequence-neighbour clusters are edit-distance neighbourhoods and are not labelled TCRdist clusters",
        "HLA compatibility remains not available without patient HLA typing",
        "БА1_071 T2 replicate libraries are rarefied separately and averaged only after within-library normalization",
    ]
)
index.write_text("\n".join(lines) + "\n", encoding="utf-8")

# Remove the internal staging tree after successful publication so figures/ stays clean.
shutil.rmtree(CURRENT)

for section in published:
    print(section, len(published[section]), "files")
print("Published compact final figure tree")
