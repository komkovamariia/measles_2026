#!/usr/bin/env python3
"""Patient level panels for expanded TRB clonotypes

The primary patient comparison is T1 to T2 when both points pass quality control
If T1 is unavailable but T2 and T3 are valid the panel shows T2 to T3 as follow up
Overlap is never calculated on the complete low frequency repertoire
"""
from pathlib import Path
import os
import shutil
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FIG = ROOT / "figures" / "current"
PATFIG = FIG / "patients"
PATOUT = CURRENT / "patients"
if PATFIG.exists():
    shutil.rmtree(PATFIG)
PATFIG.mkdir(parents=True, exist_ok=True)
PATOUT.mkdir(parents=True, exist_ok=True)

EXPANDED_FRACTION = 0.001
meta = pd.read_csv(CURRENT / "sample_metadata_current.csv")
metrics = pd.read_csv(RESULTS / "diversity_metrics.csv")
clones = pd.read_csv(RESULTS / "combined_trb_clonotypes.csv")
for c in ["analysis_include_final", "is_control"]:
    if c in meta.columns:
        meta[c] = meta[c].astype(str).str.lower().isin(["true", "1", "yes"])
clones["fraction"] = pd.to_numeric(clones["fraction"], errors="coerce").fillna(0.0)
clones["expanded"] = clones["fraction"] >= EXPANDED_FRACTION

specificity = {}
match_path = RESULTS / "antigen_specific_clonotypes.csv"
if match_path.exists() and match_path.stat().st_size > 0:
    m = pd.read_csv(match_path)
    if {"sample_id", "cdr3aa", "pathogen"}.issubset(m.columns):
        m = m.dropna(subset=["sample_id", "cdr3aa", "pathogen"])
        for key, z in m.groupby(["sample_id", "cdr3aa"]):
            specificity[key] = "; ".join(sorted(set(z.pathogen.astype(str))))

POINT_LABEL = {
    "T1": "T1 early infection",
    "T2": "T2 recovery",
    "T3": "T3 follow up",
}
POINT_ORDER = {"T1": 1, "T2": 2, "T3": 3}


def point_table(patient, point):
    z = clones[(clones.patient_id == patient) & clones.clinical_timepoint.eq(point)].copy()
    z["clone_key"] = list(zip(z.cdr3aa.astype(str), z.v_gene.astype(str)))
    return z.set_index("clone_key")


def choose_pair(points):
    s = set(points)
    if {"T1", "T2"}.issubset(s):
        return "T1", "T2"
    if {"T2", "T3"}.issubset(s):
        return "T2", "T3"
    return None


def draw_qc(ax, patient):
    ax.axis("off")
    ax.set_title("Library quality control and aggregation")
    z = meta[meta.patient_id.eq(patient)].copy()
    z["point_order"] = z.clinical_timepoint.map(POINT_ORDER).fillna(99)
    z = z.sort_values(["point_order", "library_id"])
    lines = []
    for r in z.itertuples(index=False):
        state = "included" if r.analysis_include_final else "excluded"
        lines.append(
            f"{str(r.clinical_timepoint):<4} {str(r.library_id):<16} clones {int(r.productive_trb_clonotypes):>6}  {state}"
        )
    t2 = z[z.clinical_timepoint.eq("T2") & z.analysis_include_final]
    if len(t2) > 1:
        lines.append("")
        lines.append(f"T2 has {len(t2)} accepted libraries")
        lines.append("Normalized clonotype frequencies are averaged")
        lines.append("The patient contributes one T2 observation")
    ax.text(0.01, 0.96, "\n".join(lines), transform=ax.transAxes, va="top", ha="left", family="monospace", fontsize=8.5)


def draw_metrics(ax, patient, a, b):
    z = metrics[(metrics.patient_id == patient) & metrics.clinical_timepoint.isin([a, b])].set_index("clinical_timepoint")
    ax2 = ax.twinx()
    x = [1, 2]
    clonality = [float(z.loc[a, "clonality"]), float(z.loc[b, "clonality"])]
    shannon = [float(z.loc[a, "shannon"]), float(z.loc[b, "shannon"])]
    ax.plot(x, clonality, marker="o", linewidth=2, label="Clonality")
    ax2.plot(x, shannon, marker="s", linewidth=2, label="Shannon diversity")
    ax.set_xticks(x, [POINT_LABEL[a], POINT_LABEL[b]])
    ax.set_ylabel("Clonality")
    ax2.set_ylabel("Shannon diversity")
    ax.set_title("Complete repertoire diversity")
    ax.grid(axis="y", alpha=0.2)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], fontsize=8, loc="best")


def expanded_sets(patient, a, b):
    x = point_table(patient, a)
    y = point_table(patient, b)
    sx = set(x.index[x.expanded])
    sy = set(y.index[y.expanded])
    return x, y, sx, sy


def draw_venn(ax, patient, a, b):
    x, y, sx, sy = expanded_sets(patient, a, b)
    v = venn2([sx, sy], set_labels=[POINT_LABEL[a], POINT_LABEL[b]], ax=ax)
    for lab in getattr(v, "subset_labels", []) or []:
        if lab is not None:
            lab.set_fontweight("bold")
            lab.set_bbox({"facecolor": "white", "edgecolor": "0.8", "alpha": 0.9, "boxstyle": "round,pad=0.2"})
    ax.set_title("Expanded clonotype overlap")
    return x, y, sx, sy


def fate_table(patient, a, b):
    x, y, sx, sy = expanded_sets(patient, a, b)
    rows = []
    for key in sorted(sx):
        fa = float(x.at[key, "fraction"])
        fb = float(y.at[key, "fraction"]) if key in y.index else 0.0
        if key in sy:
            fate = "retained expanded"
        elif key in y.index:
            fate = "detected below threshold"
        else:
            fate = "not detected later"
        rows.append({"cdr3aa": key[0], "v_gene": key[1], f"{a}_fraction": fa, f"{b}_fraction": fb, "fate": fate})
    return pd.DataFrame(rows)


def draw_fate(ax, fate, a):
    cats = ["retained expanded", "detected below threshold", "not detected later"]
    counts = [int((fate.fate == c).sum()) if not fate.empty else 0 for c in cats]
    labels = ["Retained expanded", "Detected below threshold", "Not detected later"]
    ax.bar(np.arange(3), counts)
    ax.set_xticks(np.arange(3), labels, rotation=20, ha="right")
    ax.set_ylabel(f"Number of expanded {a} clonotypes")
    ax.set_title(f"Fate of expanded {a} clonotypes")
    ax.grid(axis="y", alpha=0.2)


def new_expansion_table(patient, a, b):
    x, y, sx, sy = expanded_sets(patient, a, b)
    rows = []
    for key in sorted(sy - sx):
        fa = float(x.at[key, "fraction"]) if key in x.index else 0.0
        fb = float(y.at[key, "fraction"])
        origin = "present below threshold earlier" if key in x.index else "not detected earlier"
        match = specificity.get((f"{patient}__{b}", key[0]), "")
        rows.append({
            "cdr3aa": key[0], "v_gene": key[1], f"{a}_fraction": fa, f"{b}_fraction": fb,
            "origin": origin, "exact_reference_match": match,
        })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(f"{b}_fraction", ascending=False)


def short_label(row):
    cdr3 = str(row.cdr3aa)
    if len(cdr3) > 18:
        cdr3 = cdr3[:9] + "…" + cdr3[-5:]
    base = f"{cdr3}  {row.v_gene}"
    match = str(row.exact_reference_match) if pd.notna(row.exact_reference_match) else ""
    return base if not match else f"{base}  {match}"


def draw_new_expansions(ax, new, a, b):
    if new.empty:
        ax.axis("off")
        ax.set_title(f"New expanded clonotypes at {b}")
        ax.text(0.5, 0.5, "No new expanded clonotypes", ha="center", va="center")
        return
    top = new.head(10).iloc[::-1]
    y = np.arange(len(top))
    ax.barh(y, 100 * top[f"{b}_fraction"].to_numpy(float))
    ax.set_yticks(y, [short_label(r) for r in top.itertuples(index=False)])
    ax.set_xlabel(f"Share of productive TRB reads at {b} percent")
    ax.set_title(f"Largest new expanded clonotypes at {b}")
    ax.grid(axis="x", alpha=0.2)

valid = metrics[metrics.disease_status.eq("measles")].copy()
patient_ids = sorted(valid.patient_id.unique())
summary = []

for patient in patient_ids:
    points = sorted(valid.loc[valid.patient_id.eq(patient), "clinical_timepoint"].unique(), key=lambda p: POINT_ORDER.get(p, 99))
    pair = choose_pair(points)
    if pair is None:
        continue
    a, b = pair
    fate = fate_table(patient, a, b)
    new = new_expansion_table(patient, a, b)
    fate.to_csv(PATOUT / f"{patient}_{a}_to_{b}_expanded_fates.csv", index=False)
    new.to_csv(PATOUT / f"{patient}_{a}_to_{b}_new_expansions.csv", index=False)

    fig = plt.figure(figsize=(16, 11.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1.0, 1.25], hspace=0.38, wspace=0.38)
    ax_metrics = fig.add_subplot(gs[0, 0])
    ax_venn = fig.add_subplot(gs[0, 1])
    ax_qc = fig.add_subplot(gs[0, 2])
    ax_fate = fig.add_subplot(gs[1, 0:2])
    ax_new = fig.add_subplot(gs[1, 2])

    draw_metrics(ax_metrics, patient, a, b)
    x, y, sx, sy = draw_venn(ax_venn, patient, a, b)
    draw_qc(ax_qc, patient)
    draw_fate(ax_fate, fate, a)
    draw_new_expansions(ax_new, new, a, b)

    age = valid.loc[valid.patient_id.eq(patient), "age_group"].iloc[0]
    age_text = "child younger than 18 years" if age == "child" else "adult 18 years and older"
    fig.suptitle(f"Patient {patient}: expanded TRB clonotypes from {a} to {b}, {age_text}", fontsize=16, fontweight="bold")
    fig.text(
        0.5, 0.008,
        "Expanded clonotype means aaV frequency at least 0.10 percent of productive TRB reads. Low frequency background is excluded from overlap",
        ha="center", fontsize=8.5,
    )
    fig.savefig(PATFIG / f"{patient}_expanded_repertoire.png", dpi=210, bbox_inches="tight")
    fig.savefig(PATFIG / f"{patient}_expanded_repertoire.pdf", bbox_inches="tight")
    plt.close(fig)

    union = sx | sy
    summary.append({
        "patient_id": patient,
        "age_group": age,
        "comparison": f"{a} to {b}",
        "n_expanded_earlier": len(sx),
        "n_expanded_later": len(sy),
        "n_shared_expanded": len(sx & sy),
        "expanded_jaccard": len(sx & sy) / len(union) if union else np.nan,
        "n_new_expanded_later": len(sy - sx),
        "figure": f"figures/current/patients/{patient}_expanded_repertoire.png",
    })

index = pd.DataFrame(summary)
index.to_csv(CURRENT / "patient_panel_index.csv", index=False)
print(f"Patient panels generated: {len(index)}")
print("Primary overlap uses expanded clonotypes only")
