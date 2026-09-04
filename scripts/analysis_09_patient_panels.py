#!/usr/bin/env python3
"""Final patient-level longitudinal reporting.

Creates one clinician/biologist-facing composite figure per measles patient with
at least two valid longitudinal clinical points after final QC.

Each figure contains:
1) aaV Venn overlap with biological time-point labels and readable boxed counts;
2) QC summary for every recorded clinical point, including excluded libraries;
3) top 10 clonotypes from the Venn intersection, shown by read share in each
   original repertoire;
4) top 10 TRBV genes by change in summed read share for the primary transition
   (T1->T2 when available; otherwise T2->T3).

Clinical meaning used in plot labels:
T1 = acute measles / admission
T2 = clinical recovery
T3 = follow-up after recovery

There is no T0 in this study. For new sequencing libraries, the suffix is the
authoritative time point: .1 -> T1, .2/.2a -> T2, .3 -> T3.
"""
from __future__ import annotations

from pathlib import Path
import shutil
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results"
OUT = BASE / "current"
PATOUT = OUT / "patients"
FIG = ROOT / "figures" / "current"
PATFIG = FIG / "patients"

META = OUT / "sample_metadata_current.csv"
CLONO = BASE / "library_trb_clonotypes.csv"

if not META.exists():
    raise RuntimeError(f"Missing final metadata: {META}. Run analysis_08_final_reporting.py first.")
if not CLONO.exists():
    raise RuntimeError(f"Missing clonotype table: {CLONO}.")

# analysis_08 may have created older patient figures. Replace them entirely so
# figures/current/patients contains only the authoritative patient panels.
if PATFIG.exists():
    shutil.rmtree(PATFIG)
PATFIG.mkdir(parents=True, exist_ok=True)
PATOUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(META)
clono = pd.read_csv(CLONO)

for c in ["analysis_include_final", "is_control"]:
    if c in meta.columns:
        meta[c] = meta[c].astype(str).str.lower().isin(["true", "1", "yes"])
for c in ["productive_trb_reads", "productive_trb_clonotypes", "source_numeric_timepoint"]:
    if c in meta.columns:
        meta[c] = pd.to_numeric(meta[c], errors="coerce")
clono["count"] = pd.to_numeric(clono["count"], errors="coerce").fillna(0.0)
for c in ["cdr3aa", "v_gene"]:
    clono[c] = clono[c].fillna("").astype(str)

PHASE_ORDER = {"T1": 1, "T2": 2, "T3": 3}
PHASE_LABEL = {
    "T1": "T1 — acute measles / admission",
    "T2": "T2 — clinical recovery",
    "T3": "T3 — follow-up after recovery",
}
PHASE_SHORT = {
    "T1": "Acute measles\n(T1)",
    "T2": "Clinical recovery\n(T2)",
    "T3": "Follow-up\n(T3)",
}

plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10})


def short_clone_label(cdr3: str, vgene: str) -> str:
    cdr3 = str(cdr3)
    if len(cdr3) > 19:
        cdr3 = cdr3[:9] + "…" + cdr3[-6:]
    return f"{cdr3} | {vgene}"


def phase_sort(phases):
    return sorted(phases, key=lambda x: PHASE_ORDER.get(x, 99))


def draw_venn(ax, point_tables: dict[str, pd.DataFrame], phases: list[str]):
    sets = {
        p: set(zip(point_tables[p]["cdr3aa"], point_tables[p]["v_gene"]))
        for p in phases
    }
    labels = [PHASE_SHORT.get(p, p) for p in phases]
    if len(phases) == 2:
        v = venn2([sets[phases[0]], sets[phases[1]]], set_labels=labels, ax=ax)
    elif len(phases) == 3:
        v = venn3([sets[phases[0]], sets[phases[1]], sets[phases[2]]], set_labels=labels, ax=ax)
    else:
        raise ValueError("Venn supports exactly 2 or 3 phases in this report.")

    # Make subset counts readable regardless of circle color/overlap.
    for lab in getattr(v, "subset_labels", []) or []:
        if lab is None:
            continue
        lab.set_fontsize(11)
        lab.set_fontweight("bold")
        lab.set_bbox({
            "facecolor": "white",
            "edgecolor": "0.75",
            "alpha": 0.88,
            "boxstyle": "round,pad=0.25",
        })
    for lab in getattr(v, "set_labels", []) or []:
        if lab is not None:
            lab.set_fontsize(10)
            lab.set_fontweight("bold")

    ax.set_title("Unique TRB aaV clonotypes\n(CDR3β amino-acid sequence + TRBV gene)")
    legend_lines = [PHASE_LABEL[p] for p in phases if p in PHASE_LABEL]
    ax.text(
        0.02, -0.12, "\n".join(legend_lines), transform=ax.transAxes,
        ha="left", va="top", fontsize=8.5,
        bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.92, boxstyle="round,pad=0.35"),
    )
    return sets


def draw_qc(ax, patient_meta: pd.DataFrame):
    ax.axis("off")
    ax.set_title("Sample QC and clinical interpretation")
    lines = []
    z = patient_meta.copy()
    z["_ord"] = z["clinical_timepoint"].map(PHASE_ORDER).fillna(99)
    z = z.sort_values(["_ord", "source_numeric_timepoint", "library_id"])
    for r in z.itertuples(index=False):
        phase = r.clinical_timepoint if pd.notna(r.clinical_timepoint) else "unresolved"
        reads = int(r.productive_trb_reads) if pd.notna(r.productive_trb_reads) else 0
        clones = int(r.productive_trb_clonotypes) if pd.notna(r.productive_trb_clonotypes) else 0
        include = bool(r.analysis_include_final)
        state = "INCLUDED" if include else "EXCLUDED"
        reason = ""
        if not include and hasattr(r, "final_exclusion_reason"):
            reason = str(r.final_exclusion_reason) if pd.notna(r.final_exclusion_reason) else ""
        lines.append(
            f"{phase:<10} {str(r.raw_sample_id):<20} "
            f"reads={reads:,}  clonotypes={clones:,}  {state}"
        )
        if reason:
            lines.append(f"            reason: {reason}")
    ax.text(
        0.01, 0.96, "\n".join(lines), transform=ax.transAxes,
        va="top", ha="left", family="monospace", fontsize=8.2,
        bbox=dict(facecolor="white", edgecolor="0.85", alpha=0.95, boxstyle="round,pad=0.4"),
    )


def shared_clone_table(point_tables: dict[str, pd.DataFrame], phases: list[str]) -> pd.DataFrame:
    sets = {
        p: set(zip(point_tables[p]["cdr3aa"], point_tables[p]["v_gene"]))
        for p in phases
    }
    shared = set.intersection(*(sets[p] for p in phases)) if phases else set()
    fraction_maps = {
        p: dict(zip(zip(point_tables[p]["cdr3aa"], point_tables[p]["v_gene"]), point_tables[p]["fraction"]))
        for p in phases
    }
    rows = []
    for clone in shared:
        cdr3, vgene = clone
        rec = {"cdr3aa": cdr3, "v_gene": vgene}
        score = 0.0
        for p in phases:
            frac = float(fraction_maps[p].get(clone, 0.0))
            rec[f"{p}_read_share_percent"] = 100.0 * frac
            score += frac
        rec["combined_fraction_score"] = score
        rows.append(rec)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values("combined_fraction_score", ascending=False)


def draw_shared_clones(ax, shared: pd.DataFrame, phases: list[str]):
    ax.set_title("Top shared clonotypes: read representation")
    if shared.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No clonotypes are shared across all displayed points.", ha="center", va="center")
        return
    top = shared.head(10).copy().iloc[::-1]
    y = np.arange(len(top))
    n = max(len(phases), 1)
    height = 0.75 / n
    offsets = (np.arange(n) - (n - 1) / 2) * height
    for off, p in zip(offsets, phases):
        vals = top[f"{p}_read_share_percent"].to_numpy()
        ax.barh(y + off, vals, height=height * 0.9, label=PHASE_SHORT.get(p, p).replace("\n", " "))
    ax.set_yticks(y, [short_clone_label(r.cdr3aa, r.v_gene) for r in top.itertuples(index=False)])
    ax.set_xlabel("Share of productive TRB reads in that repertoire (%)")
    ax.legend(fontsize=8)
    ax.grid(axis="x", alpha=0.2)
    ax.text(
        0.01, -0.15,
        "Clones are selected from the Venn intersection and ranked by summed read share across displayed points.",
        transform=ax.transAxes, fontsize=8, ha="left", va="top",
    )


def choose_primary_transition(phases: list[str]) -> tuple[str, str]:
    phase_set = set(phases)
    for a, b in [("T1", "T2"), ("T2", "T3")]:
        if a in phase_set and b in phase_set:
            return a, b
    ps = phase_sort(phases)
    if len(ps) < 2:
        raise ValueError("Need at least two valid phases for a transition.")
    return ps[0], ps[1]


def trbv_change_table(point_tables: dict[str, pd.DataFrame], a: str, b: str) -> pd.DataFrame:
    def one(t):
        x = t.groupby("v_gene", as_index=False)["count"].sum()
        total = float(x["count"].sum())
        x["share"] = x["count"] / total if total > 0 else 0.0
        return x.set_index("v_gene")["share"]
    x = one(point_tables[a])
    y = one(point_tables[b])
    genes = x.index.union(y.index)
    out = pd.DataFrame({
        "v_gene": genes,
        f"{a}_read_share_percent": 100.0 * x.reindex(genes, fill_value=0).to_numpy(),
        f"{b}_read_share_percent": 100.0 * y.reindex(genes, fill_value=0).to_numpy(),
    })
    out["delta_percentage_points"] = out[f"{b}_read_share_percent"] - out[f"{a}_read_share_percent"]
    out["abs_delta"] = out["delta_percentage_points"].abs()
    return out.sort_values("abs_delta", ascending=False)


def draw_trbv_change(ax, table: pd.DataFrame, a: str, b: str):
    top = table.head(10).sort_values("delta_percentage_points")
    y = np.arange(len(top))
    vals = top["delta_percentage_points"].to_numpy()
    colors = ["C1" if v < 0 else "C0" for v in vals]
    ax.barh(y, vals, color=colors)
    ax.axvline(0, color="black", linewidth=0.9)
    ax.set_yticks(y, top["v_gene"])
    ax.set_xlabel(f"Change in TRBV share, {b} − {a} (percentage points)")
    ax.set_title(
        "Top 10 TRBV shifts by summed clone read share\n"
        f"{PHASE_SHORT.get(a, a).replace(chr(10), ' ')} → {PHASE_SHORT.get(b, b).replace(chr(10), ' ')}"
    )
    span = max(np.max(np.abs(vals)) if len(vals) else 0, 0.1)
    for yi, v in enumerate(vals):
        ax.text(v + (0.025 * span if v >= 0 else -0.025 * span), yi, f"{v:+.2f}",
                va="center", ha="left" if v >= 0 else "right", fontsize=8)
    ax.grid(axis="x", alpha=0.2)


valid_meta = meta[meta["analysis_include_final"]].copy()
case_meta = meta[meta["disease_status"].eq("measles")].copy()

# Build all valid patient/phase repertoires once. This avoids repeatedly scanning
# the large clonotype table for every patient panel.
valid_map = valid_meta[["library_id", "patient_id", "clinical_timepoint"]].drop_duplicates()
cl_valid = clono.merge(valid_map, on="library_id", how="inner", suffixes=("", "_meta"))
if "patient_id_meta" in cl_valid.columns:
    cl_valid["patient_id"] = cl_valid["patient_id_meta"]
if "clinical_timepoint_meta" in cl_valid.columns:
    cl_valid["clinical_timepoint"] = cl_valid["clinical_timepoint_meta"]
cl_valid = cl_valid[(cl_valid["cdr3aa"] != "") & (cl_valid["v_gene"] != "")].copy()
phase_clones = cl_valid.groupby(
    ["patient_id", "clinical_timepoint", "cdr3aa", "v_gene"], as_index=False
)["count"].sum()
point_cache = {}
for (pid, phase), z in phase_clones.groupby(["patient_id", "clinical_timepoint"], sort=False):
    z = z[["cdr3aa", "v_gene", "count"]].copy()
    total = float(z["count"].sum())
    z["fraction"] = z["count"] / total if total > 0 else 0.0
    point_cache[(pid, phase)] = z


def aggregate_point_cached(patient_id: str, phase: str) -> pd.DataFrame:
    return point_cache.get((patient_id, phase), pd.DataFrame(columns=["cdr3aa", "v_gene", "count", "fraction"]))


# Patient list is based on at least two distinct valid clinical phases.
patient_ids = []
for pid, g in valid_meta[valid_meta["disease_status"].eq("measles")].groupby("patient_id"):
    phases = [p for p in g["clinical_timepoint"].dropna().unique() if p in PHASE_ORDER]
    if len(set(phases)) >= 2:
        patient_ids.append(pid)

summary_rows = []
for pid in sorted(patient_ids):
    pm_valid = valid_meta[valid_meta["patient_id"].eq(pid)].copy()
    phases_all = phase_sort([p for p in pm_valid["clinical_timepoint"].dropna().unique() if p in PHASE_ORDER])

    # Venn uses all valid study phases; the study has at most T1/T2/T3.
    phases = phases_all
    if len(phases) < 2:
        continue

    point_tables = {p: aggregate_point_cached(pid, p) for p in phases}
    if any(point_tables[p].empty for p in phases):
        continue

    shared = shared_clone_table(point_tables, phases)
    a, b = choose_primary_transition(phases_all)
    trbv = trbv_change_table(point_tables, a, b)

    # Persist underlying patient-level tables for interpretation/audit.
    shared_csv = PATOUT / f"{pid}_shared_clonotypes.csv"
    trbv_csv = PATOUT / f"{pid}_{a}_to_{b}_TRBV_changes.csv"
    shared.to_csv(shared_csv, index=False)
    trbv.drop(columns=["abs_delta"]).to_csv(trbv_csv, index=False)

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], wspace=0.32, hspace=0.42)
    ax_venn = fig.add_subplot(gs[0, 0])
    ax_qc = fig.add_subplot(gs[0, 1])
    ax_shared = fig.add_subplot(gs[1, 0])
    ax_trbv = fig.add_subplot(gs[1, 1])

    draw_venn(ax_venn, point_tables, phases)
    draw_qc(ax_qc, case_meta[case_meta["patient_id"].eq(pid)])
    draw_shared_clones(ax_shared, shared, phases)
    draw_trbv_change(ax_trbv, trbv, a, b)

    age = pm_valid["age_group"].iloc[0] if "age_group" in pm_valid.columns else ""
    fig.suptitle(
        f"{pid}: longitudinal TRB repertoire summary ({age})",
        fontsize=17, fontweight="bold", y=0.995,
    )
    fig.text(
        0.5, 0.005,
        "Venn counts are numbers of unique aaV clonotypes, not read counts. "
        "Shared-clone and TRBV panels quantify read representation.",
        ha="center", va="bottom", fontsize=9,
    )

    out_png = PATFIG / f"{pid}_repertoire_QC_and_overlap.png"
    out_pdf = PATFIG / f"{pid}_repertoire_QC_and_overlap.pdf"
    fig.savefig(out_png, dpi=210, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    summary_rows.append({
        "patient_id": pid,
        "age_group": age,
        "venn_phases": ";".join(phases),
        "primary_transition": f"{a}->{b}",
        "n_shared_clonotypes": int(len(shared)),
        "figure": str(out_png.relative_to(ROOT)),
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT / "patient_panel_index.csv", index=False)

# Hard CA1 check: T1 must be visible in QC but excluded, while T2/T3 are valid.
ca = meta[meta["patient_id"].eq("СА1_013")].copy()
if not ca.empty:
    required = {"T1", "T2", "T3"}
    observed = set(ca["clinical_timepoint"].dropna())
    if not required.issubset(observed):
        raise RuntimeError(f"СА1_013 metadata lost a clinical point. Observed: {sorted(observed)}")
    ca_t1 = ca[ca["clinical_timepoint"].eq("T1")]
    ca_t23 = ca[ca["clinical_timepoint"].isin(["T2", "T3"])]
    if ca_t1.empty or ca_t1["analysis_include_final"].any():
        raise RuntimeError("СА1_013 T1 must exist but remain excluded after QC.")
    if set(ca_t23.loc[ca_t23["analysis_include_final"], "clinical_timepoint"]) != {"T2", "T3"}:
        raise RuntimeError("СА1_013 must have valid T2 and T3 for final overlap reporting.")
    ca_fig = PATFIG / "СА1_013_repertoire_QC_and_overlap.png"
    if not ca_fig.exists() or ca_fig.stat().st_size == 0:
        raise RuntimeError("СА1_013 final patient panel was not generated.")

print(f"Final patient panels generated: {len(summary_rows)}")
print(f"Figures: {PATFIG}")
print(f"Patient tables: {PATOUT}")
