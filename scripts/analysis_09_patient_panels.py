#!/usr/bin/env python3
"""Create patient level TRB repertoire panels with complete and expanded clone overlap.

The primary comparison is T1 to T2 when both points pass final quality control.
If T1 is unavailable but T2 and T3 are valid, the panel shows T2 to T3.
Complete repertoire overlap is shown for biological context.
The immune amnesia overlap analysis is restricted to expanded aaV clonotypes.
"""
from pathlib import Path
import os
import shutil
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from plot_style import COLORS, FATE_COLORS, TIME_COLORS, finish_axis

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

for col in ["analysis_include_final", "is_control"]:
    if col in meta.columns:
        meta[col] = meta[col].astype(str).str.lower().isin(["true", "1", "yes"])

clones["fraction"] = pd.to_numeric(clones["fraction"], errors="coerce").fillna(0.0)
clones["count"] = pd.to_numeric(clones["count"], errors="coerce").fillna(0.0)
clones["expanded"] = clones["fraction"] >= EXPANDED_FRACTION

POINT_LABEL = {
    "T1": "T1 early infection",
    "T2": "T2 recovery",
    "T3": "T3 follow up",
}
POINT_ORDER = {"T1": 1, "T2": 2, "T3": 3}


def choose_pair(points):
    point_set = set(points)
    if {"T1", "T2"}.issubset(point_set):
        return "T1", "T2"
    if {"T2", "T3"}.issubset(point_set):
        return "T2", "T3"
    return None


def point_table(patient, point):
    z = clones[
        clones["patient_id"].eq(patient)
        & clones["clinical_timepoint"].eq(point)
    ].copy()
    z["clone_key"] = list(zip(z["cdr3aa"].astype(str), z["v_gene"].astype(str)))
    return z.set_index("clone_key", drop=False)


def point_metric(patient, point, column):
    z = metrics[
        metrics["patient_id"].eq(patient)
        & metrics["clinical_timepoint"].eq(point)
    ]
    if z.empty:
        return np.nan
    return float(z.iloc[0][column])


def point_summary(patient, point):
    z = point_table(patient, point)
    expanded = z[z["expanded"]]
    return {
        "all_clones": int(len(z)),
        "all_read_equivalents": float(z["count"].sum()),
        "expanded_clones": int(len(expanded)),
        "expanded_read_equivalents": float(expanded["count"].sum()),
        "expanded_read_fraction": float(expanded["fraction"].sum()),
        "clonality": point_metric(patient, point, "clonality"),
        "shannon": point_metric(patient, point, "shannon"),
        "top10_fraction": point_metric(patient, point, "top10_fraction"),
    }


def expanded_sets(patient, earlier, later):
    x = point_table(patient, earlier)
    y = point_table(patient, later)
    sx = set(x.index[x["expanded"]])
    sy = set(y.index[y["expanded"]])
    return x, y, sx, sy


def format_read_equivalents(value):
    return f"{int(round(value)):,}"


def style_venn(v):
    for patch, color in zip(getattr(v, "patches", []) or [], [COLORS["blue"], COLORS["orange"], COLORS["teal"]]):
        if patch is not None:
            patch.set_facecolor(color)
            patch.set_edgecolor("white")
            patch.set_alpha(0.72)
    for label in getattr(v, "subset_labels", []) or []:
        if label is not None:
            label.set_fontsize(8.2)
            label.set_fontweight("bold")
            label.set_bbox({
                "facecolor": "white",
                "edgecolor": "0.75",
                "alpha": 0.94,
                "boxstyle": "round,pad=0.22",
            })
    for label in getattr(v, "set_labels", []) or []:
        if label is not None:
            label.set_fontsize(9.5)
            label.set_fontweight("bold")


def draw_summary(ax, patient, earlier, later):
    ax.axis("off")
    ax.set_title("Biological interpretation guide and library quality", fontsize=12, pad=10)

    a = point_summary(patient, earlier)
    b = point_summary(patient, later)

    age = metrics.loc[metrics["patient_id"].eq(patient), "age_group"].iloc[0]
    age_text = "child younger than 18 years" if age == "child" else "adult 18 years and older"

    z = meta[meta["patient_id"].eq(patient)].copy()
    z["point_order"] = z["clinical_timepoint"].map(POINT_ORDER).fillna(99)
    z = z.sort_values(["point_order", "library_id"])

    lines = [
        f"Patient: {patient}",
        f"Age group: {age_text}",
        "",
        "Definitions",
        "Clonality = 1 - Shannon entropy / ln(number of clonotypes)",
        "Higher clonality means stronger dominance by a smaller number of clones",
        "Shannon diversity increases when the repertoire is richer and more even",
        "Repertoire fraction = fraction of productive TRB reads assigned to a clonotype",
        "Expanded clonotype = aaV clone with repertoire fraction at least 0.10%",
        "aaV clone = CDR3 beta amino acid sequence together with TRBV gene",
        "Read equivalents = repertoire fraction multiplied by mean productive reads",
        "For replicated points read equivalents are normalized estimates, not raw summed reads",
        "",
        f"{POINT_LABEL[earlier]}",
        f"  all aaV clonotypes: {a['all_clones']:,}",
        f"  productive read equivalents: {format_read_equivalents(a['all_read_equivalents'])}",
        f"  expanded clonotypes: {a['expanded_clones']}",
        f"  expanded read equivalents: {format_read_equivalents(a['expanded_read_equivalents'])}",
        f"  expanded repertoire fraction: {100 * a['expanded_read_fraction']:.1f}%",
        f"  clonality: {a['clonality']:.3f}",
        f"  Shannon diversity: {a['shannon']:.3f}",
        "",
        f"{POINT_LABEL[later]}",
        f"  all aaV clonotypes: {b['all_clones']:,}",
        f"  productive read equivalents: {format_read_equivalents(b['all_read_equivalents'])}",
        f"  expanded clonotypes: {b['expanded_clones']}",
        f"  expanded read equivalents: {format_read_equivalents(b['expanded_read_equivalents'])}",
        f"  expanded repertoire fraction: {100 * b['expanded_read_fraction']:.1f}%",
        f"  clonality: {b['clonality']:.3f}",
        f"  Shannon diversity: {b['shannon']:.3f}",
        "",
        "Accepted libraries",
    ]

    for row in z.itertuples(index=False):
        state = "included" if bool(row.analysis_include_final) else "excluded"
        lines.append(
            f"  {row.clinical_timepoint}: {row.library_id}, "
            f"{int(row.productive_trb_clonotypes):,} productive clonotypes, {state}"
        )

    later_reps = z[
        z["clinical_timepoint"].eq(later)
        & z["analysis_include_final"]
    ]
    if len(later_reps) > 1:
        lines.extend([
            "",
            f"{later} has {len(later_reps)} accepted libraries",
            "Clonotype frequencies are normalized within each library and then averaged",
            f"The patient contributes one biological {later} observation",
        ])

    ax.text(
        0.01,
        0.98,
        "\n".join(lines),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.15,
        linespacing=1.22,
        bbox={
            "facecolor": "white",
            "edgecolor": "0.8",
            "alpha": 0.98,
            "boxstyle": "round,pad=0.45",
        },
    )


def draw_all_clone_venn(ax, patient, earlier, later):
    x = point_table(patient, earlier)
    y = point_table(patient, later)

    sx = set(x.index)
    sy = set(y.index)
    left = sx - sy
    right = sy - sx
    shared = sx & sy

    exp_x = set(x.index[x["expanded"]])
    exp_y = set(y.index[y["expanded"]])

    v = venn2(
        [sx, sy],
        set_labels=[POINT_LABEL[earlier], POINT_LABEL[later]],
        ax=ax,
    )
    style_venn(v)

    left_fraction = float(x.loc[list(left), "fraction"].sum()) if left else 0.0
    left_reads = float(x.loc[list(left), "count"].sum()) if left else 0.0
    right_fraction = float(y.loc[list(right), "fraction"].sum()) if right else 0.0
    right_reads = float(y.loc[list(right), "count"].sum()) if right else 0.0

    shared_x_fraction = float(x.loc[list(shared), "fraction"].sum()) if shared else 0.0
    shared_y_fraction = float(y.loc[list(shared), "fraction"].sum()) if shared else 0.0
    shared_x_reads = float(x.loc[list(shared), "count"].sum()) if shared else 0.0
    shared_y_reads = float(y.loc[list(shared), "count"].sum()) if shared else 0.0

    left_label = v.get_label_by_id("10")
    right_label = v.get_label_by_id("01")
    shared_label = v.get_label_by_id("11")

    if left_label is not None:
        left_label.set_text(
            f"{len(left):,} clones\n"
            f"{100 * left_fraction:.1f}% reads\n"
            f"~{format_read_equivalents(left_reads)} read eq\n"
            f"expanded: {len(left & exp_x)}"
        )
    if right_label is not None:
        right_label.set_text(
            f"{len(right):,} clones\n"
            f"{100 * right_fraction:.1f}% reads\n"
            f"~{format_read_equivalents(right_reads)} read eq\n"
            f"expanded: {len(right & exp_y)}"
        )
    if shared_label is not None:
        shared_label.set_text(
            f"{len(shared):,} shared\n"
            f"{earlier}: {100 * shared_x_fraction:.1f}%  ~{format_read_equivalents(shared_x_reads)}\n"
            f"{later}: {100 * shared_y_fraction:.1f}%  ~{format_read_equivalents(shared_y_reads)}\n"
            f"expanded in both: {len(shared & exp_x & exp_y)}"
        )

    ax.set_title(
        "Complete productive aaV repertoire overlap\n"
        "Numbers show unique clonotypes, repertoire fraction and read equivalents",
        fontsize=11.5,
        pad=10,
    )


def draw_expanded_venn(ax, patient, earlier, later):
    x, y, sx, sy = expanded_sets(patient, earlier, later)

    left = sx - sy
    right = sy - sx
    shared = sx & sy

    v = venn2(
        [sx, sy],
        set_labels=[POINT_LABEL[earlier], POINT_LABEL[later]],
        ax=ax,
    )
    style_venn(v)

    left_fraction = float(x.loc[list(left), "fraction"].sum()) if left else 0.0
    left_reads = float(x.loc[list(left), "count"].sum()) if left else 0.0
    right_fraction = float(y.loc[list(right), "fraction"].sum()) if right else 0.0
    right_reads = float(y.loc[list(right), "count"].sum()) if right else 0.0

    shared_x_fraction = float(x.loc[list(shared), "fraction"].sum()) if shared else 0.0
    shared_y_fraction = float(y.loc[list(shared), "fraction"].sum()) if shared else 0.0
    shared_x_reads = float(x.loc[list(shared), "count"].sum()) if shared else 0.0
    shared_y_reads = float(y.loc[list(shared), "count"].sum()) if shared else 0.0

    left_label = v.get_label_by_id("10")
    right_label = v.get_label_by_id("01")
    shared_label = v.get_label_by_id("11")

    if left_label is not None:
        left_label.set_text(
            f"{len(left)} clones\n"
            f"{100 * left_fraction:.1f}% reads\n"
            f"~{format_read_equivalents(left_reads)} read eq"
        )
    if right_label is not None:
        right_label.set_text(
            f"{len(right)} clones\n"
            f"{100 * right_fraction:.1f}% reads\n"
            f"~{format_read_equivalents(right_reads)} read eq"
        )
    if shared_label is not None:
        shared_label.set_text(
            f"{len(shared)} shared\n"
            f"{earlier}: {100 * shared_x_fraction:.1f}%  ~{format_read_equivalents(shared_x_reads)}\n"
            f"{later}: {100 * shared_y_fraction:.1f}%  ~{format_read_equivalents(shared_y_reads)}"
        )

    ax.set_title(
        "Expanded clonotype overlap\n"
        "Expanded = aaV clonotype representing at least 0.10% of productive TRB reads",
        fontsize=11.5,
        pad=10,
    )
    return x, y, sx, sy


def fate_table(patient, earlier, later):
    x, y, sx, sy = expanded_sets(patient, earlier, later)
    rows = []

    for key in sorted(sx):
        earlier_fraction = float(x.at[key, "fraction"])
        later_fraction = float(y.at[key, "fraction"]) if key in y.index else 0.0
        earlier_reads = float(x.at[key, "count"])
        later_reads = float(y.at[key, "count"]) if key in y.index else 0.0

        if key in sy:
            fate = "Retained expanded"
        elif key in y.index:
            fate = "Detected below threshold"
        else:
            fate = "Not detected later"

        rows.append({
            "cdr3aa": key[0],
            "v_gene": key[1],
            f"{earlier}_fraction": earlier_fraction,
            f"{later}_fraction": later_fraction,
            f"{earlier}_read_equivalents": earlier_reads,
            f"{later}_read_equivalents": later_reads,
            "fate": fate,
        })

    return pd.DataFrame(rows)


def draw_fate(ax, fate, earlier):
    categories = [
        "Retained expanded",
        "Detected below threshold",
        "Not detected later",
    ]
    counts = [
        int((fate["fate"] == category).sum()) if not fate.empty else 0
        for category in categories
    ]

    positions = np.arange(len(categories))
    ax.bar(positions, counts, color=[FATE_COLORS[x] for x in categories], width=0.68)
    ax.set_xticks(
        positions,
        ["Retained expanded\nlater ≥ 0.10%",
         "Detected below threshold\n0 < later < 0.10%",
         "Not detected later\nzero observed reads"],
    )
    ax.set_ylabel(f"Number of expanded {earlier} clonotypes")
    ax.set_title(
        f"Fate of expanded {earlier} clonotypes\n"
        "This panel asks whether pre-existing large clones persist after measles",
        fontsize=11.5,
        pad=10,
    )
    finish_axis(ax, "y")

    total_count = sum(counts)
    ymax = max(counts) if counts else 0
    for pos, value in zip(positions, counts):
        pct = 100 * value / total_count if total_count else 0.0
        ax.text(
            pos,
            value + 0.025 * max(1, ymax),
            f"{value}\n{pct:.1f}% of clones",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )

    if not fate.empty:
        total_early = float(fate[f"{earlier}_fraction"].sum())
        shares = []
        for category in categories:
            category_fraction = float(
                fate.loc[fate["fate"].eq(category), f"{earlier}_fraction"].sum()
            )
            pct = 100 * category_fraction / total_early if total_early > 0 else 0.0
            shares.append(pct)

        ax.text(
            0.02,
            0.96,
            "Share of the early expanded repertoire by read fraction\n"
            f"retained {shares[0]:.1f}%   below threshold {shares[1]:.1f}%   "
            f"not detected {shares[2]:.1f}%",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8.2,
            bbox={
                "facecolor": "white",
                "edgecolor": "0.8",
                "alpha": 0.94,
                "boxstyle": "round,pad=0.3",
            },
        )


def short_cdr3(value):
    text = str(value)
    if len(text) <= 20:
        return text
    return text[:9] + "..." + text[-6:]


def shared_expanded_table(patient, earlier, later):
    x, y, sx, sy = expanded_sets(patient, earlier, later)
    shared = sx & sy
    rows = []

    for key in shared:
        fa = float(x.at[key, "fraction"])
        fb = float(y.at[key, "fraction"])
        rows.append({
            "cdr3aa": key[0],
            "v_gene": key[1],
            f"{earlier}_fraction": fa,
            f"{later}_fraction": fb,
            "mean_fraction": (fa + fb) / 2.0,
        })

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows).sort_values("mean_fraction", ascending=False)


def draw_shared_expanded(ax, shared, earlier, later):
    if shared.empty:
        ax.axis("off")
        ax.set_title("Top shared expanded clonotypes")
        ax.text(
            0.5,
            0.5,
            "No expanded clonotypes are shared between these points",
            ha="center",
            va="center",
        )
        return

    top = shared.head(10).iloc[::-1].copy()
    y = np.arange(len(top))
    height = 0.36

    ax.barh(
        y - height / 2,
        100 * top[f"{earlier}_fraction"].to_numpy(float),
        height=height,
        label=POINT_LABEL[earlier], color=TIME_COLORS["T1"],
    )
    ax.barh(
        y + height / 2,
        100 * top[f"{later}_fraction"].to_numpy(float),
        height=height,
        label=POINT_LABEL[later], color=TIME_COLORS["T2"],
    )

    labels = [
        f"{short_cdr3(row.cdr3aa)}  {row.v_gene}"
        for row in top.itertuples(index=False)
    ]
    ax.set_yticks(y, labels)
    ax.tick_params(axis="y", labelsize=7.8)
    ax.set_xlabel("Repertoire fraction, percent of productive TRB reads")
    ax.set_title(
        "Top 10 shared expanded clonotypes\n"
        "Ranked by mean repertoire fraction across the two points",
        fontsize=11.5,
        pad=10,
    )
    ax.legend(fontsize=8, loc="lower right")
    finish_axis(ax, "x")


def trbv_delta_table(patient, earlier, later):
    x = point_table(patient, earlier).reset_index(drop=True)
    y = point_table(patient, later).reset_index(drop=True)

    xa = (
        x.groupby("v_gene", as_index=False)["fraction"]
        .sum()
        .rename(columns={"fraction": f"{earlier}_fraction"})
    )
    yb = (
        y.groupby("v_gene", as_index=False)["fraction"]
        .sum()
        .rename(columns={"fraction": f"{later}_fraction"})
    )

    z = xa.merge(yb, on="v_gene", how="outer").fillna(0.0)
    z["delta"] = z[f"{later}_fraction"] - z[f"{earlier}_fraction"]
    z["abs_delta"] = z["delta"].abs()
    return z.sort_values("abs_delta", ascending=False)


def draw_trbv_delta(ax, trbv, earlier, later):
    if trbv.empty:
        ax.axis("off")
        ax.set_title("TRBV repertoire change")
        ax.text(0.5, 0.5, "No TRBV data available", ha="center", va="center")
        return

    top = trbv.head(10).sort_values("delta")
    y = np.arange(len(top))
    values = 100 * top["delta"].to_numpy(float)

    ax.barh(y, values, color=np.where(values >= 0, COLORS["teal"], COLORS["orange"]))
    ax.axvline(0, color=COLORS["ink"], linewidth=0.9)
    ax.set_yticks(y, top["v_gene"].astype(str).tolist())
    ax.tick_params(axis="y", labelsize=8.5)
    ax.set_xlabel(
        f"Change in productive TRB read fraction from {earlier} to {later}, percentage points"
    )
    ax.set_title(
        "Top 10 TRBV changes\n"
        "Positive values indicate higher representation at the later point",
        fontsize=11.5,
        pad=10,
    )
    finish_axis(ax, "x")


valid = metrics[metrics["disease_status"].eq("measles")].copy()
patient_ids = sorted(valid["patient_id"].unique())
summary_rows = []

for patient in patient_ids:
    points = sorted(
        valid.loc[
            valid["patient_id"].eq(patient),
            "clinical_timepoint",
        ].unique(),
        key=lambda point: POINT_ORDER.get(point, 99),
    )
    pair = choose_pair(points)
    if pair is None:
        continue

    earlier, later = pair

    fate = fate_table(patient, earlier, later)
    shared = shared_expanded_table(patient, earlier, later)
    trbv = trbv_delta_table(patient, earlier, later)

    fate.to_csv(
        PATOUT / f"{patient}_{earlier}_to_{later}_expanded_fates.csv",
        index=False,
    )
    shared.to_csv(
        PATOUT / f"{patient}_{earlier}_to_{later}_shared_expanded.csv",
        index=False,
    )
    trbv.to_csv(
        PATOUT / f"{patient}_{earlier}_to_{later}_trbv_delta.csv",
        index=False,
    )

    fig = plt.figure(figsize=(19, 16), constrained_layout=True)
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[1.35, 1.0, 1.0],
        width_ratios=[1.0, 1.0],
    )

    ax_summary = fig.add_subplot(grid[0, 0])
    ax_all = fig.add_subplot(grid[0, 1])
    ax_expanded = fig.add_subplot(grid[1, 0])
    ax_fate = fig.add_subplot(grid[1, 1])
    ax_shared = fig.add_subplot(grid[2, 0])
    ax_trbv = fig.add_subplot(grid[2, 1])

    draw_summary(ax_summary, patient, earlier, later)
    draw_all_clone_venn(ax_all, patient, earlier, later)
    x, y, sx, sy = draw_expanded_venn(ax_expanded, patient, earlier, later)
    draw_fate(ax_fate, fate, earlier)
    draw_shared_expanded(ax_shared, shared, earlier, later)
    draw_trbv_delta(ax_trbv, trbv, earlier, later)

    age = valid.loc[valid["patient_id"].eq(patient), "age_group"].iloc[0]
    age_text = (
        "child younger than 18 years"
        if age == "child"
        else "adult 18 years and older"
    )

    fig.suptitle(
        f"Patient {patient}: TRB repertoire change from {POINT_LABEL[earlier]} "
        f"to {POINT_LABEL[later]}, {age_text}",
        fontsize=17,
        fontweight="bold",
    )

    fig.text(
        0.5,
        -0.005,
        "Expanded clonotype definition: aaV clone with at least 0.10% of productive TRB reads. "
        "aaV means CDR3 beta amino acid sequence together with TRBV gene. "
        "Complete Venn shows all productive aaV clonotypes; expanded Venn is the primary overlap view for immune amnesia.",
        ha="center",
        va="top",
        fontsize=9,
        wrap=True,
    )

    png_path = PATFIG / f"{patient}_expanded_repertoire.png"
    pdf_path = PATFIG / f"{patient}_expanded_repertoire.pdf"
    fig.savefig(png_path, dpi=220, bbox_inches="tight", pad_inches=0.25)
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)

    all_x = set(x.index)
    all_y = set(y.index)
    expanded_union = sx | sy

    summary_rows.append({
        "patient_id": patient,
        "age_group": age,
        "comparison": f"{earlier} to {later}",
        "n_all_clones_earlier": len(all_x),
        "n_all_clones_later": len(all_y),
        "n_shared_all_clones": len(all_x & all_y),
        "n_expanded_earlier": len(sx),
        "n_expanded_later": len(sy),
        "n_shared_expanded": len(sx & sy),
        "expanded_jaccard": (
            len(sx & sy) / len(expanded_union)
            if expanded_union
            else np.nan
        ),
        "figure": f"figures/current/patients/{patient}_expanded_repertoire.png",
    })

index = pd.DataFrame(summary_rows)
index.to_csv(CURRENT / "patient_panel_index.csv", index=False)

print(f"Patient panels generated: {len(index)}")
print("Each panel includes complete repertoire overlap and expanded clonotype overlap")
print("Expanded clone threshold: 0.10 percent of productive TRB reads")
