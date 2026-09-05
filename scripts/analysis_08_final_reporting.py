#!/usr/bin/env python3
"""Final group reporting focused on expanded clonotypes and longitudinal overlap summaries.

This script is self contained for downstream use. If patient level repertoire tables are
missing it rebuilds them from the accepted MiXCR clonotype exports by running
analysis_01_prepare.py first.

The summary overlap figures average patient level quantities. Each patient therefore
contributes one observation regardless of sequencing depth or repertoire size.
"""
from pathlib import Path
import os
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FIG = ROOT / "figures" / "current"
FIG.mkdir(parents=True, exist_ok=True)
CURRENT.mkdir(parents=True, exist_ok=True)

EXPANDED_FRACTION = 0.001


def ensure_patient_level_inputs():
    required = [
        RESULTS / "combined_trb_clonotypes.csv",
        RESULTS / "diversity_metrics.csv",
        CURRENT / "sample_metadata_current.csv",
    ]
    missing = [p for p in required if not p.exists()]
    if missing:
        print("Patient level repertoire tables are missing")
        print("Rebuilding them with analysis_01_prepare.py before final reporting")
        subprocess.run([sys.executable, str(ROOT / "scripts" / "analysis_01_prepare.py")], check=True)
    still_missing = [p for p in required if not p.exists()]
    if still_missing:
        names = ", ".join(str(p.relative_to(ROOT)) for p in still_missing)
        raise FileNotFoundError(f"Required patient level tables were not created: {names}")


ensure_patient_level_inputs()

clones = pd.read_csv(RESULTS / "combined_trb_clonotypes.csv")
metrics = pd.read_csv(RESULTS / "diversity_metrics.csv")
meta = pd.read_csv(CURRENT / "sample_metadata_current.csv")

clones["fraction"] = pd.to_numeric(clones["fraction"], errors="coerce").fillna(0.0)
clones["count"] = pd.to_numeric(clones["count"], errors="coerce").fillna(0.0)
clones["expanded"] = clones["fraction"] >= EXPANDED_FRACTION

# Expanded aaV clone table used by all downstream overlap work
expanded = clones[clones.expanded].copy()
expanded.to_csv(CURRENT / "expanded_clonotypes.csv", index=False)


def gene_summary(gene_col):
    x = expanded.groupby(
        ["sample_id", "patient_id", "clinical_timepoint", "age_group", gene_col],
        as_index=False,
    )["fraction"].sum()
    rows = []
    for age in ["child", "adult"]:
        for point in ["CONTROL", "T1", "T2"]:
            info = metrics[(metrics.age_group == age) & metrics.clinical_timepoint.eq(point)]
            sample_ids = info.sample_id.tolist()
            if not sample_ids:
                continue
            z = x[x.sample_id.isin(sample_ids)]
            genes = sorted(g for g in expanded[gene_col].dropna().astype(str).unique() if g)
            for gene in genes:
                values = (
                    z.loc[z[gene_col].eq(gene)]
                    .set_index("sample_id")["fraction"]
                    .reindex(sample_ids, fill_value=0.0)
                )
                rows.append({
                    "age_group": age,
                    "clinical_group": point,
                    "gene": gene,
                    "n_patients": len(sample_ids),
                    "n_patients_with_expanded_gene": int((values > 0).sum()),
                    "mean_repertoire_share_percent": 100 * float(values.mean()),
                    "median_repertoire_share_percent": 100 * float(values.median()),
                })
    return pd.DataFrame(rows)


trbv = gene_summary("v_gene")
trbv.to_csv(CURRENT / "expanded_TRBV_usage.csv", index=False)

(CURRENT / "expanded_gene_reporting_note.txt").write_text(
    "Expanded gene usage is reported for TRBV because the patient level aaV definition is CDR3 beta amino acid sequence plus TRBV. "
    "TRBJ can be added from library level tables if required, but it is not used in the expanded overlap definition.\n",
    encoding="utf-8",
)


def heatmap(age, filename):
    z = trbv[trbv.age_group.eq(age)].copy()
    if z.empty:
        return
    group_order = [g for g in ["CONTROL", "T1", "T2"] if g in set(z.clinical_group)]
    score = z.groupby("gene")["mean_repertoire_share_percent"].mean().sort_values(ascending=False)
    genes = score.head(12).index.tolist()
    mat = np.zeros((len(genes), len(group_order)))
    ann = [["" for _ in group_order] for _ in genes]
    for i, gene in enumerate(genes):
        for j, group in enumerate(group_order):
            r = z[(z.gene == gene) & (z.clinical_group == group)]
            if r.empty:
                continue
            rr = r.iloc[0]
            mat[i, j] = rr.mean_repertoire_share_percent
            ann[i][j] = f"{int(rr.n_patients_with_expanded_gene)} of {int(rr.n_patients)}\nmean {rr.mean_repertoire_share_percent:.2f}%"
    fig, ax = plt.subplots(figsize=(9.0, max(6.2, 0.52 * len(genes) + 2.5)))
    im = ax.imshow(mat, aspect="auto")
    labels = {"CONTROL": "Healthy controls", "T1": "T1 early infection", "T2": "T2 recovery"}
    ax.set_xticks(np.arange(len(group_order)), [labels[g] for g in group_order])
    ax.set_yticks(np.arange(len(genes)), genes)
    ax.set_ylabel("TRBV gene")
    age_title = "Children younger than 18 years" if age == "child" else "Adults 18 years and older"
    ax.set_title(f"{age_title}: TRBV usage among expanded clonotypes")
    vmax = float(np.nanmax(mat)) if mat.size else 0.0
    for i in range(len(genes)):
        for j in range(len(group_order)):
            color = "white" if vmax and mat[i, j] > 0.55 * vmax else "black"
            ax.text(j, i, ann[i][j], ha="center", va="center", fontsize=8, color=color)
    fig.colorbar(im, ax=ax, label="Mean share of complete productive TRB repertoire percent")
    fig.text(0.5, 0.015, "Only clonotypes expanded to at least 0.10 percent are included", ha="center", fontsize=8.5)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(FIG / f"{filename}.png", dpi=210, bbox_inches="tight")
    fig.savefig(FIG / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)


heatmap("child", "07_expanded_TRBV_children")
heatmap("adult", "08_expanded_TRBV_adults")

# Specificity reference coverage figure so absence of matches is not confused with
# biological absence of a clone
ref = RESULTS / "antigen_reference_used.csv"
if ref.exists():
    r = pd.read_csv(ref)
    if "pathogen" in r.columns:
        coverage = r.groupby("pathogen", as_index=False).agg(reference_records=("cdr3_trimmed", "nunique"))
        coverage.to_csv(CURRENT / "specificity_reference_coverage.csv", index=False)
        fig, ax = plt.subplots(figsize=(7.6, 5.2))
        ax.bar(coverage.pathogen.astype(str), coverage.reference_records.to_numpy())
        ax.set_ylabel("Unique CDR3 sequences in reference")
        ax.set_title("Specificity reference coverage used for exact matching")
        ax.tick_params(axis="x", rotation=25)
        fig.tight_layout()
        fig.savefig(FIG / "09_specificity_reference_coverage.png", dpi=210, bbox_inches="tight")
        fig.savefig(FIG / "09_specificity_reference_coverage.pdf", bbox_inches="tight")
        plt.close(fig)

# -----------------------------------------------------------------------------
# Patient averaged overlap summaries for the paired T1 and T2 cohort
# -----------------------------------------------------------------------------

def point_table(patient, point):
    z = clones[
        clones["patient_id"].eq(patient)
        & clones["clinical_timepoint"].eq(point)
        & clones["disease_status"].eq("measles")
    ].copy()
    z["clone_key"] = list(zip(z["cdr3aa"].astype(str), z["v_gene"].astype(str)))
    return z.set_index("clone_key", drop=False)


pair_rows = []
case_metrics = metrics[metrics["disease_status"].eq("measles")].copy()
for age in ["child", "adult"]:
    z = case_metrics[case_metrics.age_group.eq(age)]
    t1 = set(z.loc[z.clinical_timepoint.eq("T1"), "patient_id"])
    t2 = set(z.loc[z.clinical_timepoint.eq("T2"), "patient_id"])
    patients = sorted(t1 & t2)
    for patient in patients:
        a = point_table(patient, "T1")
        b = point_table(patient, "T2")
        all1 = set(a.index)
        all2 = set(b.index)
        all_union = all1 | all2
        all_inter = all1 & all2

        exp1 = set(a.index[a["expanded"]])
        exp2 = set(b.index[b["expanded"]])
        exp_union = exp1 | exp2
        exp_inter = exp1 & exp2

        retained = 0
        below = 0
        lost = 0
        retained_mass = 0.0
        below_mass = 0.0
        lost_mass = 0.0
        t1_expanded_mass = float(a.loc[list(exp1), "fraction"].sum()) if exp1 else 0.0

        for key in exp1:
            t1f = float(a.at[key, "fraction"])
            if key in exp2:
                retained += 1
                retained_mass += t1f
            elif key in b.index:
                below += 1
                below_mass += t1f
            else:
                lost += 1
                lost_mass += t1f

        n_t1_exp = len(exp1)
        pair_rows.append({
            "patient_id": patient,
            "age_group": age,
            "all_T1_only_count": len(all1 - all2),
            "all_shared_count": len(all_inter),
            "all_T2_only_count": len(all2 - all1),
            "all_T1_only_union_percent": 100 * len(all1 - all2) / len(all_union) if all_union else np.nan,
            "all_shared_union_percent": 100 * len(all_inter) / len(all_union) if all_union else np.nan,
            "all_T2_only_union_percent": 100 * len(all2 - all1) / len(all_union) if all_union else np.nan,
            "all_jaccard": len(all_inter) / len(all_union) if all_union else np.nan,
            "expanded_T1_only_count": len(exp1 - exp2),
            "expanded_shared_count": len(exp_inter),
            "expanded_T2_only_count": len(exp2 - exp1),
            "expanded_T1_only_union_percent": 100 * len(exp1 - exp2) / len(exp_union) if exp_union else np.nan,
            "expanded_shared_union_percent": 100 * len(exp_inter) / len(exp_union) if exp_union else np.nan,
            "expanded_T2_only_union_percent": 100 * len(exp2 - exp1) / len(exp_union) if exp_union else np.nan,
            "expanded_jaccard": len(exp_inter) / len(exp_union) if exp_union else np.nan,
            "retained_T1_expanded_count_percent": 100 * retained / n_t1_exp if n_t1_exp else np.nan,
            "below_T1_expanded_count_percent": 100 * below / n_t1_exp if n_t1_exp else np.nan,
            "lost_T1_expanded_count_percent": 100 * lost / n_t1_exp if n_t1_exp else np.nan,
            "retained_T1_expanded_read_percent": 100 * retained_mass / t1_expanded_mass if t1_expanded_mass else np.nan,
            "below_T1_expanded_read_percent": 100 * below_mass / t1_expanded_mass if t1_expanded_mass else np.nan,
            "lost_T1_expanded_read_percent": 100 * lost_mass / t1_expanded_mass if t1_expanded_mass else np.nan,
        })

pair_overlap = pd.DataFrame(pair_rows)
pair_overlap.to_csv(CURRENT / "patient_venn_overlap_summary.csv", index=False)

AGE_LABEL = {
    "child": "Children younger than 18 years",
    "adult": "Adults 18 years and older",
}


def mean_with_patient_dots(ax, data, columns, labels, title, ylabel):
    x = np.arange(len(columns), dtype=float)
    offsets = {"child": -0.16, "adult": 0.16}
    width = 0.28
    for age in ["child", "adult"]:
        z = data[data.age_group.eq(age)]
        means = [float(pd.to_numeric(z[c], errors="coerce").mean()) for c in columns]
        ax.bar(x + offsets[age], means, width=width, alpha=0.72, label=f"{AGE_LABEL[age]}  n {len(z)}")
        for j, c in enumerate(columns):
            vals = pd.to_numeric(z[c], errors="coerce").dropna().to_numpy(float)
            if len(vals):
                jitter = np.linspace(-0.045, 0.045, len(vals)) if len(vals) > 1 else np.array([0.0])
                ax.scatter(np.full(len(vals), x[j] + offsets[age]) + jitter, vals, s=23, alpha=0.72)
    ax.set_xticks(x, labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(fontsize=8)


if not pair_overlap.empty:
    # Complete repertoire Venn composition averaged at patient level
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    mean_with_patient_dots(
        ax,
        pair_overlap,
        ["all_T1_only_union_percent", "all_shared_union_percent", "all_T2_only_union_percent"],
        ["T1 only", "Shared", "T2 only"],
        "Complete productive aaV repertoire overlap across paired patients",
        "Percent of each patient's T1 and T2 clonotype union",
    )
    ax.set_ylim(bottom=0)
    fig.text(
        0.5, 0.015,
        "Percentages are calculated within each patient before averaging, so every patient has equal weight",
        ha="center", fontsize=8.5,
    )
    fig.tight_layout(rect=[0, 0.045, 1, 1])
    fig.savefig(FIG / "10_complete_clonotype_overlap_summary.png", dpi=210, bbox_inches="tight")
    fig.savefig(FIG / "10_complete_clonotype_overlap_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    # Expanded repertoire Venn composition averaged at patient level
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    mean_with_patient_dots(
        ax,
        pair_overlap,
        ["expanded_T1_only_union_percent", "expanded_shared_union_percent", "expanded_T2_only_union_percent"],
        ["T1 only expanded", "Shared expanded", "T2 only expanded"],
        "Expanded clonotype overlap across paired patients",
        "Percent of each patient's expanded clonotype union",
    )
    ax.set_ylim(bottom=0)
    fig.text(
        0.5, 0.015,
        "Expanded clonotype means aaV clone representing at least 0.10 percent of productive TRB reads. Each dot is one patient",
        ha="center", fontsize=8.5,
    )
    fig.tight_layout(rect=[0, 0.045, 1, 1])
    fig.savefig(FIG / "11_expanded_clonotype_overlap_summary.png", dpi=210, bbox_inches="tight")
    fig.savefig(FIG / "11_expanded_clonotype_overlap_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fate of T1 expanded clones summarized both by clone count and by T1 read mass
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.6), constrained_layout=True)
    mean_with_patient_dots(
        axes[0],
        pair_overlap,
        ["retained_T1_expanded_count_percent", "below_T1_expanded_count_percent", "lost_T1_expanded_count_percent"],
        ["Retained expanded", "Below threshold", "Not detected at T2"],
        "Fate of T1 expanded clones by clone count",
        "Percent of expanded T1 clonotypes",
    )
    mean_with_patient_dots(
        axes[1],
        pair_overlap,
        ["retained_T1_expanded_read_percent", "below_T1_expanded_read_percent", "lost_T1_expanded_read_percent"],
        ["Retained expanded", "Below threshold", "Not detected at T2"],
        "Fate of T1 expanded clones by T1 read representation",
        "Percent of T1 expanded repertoire read fraction",
    )
    for ax in axes:
        ax.set_ylim(0, 100)
        ax.tick_params(axis="x", rotation=12)
    fig.suptitle(
        "Patient averaged fate of clonotypes expanded during early measles infection",
        fontsize=14, fontweight="bold",
    )
    fig.savefig(FIG / "12_expanded_T1_fate_summary.png", dpi=210, bbox_inches="tight")
    fig.savefig(FIG / "12_expanded_T1_fate_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    # Complete and expanded Jaccard similarity
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.9), constrained_layout=True)
    for ax, column, title in [
        (axes[0], "all_jaccard", "Complete repertoire Jaccard similarity"),
        (axes[1], "expanded_jaccard", "Expanded repertoire Jaccard similarity"),
    ]:
        for xpos, age in enumerate(["child", "adult"], start=1):
            vals = pd.to_numeric(pair_overlap.loc[pair_overlap.age_group.eq(age), column], errors="coerce").dropna().to_numpy(float)
            jitter = np.linspace(-0.07, 0.07, len(vals)) if len(vals) > 1 else np.array([0.0])
            ax.scatter(np.full(len(vals), xpos) + jitter, vals, s=38, alpha=0.75)
            if len(vals):
                ax.hlines(np.median(vals), xpos - 0.18, xpos + 0.18, linewidth=2.4)
        ax.set_xticks([1, 2], ["Children", "Adults"])
        ax.set_ylim(0, 1)
        ax.set_ylabel("Jaccard similarity")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Similarity of paired T1 and T2 repertoires across patients", fontsize=14, fontweight="bold")
    fig.savefig(FIG / "13_patient_overlap_similarity_summary.png", dpi=210, bbox_inches="tight")
    fig.savefig(FIG / "13_patient_overlap_similarity_summary.pdf", bbox_inches="tight")
    plt.close(fig)

marker = FIG / "GENERATED_CURRENT.txt"
marker.write_text(
    "Final figures use QC threshold 4400 productive TRB clonotypes, child age below 18 years, adult age 18 years and older, controls labelled CONTROL, and expanded clone threshold 0.10 percent.\n"
    "Patient overlap summaries first calculate percentages within each patient and then average across patients.\n",
    encoding="utf-8",
)

print("Final expanded clone group reporting completed")
print(f"Paired T1 and T2 patients summarized: {len(pair_overlap)}")
print("If patient level repertoire tables were absent they were rebuilt automatically")
