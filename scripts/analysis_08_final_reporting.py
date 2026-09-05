#!/usr/bin/env python3
"""Final group reporting focused on expanded clonotypes

TRBV and TRBJ summaries use only clonotypes that reach the expanded threshold
within a biological patient observation. The denominator remains the complete
productive TRB repertoire so the values show how much of the repertoire is
occupied by expanded clonotypes carrying each gene
"""
from pathlib import Path
import os
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
clones = pd.read_csv(RESULTS / "combined_trb_clonotypes.csv")
clones["fraction"] = pd.to_numeric(clones["fraction"], errors="coerce").fillna(0.0)
clones["expanded"] = clones["fraction"] >= EXPANDED_FRACTION
metrics = pd.read_csv(RESULTS / "diversity_metrics.csv")

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
            for gene in sorted(g for g in expanded[gene_col].dropna().astype(str).unique() if g):
                values = z.loc[z[gene_col].eq(gene)].set_index("sample_id")["fraction"].reindex(sample_ids, fill_value=0.0)
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

# TRBJ annotations are not retained in the combined aaV table, so report TRBV as the
# primary gene result and preserve a clear note instead of silently mixing definitions
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

marker = FIG / "GENERATED_CURRENT.txt"
marker.write_text(
    "Final figures use QC threshold 4400 productive TRB clonotypes, child age below 18 years, adult age 18 years and older, controls labelled CONTROL, and expanded clone threshold 0.10 percent.\n",
    encoding="utf-8",
)

print("Final expanded clone group reporting completed")
