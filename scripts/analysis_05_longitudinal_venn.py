from pathlib import Path
from itertools import combinations
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
LONG_RESULTS = RESULTS / "longitudinal"
PATIENT_RESULTS = LONG_RESULTS / "patients"
PATIENT_FIGURES = FIGURES / "patients"
for p in [LONG_RESULTS, PATIENT_RESULTS, PATIENT_FIGURES]:
    p.mkdir(parents=True, exist_ok=True)

clones_path = RESULTS / "combined_trb_clonotypes.csv"
if not clones_path.exists():
    raise FileNotFoundError("Run analysis preparation first: combined_trb_clonotypes.csv is missing.")

clones = pd.read_csv(clones_path)
clones["analysis_timepoint"] = pd.to_numeric(clones["analysis_timepoint"], errors="coerce")
clones["count"] = pd.to_numeric(clones["count"], errors="coerce").fillna(0)
clones["is_control"] = clones["is_control"].astype(str).str.lower().isin(["true", "1", "yes"])
clones = clones[clones["count"] > 0].copy()

# Only clinical time points .1/.2/.3 are used in the patient longitudinal analysis.
long = clones[clones["analysis_timepoint"].isin([1, 2, 3])].copy()


def safe_name(text):
    return str(text).replace("/", "_").replace("\\", "_").replace(" ", "_")


def save_figure(fig, stem):
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def style_venn(v):
    """Use the same simple visual contract as the mouse-project Venn figures."""
    if v is None:
        return
    for patch in v.patches:
        if patch is not None:
            patch.set_alpha(0.55)
            patch.set_edgecolor("black")
            patch.set_linewidth(1.2)
    for label in v.set_labels or []:
        if label is not None:
            label.set_fontsize(10)
            label.set_fontweight("bold")
    for label in v.subset_labels or []:
        if label is not None:
            label.set_fontsize(9)
            label.set_bbox(dict(facecolor="white", edgecolor="none", alpha=0.75, boxstyle="round,pad=0.15"))


# -----------------------------------------------------------------------------
# 1. Patient Venn diagrams
# -----------------------------------------------------------------------------
# The mouse project used aaV features: amino-acid CDR3 + V gene. We keep the
# same feature definition here for direct comparability. A CDR3-only overlap
# table is also written because it is easy to interpret biologically.
overlap_cdr3_rows = []
overlap_aav_rows = []
patients_with_followup = []

for patient, pdf in long.groupby("patient_id"):
    timepoints = sorted(int(x) for x in pdf["analysis_timepoint"].dropna().unique())
    if len(timepoints) < 2:
        continue
    patients_with_followup.append(patient)

    cdr3_sets = {
        tp: set(pdf.loc[pdf.analysis_timepoint.eq(tp), "cdr3aa"].dropna().astype(str))
        for tp in timepoints
    }
    aav_sets = {}
    for tp in timepoints:
        x = pdf[pdf.analysis_timepoint.eq(tp) & pdf["cdr3aa"].notna() & pdf["v_gene"].notna()]
        aav_sets[tp] = set(zip(x["cdr3aa"].astype(str), x["v_gene"].astype(str)))

    for a, b in combinations(timepoints, 2):
        for sets, target in [(cdr3_sets, overlap_cdr3_rows), (aav_sets, overlap_aav_rows)]:
            sa, sb = sets[a], sets[b]
            union = sa | sb
            target.append({
                "patient_id": patient,
                "timepoint_a": a,
                "timepoint_b": b,
                "n_a": len(sa),
                "n_b": len(sb),
                "shared": len(sa & sb),
                "unique_a": len(sa - sb),
                "unique_b": len(sb - sa),
                "jaccard": len(sa & sb) / len(union) if union else np.nan,
            })

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    if len(timepoints) == 2:
        a, b = timepoints
        v = venn2([aav_sets[a], aav_sets[b]], set_labels=(f"Time point {a}", f"Time point {b}"), ax=ax)
    else:
        chosen = [tp for tp in [1, 2, 3] if tp in aav_sets][:3]
        if len(chosen) < 3:
            chosen = timepoints[:3]
        v = venn3(
            [aav_sets[chosen[0]], aav_sets[chosen[1]], aav_sets[chosen[2]]],
            set_labels=tuple(f"Time point {tp}" for tp in chosen),
            ax=ax,
        )
    style_venn(v)
    ax.set_title(f"{patient}: TRB aaV overlap")
    save_figure(fig, PATIENT_FIGURES / f"{safe_name(patient)}_aav_venn")

pd.DataFrame(overlap_cdr3_rows).to_csv(LONG_RESULTS / "patient_cdr3_overlap.csv", index=False)
pd.DataFrame(overlap_aav_rows).to_csv(LONG_RESULTS / "patient_aav_overlap.csv", index=False)


# -----------------------------------------------------------------------------
# 2. TRBV usage: fraction of all productive TRB reads
# -----------------------------------------------------------------------------
# dropna=False retains healthy CGIMU controls that have no clinical time point.
sample_meta = (
    clones.groupby(["sample_id", "patient_id", "stage", "analysis_timepoint"], as_index=False, dropna=False)
    .agg(
        age_group=("age_group", "first"),
        disease_status=("disease_status", "first"),
        is_control=("is_control", "first"),
        total_reads=("count", "sum"),
    )
)

v_counts = (
    clones[clones["v_gene"].notna()]
    .groupby(["sample_id", "patient_id", "stage", "analysis_timepoint", "v_gene"], as_index=False, dropna=False)["count"]
    .sum()
)
v_counts = v_counts.merge(
    sample_meta[["sample_id", "age_group", "disease_status", "is_control", "total_reads"]],
    on="sample_id",
    how="left",
)
v_counts["read_fraction"] = v_counts["count"] / v_counts["total_reads"]
v_counts.to_csv(LONG_RESULTS / "trbv_read_fractions.csv", index=False)

v_matrix = v_counts.pivot_table(
    index="sample_id", columns="v_gene", values="read_fraction", aggfunc="sum", fill_value=0.0
)
all_genes = sorted(v_matrix.columns)

change_rows = []
transition_rows = []

for patient, pdf in long.groupby("patient_id"):
    timepoints = sorted(int(x) for x in pdf["analysis_timepoint"].dropna().unique())
    if len(timepoints) < 2:
        continue

    age = str(pdf["age_group"].dropna().iloc[0]) if pdf["age_group"].notna().any() else "unknown"
    disease = str(pdf["disease_status"].dropna().iloc[0]) if pdf["disease_status"].notna().any() else "unknown"

    # Compare each observed point with the next observed point for this patient.
    for old_tp, new_tp in zip(timepoints[:-1], timepoints[1:]):
        old_sample = pdf.loc[pdf.analysis_timepoint.eq(old_tp), "sample_id"].iloc[0]
        new_sample = pdf.loc[pdf.analysis_timepoint.eq(new_tp), "sample_id"].iloc[0]

        local = []
        for gene in all_genes:
            old_fraction = float(v_matrix.at[old_sample, gene]) if old_sample in v_matrix.index else 0.0
            new_fraction = float(v_matrix.at[new_sample, gene]) if new_sample in v_matrix.index else 0.0
            delta = new_fraction - old_fraction

            if old_fraction == 0 and new_fraction > 0:
                presence_change = "gained"
            elif old_fraction > 0 and new_fraction == 0:
                presence_change = "lost"
            elif old_fraction > 0 and new_fraction > 0:
                presence_change = "persistent"
            else:
                presence_change = "absent_both"

            direction = "increased" if delta > 0 else "decreased" if delta < 0 else "unchanged"
            row = {
                "patient_id": patient,
                "age_group": age,
                "disease_status": disease,
                "from_timepoint": old_tp,
                "to_timepoint": new_tp,
                "trbv": gene,
                "fraction_old": old_fraction,
                "fraction_new": new_fraction,
                "delta_fraction": delta,
                "delta_percentage_points": 100 * delta,
                "presence_change": presence_change,
                "direction": direction,
            }
            change_rows.append(row)
            local.append(row)

        local_df = pd.DataFrame(local)
        gained = local_df.loc[local_df.presence_change.eq("gained"), "trbv"].tolist()
        lost = local_df.loc[local_df.presence_change.eq("lost"), "trbv"].tolist()
        top_up = local_df.sort_values("delta_fraction", ascending=False).head(5)["trbv"].tolist()
        top_down = local_df.sort_values("delta_fraction", ascending=True).head(5)["trbv"].tolist()

        transition_rows.append({
            "patient_id": patient,
            "age_group": age,
            "disease_status": disease,
            "from_timepoint": old_tp,
            "to_timepoint": new_tp,
            "n_gained_trbv": len(gained),
            "n_lost_trbv": len(lost),
            "gained_trbv": ";".join(gained),
            "lost_trbv": ";".join(lost),
            "top_increased_trbv": ";".join(top_up),
            "top_decreased_trbv": ";".join(top_down),
        })

        patient_table = local_df.sort_values("delta_fraction", key=lambda x: x.abs(), ascending=False)
        patient_table.to_csv(
            PATIENT_RESULTS / f"{safe_name(patient)}_T{old_tp}_to_T{new_tp}_trbv_changes.csv",
            index=False,
        )

        plot_df = patient_table.head(12).sort_values("delta_percentage_points")
        fig, ax = plt.subplots(figsize=(8.0, max(4.5, 0.36 * len(plot_df) + 1.7)))
        ax.barh(plot_df["trbv"], plot_df["delta_percentage_points"])
        ax.axvline(0, linewidth=1, color="black")
        ax.set_xlabel("Change in share of productive TRB reads, percentage points")
        ax.set_ylabel("TRBV")
        ax.set_title(f"{patient}: TRBV change, time point {old_tp} → {new_tp}")
        save_figure(fig, PATIENT_FIGURES / f"{safe_name(patient)}_T{old_tp}_to_T{new_tp}_trbv_change")

changes = pd.DataFrame(change_rows)
transitions = pd.DataFrame(transition_rows)
changes.to_csv(LONG_RESULTS / "trbv_patient_changes.csv", index=False)
transitions.to_csv(LONG_RESULTS / "trbv_patient_transition_summary.csv", index=False)


# -----------------------------------------------------------------------------
# 3. Group statistics and recurrent patterns
# -----------------------------------------------------------------------------
# Individual patient files contain observed effects. Statistical significance is
# evaluated across independent biological patients.
stats_rows = []


def adjust_block(block):
    idx = [i for i, row in enumerate(block) if pd.notna(row.get("p_value"))]
    if idx:
        fdr = multipletests([block[i]["p_value"] for i in idx], method="fdr_bh")[1]
        for i, q in zip(idx, fdr):
            block[i]["fdr"] = float(q)
    return block


def vector_for(samples, gene):
    return np.asarray([
        float(v_matrix.at[s, gene]) if s in v_matrix.index else 0.0
        for s in samples
    ])


# A. Healthy controls vs acute measles, separately for adults and children.
for age in ["adult", "child"]:
    healthy_samples = sample_meta[(sample_meta.age_group == age) & (sample_meta.stage == "healthy")]["sample_id"].tolist()
    acute_samples = sample_meta[(sample_meta.age_group == age) & (sample_meta.stage == "acute")]["sample_id"].tolist()
    block = []
    for gene in all_genes:
        healthy = vector_for(healthy_samples, gene)
        acute = vector_for(acute_samples, gene)
        p_value = float(mannwhitneyu(acute, healthy, alternative="two-sided").pvalue) if len(healthy) >= 3 and len(acute) >= 3 else np.nan
        block.append({
            "comparison": f"{age}: acute vs healthy",
            "age_group": age,
            "comparison_type": "case_control",
            "trbv": gene,
            "n_1": len(healthy),
            "n_2": len(acute),
            "effect_delta_fraction": float(np.median(acute) - np.median(healthy)) if len(healthy) and len(acute) else np.nan,
            "p_value": p_value,
            "fdr": np.nan,
        })
    stats_rows.extend(adjust_block(block))

# B. Adults vs children within the same clinical status.
for stage in ["healthy", "acute"]:
    adult_samples = sample_meta[(sample_meta.age_group == "adult") & (sample_meta.stage == stage)]["sample_id"].tolist()
    child_samples = sample_meta[(sample_meta.age_group == "child") & (sample_meta.stage == stage)]["sample_id"].tolist()
    block = []
    for gene in all_genes:
        adult = vector_for(adult_samples, gene)
        child = vector_for(child_samples, gene)
        p_value = float(mannwhitneyu(adult, child, alternative="two-sided").pvalue) if len(adult) >= 3 and len(child) >= 3 else np.nan
        block.append({
            "comparison": f"{stage}: adult vs child",
            "age_group": "adult_vs_child",
            "comparison_type": "age_within_status",
            "trbv": gene,
            "n_1": len(adult),
            "n_2": len(child),
            "effect_delta_fraction": float(np.median(adult) - np.median(child)) if len(adult) and len(child) else np.nan,
            "p_value": p_value,
            "fdr": np.nan,
        })
    stats_rows.extend(adjust_block(block))

# C. Paired longitudinal changes separately in adults and children.
long_meta = sample_meta[(~sample_meta.is_control.astype(bool)) & sample_meta.analysis_timepoint.isin([1, 2, 3])].copy()
for age in ["adult", "child"]:
    age_meta = long_meta[long_meta.age_group == age]
    for old_tp, new_tp in [(1, 2), (2, 3)]:
        old_map = age_meta[age_meta.analysis_timepoint.eq(old_tp)].set_index("patient_id")["sample_id"].to_dict()
        new_map = age_meta[age_meta.analysis_timepoint.eq(new_tp)].set_index("patient_id")["sample_id"].to_dict()
        patients = sorted(set(old_map) & set(new_map))
        block = []
        for gene in all_genes:
            old = np.asarray([float(v_matrix.at[old_map[p], gene]) if old_map[p] in v_matrix.index else 0.0 for p in patients])
            new = np.asarray([float(v_matrix.at[new_map[p], gene]) if new_map[p] in v_matrix.index else 0.0 for p in patients])
            p_value = np.nan
            if len(patients) >= 3:
                diff = new - old
                if np.allclose(diff, 0):
                    p_value = 1.0
                else:
                    try:
                        p_value = float(wilcoxon(new, old, zero_method="wilcox").pvalue)
                    except ValueError:
                        p_value = np.nan
            block.append({
                "comparison": f"{age}: T{old_tp} → T{new_tp}",
                "age_group": age,
                "comparison_type": "paired_longitudinal",
                "trbv": gene,
                "n_1": len(patients),
                "n_2": len(patients),
                "effect_delta_fraction": float(np.median(new - old)) if len(patients) else np.nan,
                "p_value": p_value,
                "fdr": np.nan,
            })
        stats_rows.extend(adjust_block(block))

# D. Do adults and children recover differently? Compare patient-level deltas.
if not changes.empty:
    for old_tp, new_tp in [(1, 2), (2, 3)]:
        trans = changes[(changes.from_timepoint == old_tp) & (changes.to_timepoint == new_tp)]
        block = []
        for gene in all_genes:
            g = trans[trans.trbv == gene]
            adult = g.loc[g.age_group == "adult", "delta_fraction"].to_numpy(dtype=float)
            child = g.loc[g.age_group == "child", "delta_fraction"].to_numpy(dtype=float)
            p_value = float(mannwhitneyu(adult, child, alternative="two-sided").pvalue) if len(adult) >= 3 and len(child) >= 3 else np.nan
            block.append({
                "comparison": f"T{old_tp}→T{new_tp} delta: adult vs child",
                "age_group": "adult_vs_child",
                "comparison_type": "age_longitudinal_delta",
                "trbv": gene,
                "n_1": len(adult),
                "n_2": len(child),
                "effect_delta_fraction": float(np.median(adult) - np.median(child)) if len(adult) and len(child) else np.nan,
                "p_value": p_value,
                "fdr": np.nan,
            })
        stats_rows.extend(adjust_block(block))

stats = pd.DataFrame(stats_rows)
stats.to_csv(LONG_RESULTS / "trbv_group_statistics.csv", index=False)

# Compact list of strongest recurrent effects for reporting.
pattern_summary = (
    stats[stats["effect_delta_fraction"].notna()]
    .assign(abs_effect=lambda x: x["effect_delta_fraction"].abs())
    .sort_values(["comparison", "abs_effect"], ascending=[True, False])
    .groupby("comparison", as_index=False, group_keys=False)
    .head(10)
    .drop(columns="abs_effect")
)
pattern_summary.to_csv(LONG_RESULTS / "trbv_top_patterns.csv", index=False)

# Group figure 1: strongest TRBV effects across all requested comparisons.
if not stats.empty:
    effect_table = stats.pivot_table(index="comparison", columns="trbv", values="effect_delta_fraction", aggfunc="first")
    if not effect_table.empty:
        top_genes = effect_table.abs().max(axis=0).sort_values(ascending=False).head(20).index
        fig, ax = plt.subplots(figsize=(12, max(5.0, 0.58 * len(effect_table) + 2)))
        sns.heatmap(
            effect_table[top_genes] * 100,
            center=0,
            cmap="vlag",
            ax=ax,
            cbar_kws={"label": "Median difference, percentage points"},
        )
        ax.set_xlabel("TRBV")
        ax.set_ylabel("")
        ax.set_title("Largest TRBV read-share effects across study comparisons")
        save_figure(fig, FIGURES / "09_group_trbv_patterns")

# Group figure 2: number of TRBV genes appearing/disappearing between time points.
if not transitions.empty:
    plot_counts = transitions[transitions.age_group.isin(["adult", "child"])].copy()
    plot_counts["transition"] = (
        "T" + plot_counts.from_timepoint.astype(int).astype(str)
        + "→T" + plot_counts.to_timepoint.astype(int).astype(str)
    )
    melted = plot_counts.melt(
        id_vars=["patient_id", "age_group", "transition"],
        value_vars=["n_gained_trbv", "n_lost_trbv"],
        var_name="change_type",
        value_name="n_trbv",
    )
    melted["change_type"] = melted["change_type"].map({"n_gained_trbv": "Gained", "n_lost_trbv": "Lost"})
    melted["group"] = melted["age_group"].str.capitalize() + ", " + melted["transition"]
    fig, ax = plt.subplots(figsize=(9, 5.2))
    sns.barplot(data=melted, x="group", y="n_trbv", hue="change_type", errorbar="se", ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("Number of TRBV genes per patient")
    ax.set_title("TRBV genes gained or lost between time points")
    ax.tick_params(axis="x", rotation=25)
    save_figure(fig, FIGURES / "10_group_trbv_gain_loss")

print("Patients with at least two time points:", len(patients_with_followup))
print("Patient aaV Venn figures:", len(patients_with_followup))
print("TRBV transition tables:", len(transitions))
print("Group TRBV tests/effects:", len(stats))
