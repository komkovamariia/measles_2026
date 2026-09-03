task1 = antigen[antigen.pathogen.eq("Measles") & antigen.age_group.isin(["adult", "child"])].copy()
task1_cmp = task1[task1.stage.isin(["healthy", "acute"])].copy()
task1_cmp["matched_percent"] = 100 * task1_cmp["matched_clonotype_fraction"]

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, age in zip(axes, ["adult", "child"]):
    sub = task1_cmp[task1_cmp.age_group.eq(age)]
    if sub.empty:
        ax.text(0.5, 0.5, "No samples", ha="center", va="center")
    else:
        order = [x for x in ["healthy", "acute"] if x in sub.stage.unique()]
        sns.boxplot(data=sub, x="stage", y="matched_percent", order=order, ax=ax, showfliers=False)
        sns.stripplot(data=sub, x="stage", y="matched_percent", order=order, ax=ax, size=5)
    ax.set_title(age.capitalize())
    ax.set_xlabel("")
    ax.set_ylabel("Exact IEDB matches (% of unique TRB clonotypes)" if age == "adult" else "")
fig.suptitle("Task 1: exact IEDB measles-associated CDR3β matches")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGURES / f"01_task1_measles_group_comparison.{ext}", dpi=220, bbox_inches="tight")
plt.show()

task1_oligo = diversity[
    diversity.age_group.isin(["adult", "child"]) &
    diversity.stage.isin(["healthy", "acute"])
].copy()
task1_oligo["top10_percent"] = 100 * task1_oligo["rarefied_top10_fraction"]

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, age in zip(axes, ["adult", "child"]):
    sub = task1_oligo[task1_oligo.age_group.eq(age)]
    if sub.empty:
        ax.text(0.5, 0.5, "No samples", ha="center", va="center")
    else:
        order = [x for x in ["healthy", "acute"] if x in sub.stage.unique()]
        sns.boxplot(data=sub, x="stage", y="top10_percent", order=order, ax=ax, showfliers=False)
        sns.stripplot(data=sub, x="stage", y="top10_percent", order=order, ax=ax, size=5)
    ax.set_title(age.capitalize())
    ax.set_xlabel("")
    ax.set_ylabel("Reads represented by top 10 clones (%)" if age == "adult" else "")
fig.suptitle("Task 1: oligoclonal dominance of the TRB repertoire")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGURES / f"02_task1_oligoclonal_dominance.{ext}", dpi=220, bbox_inches="tight")
plt.show()

other_names = ["Hepatitis B", "Poliovirus", "Tetanus"]
task2 = antigen[antigen.pathogen.isin(other_names) & antigen.age_group.isin(["adult", "child"])].copy()
task2_cmp = task2[task2.stage.isin(["healthy", "acute"])].copy()
task2_cmp["matched_percent"] = 100 * task2_cmp["matched_clonotype_fraction"]

g = sns.catplot(
    data=task2_cmp,
    x="pathogen", y="matched_percent",
    hue="stage", col="age_group",
    col_order=["adult", "child"],
    kind="strip", dodge=True, height=4, aspect=1.15,
)
g.set_axis_labels("", "Exact IEDB matches (% of unique TRB clonotypes)")
g.set_titles("{col_name}")
g.fig.suptitle("Task 2: exact IEDB matches to common infections", y=1.05)
for ext in ["png", "pdf"]:
    g.fig.savefig(FIGURES / f"03_task2_other_infections_group_comparison.{ext}", dpi=220, bbox_inches="tight")
plt.show()

# Fraction of unique clonotypes matching at least one of HBV, poliovirus or tetanus.
other_match_rows = matched_unique[matched_unique.pathogen.isin(other_names)].copy()
other_union_counts = (
    other_match_rows.drop_duplicates(["sample_id", "cdr3aa"])
    .groupby("sample_id")["cdr3aa"].nunique()
    .rename("other_matched_clonotypes")
)
other_total = sample_meta.merge(other_union_counts, on="sample_id", how="left")
other_total["other_matched_clonotypes"] = other_total["other_matched_clonotypes"].fillna(0)
other_total["other_infection_percent"] = 100 * other_total["other_matched_clonotypes"] / other_total["total_clonotypes"]
other_total = other_total[
    other_total.age_group.isin(["adult", "child"]) &
    other_total.stage.isin(["acute", "recovery", "follow_up"])
].copy()
other_total["stage"] = pd.Categorical(other_total["stage"], categories=["acute", "recovery", "follow_up"], ordered=True)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, age in zip(axes, ["adult", "child"]):
    sub = other_total[other_total.age_group.eq(age)].sort_values(["patient_id", "stage"])
    n_long = sub.groupby("patient_id")["stage"].nunique()
    sub = sub[sub.patient_id.isin(n_long[n_long >= 2].index)]
    for patient, pdf in sub.groupby("patient_id"):
        ax.plot(pdf["stage"].astype(str), pdf["other_infection_percent"], marker="o", alpha=0.65)
    if sub.empty:
        ax.text(0.5, 0.5, "No longitudinal pairs", ha="center", va="center")
    ax.set_title(age.capitalize())
    ax.set_xlabel("")
    ax.set_ylabel("HBV/polio/tetanus matches (% of unique TRB clonotypes)" if age == "adult" else "")
fig.suptitle("Task 2: common-infection-associated clonotypes during recovery")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGURES / f"04_task2_other_infections_longitudinal.{ext}", dpi=220, bbox_inches="tight")
plt.show()

div_plot = diversity[diversity.age_group.isin(["adult", "child"]) & diversity.stage.isin(["healthy", "acute", "recovery", "follow_up"])].copy()
stage_order = ["healthy", "acute", "recovery", "follow_up"]

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
sns.boxplot(data=div_plot, x="stage", y="rarefied_richness", hue="age_group", order=stage_order, ax=axes[0], showfliers=False)
sns.stripplot(data=div_plot, x="stage", y="rarefied_richness", hue="age_group", order=stage_order, dodge=True, ax=axes[0], size=4, legend=False)
axes[0].set_title("CDR3β richness")
axes[0].set_xlabel("")
axes[0].set_ylabel(f"Unique clonotypes after downsampling to {rarefaction_depth:,} reads")

sns.boxplot(data=div_plot, x="stage", y="shannon", hue="age_group", order=stage_order, ax=axes[1], showfliers=False)
sns.stripplot(data=div_plot, x="stage", y="shannon", hue="age_group", order=stage_order, dodge=True, ax=axes[1], size=4, legend=False)
axes[1].set_title("Shannon diversity")
axes[1].set_xlabel("")
axes[1].set_ylabel("Shannon index")

handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles[:2], labels[:2], title="Age group")
if axes[1].get_legend() is not None:
    axes[1].get_legend().remove()
fig.suptitle("Task 3: depth-normalized TCRβ repertoire diversity")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGURES / f"05_task3_diversity.{ext}", dpi=220, bbox_inches="tight")
plt.show()

usage_source = clones[clones.age_group.isin(["adult", "child"]) & clones.stage.isin(["healthy", "acute", "recovery", "follow_up"])].copy()
usage_source = usage_source.drop_duplicates(["sample_id", "cdr3aa"])
usage_source["group"] = usage_source["age_group"].str.capitalize() + " " + usage_source["stage"].str.replace("_", " ")
group_order = [
    "Adult healthy", "Adult acute", "Adult recovery", "Adult follow up",
    "Child healthy", "Child acute", "Child recovery", "Child follow up",
]

def mean_gene_usage(df, gene_col, top_n):
    x = df.dropna(subset=[gene_col]).copy()
    per_sample = x.groupby(["sample_id", "group", gene_col]).size().rename("n").reset_index()
    totals = per_sample.groupby("sample_id")["n"].transform("sum")
    per_sample["fraction"] = per_sample["n"] / totals
    mean_usage = per_sample.groupby(["group", gene_col])["fraction"].mean().reset_index()
    top = mean_usage.groupby(gene_col)["fraction"].mean().nlargest(top_n).index
    matrix = mean_usage[mean_usage[gene_col].isin(top)].pivot(index=gene_col, columns="group", values="fraction").fillna(0)
    return matrix.reindex(columns=[g for g in group_order if g in matrix.columns])

v_matrix = mean_gene_usage(usage_source, "v_gene", 15)
j_matrix = mean_gene_usage(usage_source, "j_gene", 12)

v_matrix.to_csv(RESULTS / "v_usage_group_mean.csv")
j_matrix.to_csv(RESULTS / "j_usage_group_mean.csv")

fig, axes = plt.subplots(1, 2, figsize=(15, 7))
if not v_matrix.empty:
    sns.heatmap(v_matrix, cmap="viridis", ax=axes[0])
    axes[0].set_title("Top TRBV usage")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("TRBV gene")
else:
    axes[0].text(0.5, 0.5, "No V-gene data", ha="center", va="center")
if not j_matrix.empty:
    sns.heatmap(j_matrix, cmap="viridis", ax=axes[1])
    axes[1].set_title("Top TRBJ usage")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("TRBJ gene")
else:
    axes[1].text(0.5, 0.5, "No J-gene data", ha="center", va="center")
fig.suptitle("Task 3: unique-clonotype V/J usage")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGURES / f"06_task3_vj_usage.{ext}", dpi=220, bbox_inches="tight")
plt.show()

long_div = diversity[
    diversity.age_group.isin(["adult", "child"]) &
    diversity.stage.isin(["acute", "recovery", "follow_up"])
].copy()
long_div["stage"] = pd.Categorical(long_div["stage"], categories=["acute", "recovery", "follow_up"], ordered=True)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
for ax, age in zip(axes, ["adult", "child"]):
    sub = long_div[long_div.age_group.eq(age)].sort_values(["patient_id", "stage"])
    n_long = sub.groupby("patient_id")["stage"].nunique()
    sub = sub[sub.patient_id.isin(n_long[n_long >= 2].index)]
    for patient, pdf in sub.groupby("patient_id"):
        ax.plot(pdf["stage"].astype(str), pdf["rarefied_richness"], marker="o", alpha=0.65)
    if sub.empty:
        ax.text(0.5, 0.5, "No longitudinal pairs", ha="center", va="center")
    ax.set_title(age.capitalize())
    ax.set_xlabel("")
    ax.set_ylabel("Depth-normalized CDR3β richness")
fig.suptitle("Task 4: within-patient repertoire richness")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGURES / f"07_task4_richness_longitudinal.{ext}", dpi=220, bbox_inches="tight")
plt.show()

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
for ax, age in zip(axes, ["adult", "child"]):
    sub = long_div[long_div.age_group.eq(age)].sort_values(["patient_id", "stage"])
    n_long = sub.groupby("patient_id")["stage"].nunique()
    sub = sub[sub.patient_id.isin(n_long[n_long >= 2].index)]
    for patient, pdf in sub.groupby("patient_id"):
        ax.plot(pdf["stage"].astype(str), pdf["shannon"], marker="o", alpha=0.65)
    if sub.empty:
        ax.text(0.5, 0.5, "No longitudinal pairs", ha="center", va="center")
    ax.set_title(age.capitalize())
    ax.set_xlabel("")
    ax.set_ylabel("Shannon diversity")
fig.suptitle("Task 4: within-patient Shannon diversity")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGURES / f"08_task4_shannon_longitudinal.{ext}", dpi=220, bbox_inches="tight")
plt.show()
