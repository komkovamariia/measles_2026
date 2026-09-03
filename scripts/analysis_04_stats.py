def unpaired_test(df, value, age, label):
    sub = df[df.age_group.eq(age)]
    a = sub[sub.stage.eq("healthy")][value].dropna()
    b = sub[sub.stage.eq("acute")][value].dropna()
    result = {"comparison": label, "age_group": age, "metric": value, "n_healthy": len(a), "n_acute": len(b), "p_value": np.nan}
    if len(a) >= 3 and len(b) >= 3:
        result["p_value"] = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    return result

def paired_test(df, value, age, label):
    sub = df[df.age_group.eq(age) & df.stage.isin(["acute", "recovery"])].copy()
    piv = sub.pivot_table(index="patient_id", columns="stage", values=value, aggfunc="first").dropna(subset=["acute", "recovery"])
    result = {"comparison": label, "age_group": age, "metric": value, "n_pairs": len(piv), "p_value": np.nan}
    if len(piv) >= 3:
        try:
            result["p_value"] = float(wilcoxon(piv["acute"], piv["recovery"], alternative="two-sided").pvalue)
        except ValueError:
            result["p_value"] = 1.0
    return result

stats_rows = []
for age in ["adult", "child"]:
    stats_rows.append(unpaired_test(diversity, "rarefied_richness", age, "healthy_vs_acute"))
    stats_rows.append(unpaired_test(diversity, "shannon", age, "healthy_vs_acute"))
    stats_rows.append(unpaired_test(task1, "matched_clonotype_fraction", age, "healthy_vs_acute_measles_matches"))
    for pathogen in other_names:
        stats_rows.append(unpaired_test(task2[task2.pathogen.eq(pathogen)], "matched_clonotype_fraction", age, f"healthy_vs_acute_{pathogen}"))
    stats_rows.append(paired_test(diversity, "rarefied_richness", age, "acute_vs_recovery"))
    stats_rows.append(paired_test(diversity, "shannon", age, "acute_vs_recovery"))
    stats_rows.append(paired_test(task1, "matched_clonotype_fraction", age, "acute_vs_recovery_measles_matches"))

stats_df = pd.DataFrame(stats_rows)
mask = stats_df.p_value.notna()
stats_df["fdr"] = np.nan
if mask.any():
    stats_df.loc[mask, "fdr"] = multipletests(stats_df.loc[mask, "p_value"], method="fdr_bh")[1]
stats_df.to_csv(RESULTS / "statistical_tests.csv", index=False)
stats_df

antigen_wide = antigen.pivot_table(
    index=["sample_id", "patient_id", "age_group", "stage", "analysis_timepoint"],
    columns="pathogen",
    values="matched_clonotype_fraction",
    aggfunc="first",
).reset_index()
antigen_wide.columns.name = None
antigen_wide = antigen_wide.rename(columns={
    "Measles": "measles_match_fraction",
    "Hepatitis B": "hepatitis_b_match_fraction",
    "Poliovirus": "poliovirus_match_fraction",
    "Tetanus": "tetanus_match_fraction",
})

final_metrics = diversity.merge(
    antigen_wide,
    on=["sample_id", "patient_id", "age_group", "stage", "analysis_timepoint"],
    how="left",
)
final_metrics.to_csv(RESULTS / "sample_analysis_metrics.csv", index=False)

longitudinal = final_metrics[final_metrics.stage.isin(["acute", "recovery", "follow_up"])].copy()
counts_per_patient = longitudinal.groupby("patient_id")["stage"].nunique()
longitudinal = longitudinal[longitudinal.patient_id.isin(counts_per_patient[counts_per_patient >= 2].index)]
longitudinal.to_csv(RESULTS / "longitudinal_metrics.csv", index=False)

summary_lines = [
    "Measles TCR-beta analysis",
    f"Biological samples analysed: {final_metrics.sample_id.nunique()}",
    f"Common rarefaction depth: {rarefaction_depth}",
    f"Longitudinal patients with >=2 stages: {longitudinal.patient_id.nunique()}",
    f"Exact IEDB measles-matched sample/clonotype rows: {matched_unique[matched_unique.pathogen.eq('Measles')].shape[0]}",
    "See statistical_tests.csv for simple group and paired tests.",
]
(RESULTS / "result_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
print("\n".join(summary_lines))
