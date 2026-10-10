# Coverage of the measles clone fate specification

This checklist maps the 7 September 2026 specification to the reproducible pipeline. It distinguishes completed requirements from conclusions that cannot be produced from the available data alone.

| Specification item | Status | Implementation or limitation |
|---|---|---|
| One row per physical T1 clone using CDR3 nucleotide, TRBV and TRBJ | Complete | `paired_T1_clone_master_nt.csv`; curated candidates in `results/final/candidate_clones.csv` |
| Functional and public clone identity using CDR3 amino acid, TRBV and TRBJ | Complete | Public/private status and control prevalence are calculated at the amino-acid key |
| T1 and T2 counts, frequencies, age and sequencing source | Complete | The concise candidate table contains the T1 and T2 percentages requested in its example schema. Counts, age and source/batch fields remain in `supplementary_tables/01_full_clone_master.csv` rather than duplicating the master table. |
| Primary T1-expanded pool at 0.10% | Complete | The same pre-specified pool is used in both analytical branches |
| Threshold sensitivity | Complete and extended | 0.05%, 0.10% and every 0.10 percentage point from 0.20% through 1.50% |
| Fate as retained expanded, detected below threshold or not detected | Complete | Both clone counts and original T1 expanded mass are summarized |
| Absolute change, log2 ratio, fold contraction and T1 mass share | Complete | Zero T2 is not assigned an artificial finite fold change |
| Detection limit for zero T2 without pseudocount | Complete | A depth-based 95% upper detection bound and fold-contraction lower bound are retained |
| Pair-specific repeated rarefaction | Complete | Multivariate hypergeometric sampling without replacement; at least 500 iterations; fixed seed 1031 |
| Internal unchanged-frequency detection control | Complete | Exact same-depth detection probability is calculated after each T1 rarefaction and averaged, without an extra Bernoulli simulation |
| TRBV, TRBJ and CDR3 motif enrichment | Complete | Patient-clustered binomial GEE, adjusted for baseline T1 fraction and BH multiple testing |
| Measles, known non-measles and unknown specificity classes | Complete where annotation exists | Exact antigen annotation is data-dependent; unknown remains an explicit category |
| T3 interpretation check | Complete where T3 exists | Reappearance or continued absence is retained as `T3_fate` and included in validation |
| All-versus-all repertoire identity fingerprint | Complete | Abundance-aware Morisita similarity matrix and identity validation output |
| Separate GA1_003 cross-batch sensitivity | Complete | Retained in the validation analysis |
| Separate adult-only and child-only outputs | Complete | Parallel `figures/adults` and `figures/children` trees; subgroup enrichment is refitted rather than copied |
| HLA compatibility | Not estimable | Patient HLA typing is unavailable; the field is omitted from the concise final table rather than filled with an uninformative value |
| Candidate-by-candidate literature synthesis | Not complete | Reproducible search strings are prepared, but manual evidence review and source adjudication are still required |

## Recommended additional figures

Two additional figures are proposed for the paper.

1. Patient-level forest plot of T1-expanded mass not detected at T2. Show one estimate and the 2.5–97.5% rarefaction interval per patient, sorted within adults and children. This would expose heterogeneity, identify influential patients and keep the patient as the biological unit.
2. Compact candidate-evidence plot for the strongest contracted or lost clones. Combine T1 and T2 fraction, fold contraction or its lower bound, same-depth detection probability, specificity class and available T3 fate. This would connect statistical evidence to the final candidate table more directly than another repertoire-wide summary.

The iteration-convergence figure should remain a validation or supplementary figure. Its left panel tests stability of the abundance-weighted primary result; its right panel tests stability of expanded-set membership through Jaccard. Neither panel is intended as a biological group comparison.
