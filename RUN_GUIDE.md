# Analysis run guide

## Current contract

- metadata: `config/analysis_metadata.csv`
- library QC: at least 4400 productive TRB clonotypes
- children: younger than 18 years
- adults: 18 years and older
- controls: `CONTROL` only
- paired T1-T2 cohort: 22 patients, 6 children and 16 adults
- physical clone key: CDR3 nucleotide + TRBV + TRBJ
- functional/public key: CDR3 amino acid + TRBV + TRBJ
- primary expanded threshold: 0.10 percent
- final sensitivity thresholds: 0.05, 0.10, 0.20 and 0.50 percent

`new_106` and `new_111` are accepted T2 biological repeats for БА1_071. They are handled separately through normalization and rarefaction, then averaged into one biological T2 observation.

## Analysis order

```text
analysis_01_prepare.py
analysis_02_antigens.py
analysis_06_clinical_summary.py
analysis_07_longitudinal_expanded.py
analysis_08_final_reporting.py
analysis_09_patient_panels.py
analysis_10_clone_fate_v2.py
analysis_12_exact_rarefaction.py
analysis_11_clone_enrichment.py
analysis_14_simpson_no_downsampling.py
analysis_15_tz_final_reporting.py
```

The first established scripts are retained to build metadata, patient Venn panels and reproducible intermediate quantities. Their exploratory figures are temporary and are deleted from the finished tree.

`analysis_10_clone_fate_v2.py` builds the one-row-per-T1-nucleotide-clonotype master scaffolding.

`analysis_12_exact_rarefaction.py` is the primary sampling-controlled branch. Each patient pair is repeatedly rarefied to that pair's own minimum accepted T1/T2 library depth. At least 500 iterations are required.

`analysis_11_clone_enrichment.py` performs patient-aware TRBV, TRBJ and CDR3-motif enrichment after sampling correction.

`analysis_14_simpson_no_downsampling.py` implements the complementary no-downsampling branch using the definitions from PMID 37051229 / PMCID PMC10084475:

```text
Simpson D = sum(p_i squared)
Simpson clonality = sqrt(Simpson D)
normalized richness = unique productive clonotypes / productive TRB reads
```

The same observed clone-fate and threshold-sensitivity logic is repeated in that branch.

`analysis_15_tz_final_reporting.py` is the only final publisher. It keeps only figures and tables that map to the requested presentation logic and deletes staging output.

## Pair-specific rarefaction

For patient `i`:

```text
pair_depth_i = minimum productive-read depth among all accepted T1 and T2 libraries contributing to that patient
```

Every contributing library is sampled without replacement to `pair_depth_i` using multivariate hypergeometric sampling.

Default:

```text
iterations: 500
seed: 1031
```

Final optional run:

```bash
export MEASLES_RAREFACTION_ITERATIONS=1000
```

The internal T1 null produces empirical `P_detection_if_unchanged`. A T1-expanded clone absent from real T2 despite very high expected detectability is a sampling-robust loss candidate.

No pseudocount is introduced for T2 zero values.

## Presentation logic

Patient Venn panels are retained as the first layer.

The final rarefied branch then shows:

1. changes in representation inside the shared T1-expanded segment and overall T1 clone fate
2. sampling robustness and threshold sensitivity
3. top increased, stable and decreased nonzero candidate clonotypes with TRBV labels
4. TRBV/TRBJ/CDR3 enrichment and the measles-associated versus known non-measles versus unknown specificity fork

The no-downsampling Simpson branch repeats the same biological logic while adding Simpson D, Simpson clonality and normalized richness.

Validation keeps only identity QC, available T3 checks, GA1_003 cross-batch sensitivity and concordance between the two analytical branches.

## Final filesystem

After a successful job the temporary directories `results/current` and `figures/current` do not exist.

```text
figures/
  patients/
  rarefied/
  simpson_raw/
  validation/
  FIGURE_INDEX.txt
```

Only PNG figures are retained.

```text
results/
  clonotypes/
  cache/
  final/
```

Final presentation tables are exactly:

```text
results/final/01_clone_master.csv
results/final/02_candidates.csv
results/final/03_patient_summary.csv
results/final/04_threshold_sensitivity.csv
results/final/05_feature_enrichment.csv
results/final/06_validation.csv
```

The final candidate table contains `literature_query`. This is the handoff for a PubMed pass after the candidate list is fixed. Literature evidence is not invented inside the offline cluster workflow.

## Run

```bash
cd ~/measles_2026
git pull --ff-only
eval "$(conda shell.bash hook)"
conda activate measles-tcr-2026
python scripts/validate_repository.py
jid=$(sbatch --parsable scripts/run_analysis_refresh_slurm.sh)
echo "JOB ID: $jid"
tail -F "logs/refresh_slurm_${jid}.out"
```

After completion:

```bash
sacct -j "$jid" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
find figures -maxdepth 2 -type f | sort
find results/final -maxdepth 1 -type f | sort
```

Expected Slurm state is `COMPLETED` with exit code `0:0`.
