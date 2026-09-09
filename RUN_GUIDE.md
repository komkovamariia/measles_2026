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
- final sensitivity thresholds: 0.05, 0.10, then 0.20 through 1.50 percent in 0.10 percentage-point steps

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
analysis_16_extended_thresholds.py
analysis_15_tz_final_reporting.py
analysis_17_apply_extended_thresholds.py
```

The first established scripts are retained to build metadata, patient Venn panels and reproducible intermediate quantities. Their exploratory figures are temporary and are deleted from the finished tree.

`analysis_10_clone_fate_v2.py` builds the one-row-per-T1-nucleotide-clonotype master scaffolding.

`analysis_12_exact_rarefaction.py` is the primary sampling-controlled branch. Each patient pair is repeatedly rarefied to that pair's own minimum accepted T1/T2 library depth. At least 500 iterations are required. The complete threshold grid is evaluated inside these same rarefaction iterations, so the expensive sampling is not repeated separately for each threshold.

`analysis_11_clone_enrichment.py` performs patient-aware TRBV, TRBJ and CDR3-motif enrichment after sampling correction.

`analysis_14_simpson_no_downsampling.py` implements the complementary no-downsampling branch:

```text
Simpson D = sum(p_i squared)
Simpson clonality = sqrt(Simpson D)
normalized richness = unique productive clonotypes / productive TRB reads
```

Both analytical branches use the same exact expansion grid:

```text
0.05
0.10
0.20
0.30
0.40
0.50
0.60
0.70
0.80
0.90
1.00
1.10
1.20
1.30
1.40
1.50 percent
```

`analysis_16_extended_thresholds.py` is now only a lightweight validator and consolidator. It combines the already computed full-grid outputs from the rarefied and no-downsampling branches and does not run an additional rarefaction loop.

`analysis_15_tz_final_reporting.py` curates the final tree, creates the complete threshold and iteration figures, and moves supporting tables out of the main result folder. `analysis_17_apply_extended_thresholds.py` validates the published grid and figures without redrawing them.

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

The internal T1 null calculates exact same-depth `P_detection_if_unchanged` after each rarefied T1 draw and averages these probabilities across iterations. It does not add a second Bernoulli simulation. A T1-expanded clone absent from real T2 despite very high expected detectability is a sampling-robust loss candidate.

No pseudocount is introduced for T2 zero values.

## Presentation logic

Patient Venn panels are retained as the first layer.

The final rarefied branch then shows:

1. changes in representation inside the shared T1-expanded segment and overall T1 clone fate
2. sampling robustness and threshold sensitivity across 0.05 to 1.50 percent
3. top increased, stable and decreased nonzero candidate clonotypes with TRBV labels
4. patient-clustered TRBV/TRBJ/CDR3 enrichment with the method stated in the title

The no-downsampling branch repeats the same biological logic. Its main figure shows Simpson D only; the monotonic square-root transform and normalized richness remain in the supplementary patient table.

Validation keeps identity QC, available T3 checks, GA1_003 cross-batch sensitivity, explicit branch concordance, the 0.10 percent threshold choice and rarefaction convergence across iterations.

## Final filesystem

After a successful job the temporary directories `results/current` and `figures/current` do not exist.

```text
figures/
  adults/
    patients/
    rarefied/
    simpson_raw/
    validation/
  children/
    patients/
    rarefied/
    simpson_raw/
    validation/
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
  supplementary_tables/
```

The main result folder contains exactly one specification table:

```text
results/final/candidate_clones.csv
```

The full clone master, patient summary, complete threshold sensitivity, enrichment, validation and iteration convergence tables are stored in `results/supplementary_tables`.

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
find results/supplementary_tables -maxdepth 1 -type f | sort
```

Expected Slurm state is `COMPLETED` with exit code `0:0`.
