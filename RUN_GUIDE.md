# Analysis run guide

## Analysis settings

- metadata: `config/analysis_metadata.csv`
- library QC: at least 4400 productive TRB clonotypes
- children: younger than 18 years
- adults: 18 years and older
- controls: `CONTROL` only
- paired T1-T2 cohort: 22 patients, 6 children and 16 adults
- physical clone key: CDR3 nucleotide + TRBV + TRBJ
- functional and public key: CDR3 amino acid + TRBV + TRBJ
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

The preparation and early reporting scripts build metadata, patient Venn panels and intermediate quantities. The final reporting step removes their temporary exploratory figures.

`analysis_10_clone_fate_v2.py` builds the master table with one row per T1 nucleotide clonotype.

`analysis_12_exact_rarefaction.py` is the primary sampling-controlled branch. Each patient pair is repeatedly rarefied to that pair's own minimum accepted T1/T2 library depth. At least 500 iterations are required. The complete threshold grid uses the same rarefaction draws, so sampling runs once for all thresholds.

`analysis_11_clone_enrichment.py` performs patient-aware TRBV, TRBJ and CDR3-motif enrichment after sampling correction.

`analysis_14_simpson_no_downsampling.py` implements the complementary no-downsampling branch:

```text
Simpson D = sum(p_i squared)
Simpson clonality = sqrt(Simpson D)
normalized richness = unique productive clonotypes / productive TRB reads
```

Both branches use the expansion grid above: 0.05%, 0.10% and 0.20% through 1.50% in steps of 0.10 percentage points.

`analysis_16_extended_thresholds.py` validates and combines the completed threshold tables from both branches. It reuses those tables without repeating rarefaction.

`analysis_15_tz_final_reporting.py` exports the final tables and figures, including threshold sensitivity and iteration convergence. It moves supporting tables to the supplementary directory. `analysis_17_apply_extended_thresholds.py` checks the published grid and figures without redrawing them.

## Pair-specific rarefaction

For patient `i`:

```text
pair_depth_i = minimum productive-read depth among all accepted T1 and T2 libraries contributing to that patient
```

Every contributing library is sampled without replacement to `pair_depth_i` using multivariate hypergeometric sampling.

Each iteration draws a new random subset of reads. Rare and near-threshold clones can change their sampled count, frequency, detection, expansion status, fate and membership in the T1 and T2 overlap. The input libraries stay fixed.

Repeated draws estimate patient medians and 2.5–97.5% sampling intervals. The default is 500 iterations with seed 1031. Use the convergence figure to assess Monte Carlo precision; 1,000 iterations are available for a sensitivity run.

For expanded sets $E_{T1}$ and $E_{T2}$, Jaccard is:

$$
J = \frac{|E_{T1} \cap E_{T2}|}{|E_{T1} \cup E_{T2}|}.
$$

Jaccard compares clone membership. Its range is 0 to 1: zero means no shared clones and one means identical sets. Loss of a T1-expanded clone and appearance of a new T2-expanded clone both reduce the score. The 0.10% threshold reduces the contribution of unstable rare clones. T1 mass not detected at T2 is reported alongside Jaccard because it also accounts for abundance.

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

Clone direction labels used throughout the figures are:

```text
increased     T2/T1 >= 2
stable        0.5 < T2/T1 < 2
decreased     T2/T1 <= 0.5 and the clone is still detected at T2
not detected  zero observed T2 reads
```

The enrichment endpoint is narrower than the direction classification. Contracted or lost clones have at least five-fold contraction or zero T2 reads despite expected same-depth detection probability at least 0.99. Retained clones are detected at T2 at or above 0.10% and are outside the affected group.

## Figures

Patient Venn panels precede the final rarefied figures, which show:

1. changes in representation inside the shared T1-expanded segment and overall T1 clone fate
2. sampling robustness and threshold sensitivity across 0.05 to 1.50 percent
3. top increased, stable and decreased nonzero candidate clonotypes with TRBV labels
4. patient-clustered enrichment of TRBV, TRBJ and CDR3 motifs with the method stated in the figure

The no-downsampling branch repeats the same biological logic. Its main figure shows Simpson D only; the monotonic square-root transform and normalized richness remain in the supplementary patient table.

Validation keeps identity QC, available T3 checks, GA1_003 cross-batch sensitivity, explicit branch concordance, the 0.10 percent threshold choice and rarefaction convergence across iterations.

## Output directories

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
