# Measles TCR beta repertoire analysis

TRB repertoire analysis of the measles cohort, with paired clone tracking, depth-controlled rarefaction and observed-frequency comparisons.

See [RUN_GUIDE.md](RUN_GUIDE.md) for the analysis order and cluster commands, [SOURCE_NOTES.md](config/SOURCE_NOTES.md) for cohort assignments and [TZ_COVERAGE.md](TZ_COVERAGE.md) for specification coverage and remaining work.

## Cohort

- metadata: `config/analysis_metadata.csv`
- final library QC: at least 4400 productive TRB clonotypes
- children: younger than 18 years
- adults: 18 years and older
- controls: `CONTROL` only
- measles time points: T1, T2, T3
- accepted paired cohort: 22 T1-T2 patients, 6 children and 16 adults
- `new_106` and `new_111` are accepted T2 biological repeats for БА1_071 and are normalized separately before averaging

## Clone-level specification

Physical tracking uses CDR3 nucleotide sequence, TRBV and TRBJ. Functional comparisons and public status use CDR3 amino acid sequence, TRBV and TRBJ.

The final master table has one row per T1 nucleotide clonotype.

The primary expansion threshold is 0.10% of productive TRB reads. Sensitivity analyses use 0.05%, 0.10% and every 0.10 percentage point from 0.20% through 1.50%.

Each T1-expanded clone is classified at T2 as retained expanded, detected below the expansion threshold or not detected.

The master table also reports T1 and T2 counts and fractions, absolute change, log2 T2 over T1 for nonzero T2, fold contraction, contribution to the original T1 expanded mass, sampling robustness, public or private status, control prevalence, available antigen annotation, sequence cluster, T3 fate and candidate class.

A zero T2 value never receives an arbitrary pseudocount. A depth-based upper detection bound is retained.

## Branch 1: pair-specific repeated rarefaction

The primary branch controls sampling depth within each patient. Pair depth is the minimum productive-read depth among all accepted T1 and T2 libraries of that patient.

Each contributing library is sampled without replacement to that pair-specific depth. The default is 500 iterations and at least 500 are required.

The random subset changes between iterations, affecting the frequency, detection and threshold classification of rare or near-threshold clones. Results are summarized by medians and 95% sampling intervals with seed 1031.

Expanded-set Jaccard measures the fraction of clones shared by the two expanded sets:

$$
J = \frac{|E_{T1} \cap E_{T2}|}{|E_{T1} \cup E_{T2}|}.
$$

It decreases with loss of T1 clones or appearance of new T2 clones. The accompanying T1 mass not detected at T2 measures the change in abundance.

For БА1_071 the two T2 libraries are rarefied separately, converted to frequencies and only then averaged.

An internal T1 null produces `P_detection_if_unchanged`: the exact same-depth detection probability calculated after each T1 rarefaction and averaged across iterations. A T1-expanded clone with zero observed T2 reads and expected detection probability at least 0.99 is flagged as a sampling-robust loss candidate.

This controls finite sequencing depth but cannot remove PCR amplification bias because these libraries have no UMI.

## Branch 2: no-downsampling Simpson analysis

A complementary branch uses the observed normalized repertoire without read downsampling.

The patient table reports Simpson concentration $D = \sum_i p_i^2$, Simpson clonality $\sqrt{D}$ and normalized richness: unique productive clonotypes divided by productive TRB reads.

The main figure shows only Simpson D because `sqrt(Simpson D)` has exactly the same ordering and would be visually redundant. Normalized richness remains in the supplementary patient table because it answers a different question. The workflow also records the mass of clonotypes above 1 percent as context, while the clone-fate analysis keeps the 0.10 percent primary threshold and applies the full 0.05 to 1.50 percent sensitivity grid.

Every analysis is also rendered separately under `figures/adults` and `figures/children`; these subgroup figures never mix patients from the two age strata.

The same T1-to-T2 fate, direction, candidate, enrichment and specificity logic is repeated in this branch without downsampling.

Repeated rarefaction controls unequal read depth within each T1 and T2 pair. Simpson D measures repertoire concentration using every observed read. It supports interpretation of clone fate, while depth imbalance and the biological meaning of absence require the sampling analysis.

## Clone interpretation

Primary T1-expanded clones are described as:

- increased: detected at T2 with T2/T1 at least 2; candidate developing or persisting measles response
- stable: detected at T2 with T2/T1 between 0.5 and 2; persistent circulating clone, potentially pre-existing memory or ongoing immune surveillance
- decreased: detected at T2 with T2/T1 at most 0.5; either normal contraction of a measles response or possible loss of pre-existing memory
- not detected: zero observed T2 reads; strongest disappearance category, interpreted together with the same-depth sampling null

Specificity is classified as measles-associated, known non-measles memory or other known specificity, or unknown.

The non-measles class is the most relevant for the immune-amnesia hypothesis.

TRBV, TRBJ and CDR3 3-mer features are tested among affected versus retained clones using patient-aware GEE adjusted for baseline T1 clone frequency.

T3 is used only where available. Identity QC uses abundance-aware repertoire similarity. ГА1_003 is repeated as a cross-batch sensitivity case because its T1 is legacy and T2 is new.

## Final figures

The final figure directories are:

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

Only PNG figures are published. Temporary and duplicate exploratory figures are removed.

The patient folder retains the Venn-based patient panels. All final plots use the same typography and color system. The validation folder contains a branch-concordance matrix, an explicit 0.10 percent threshold-choice figure and a rarefaction-iteration convergence figure.

## Final tables

The main result folder contains the concise candidate table:

```text
results/final/candidate_clones.csv
```

The main table is restricted to accepted patients with both T1 and T2 and to T1-expanded candidate clonotypes. Full master, patient summary, threshold sensitivity, feature enrichment, validation and iteration-level convergence tables are placed in `results/supplementary_tables`. Large reusable computational tables are moved to `results/cache`. MiXCR clone exports remain in `results/clonotypes`. Temporary `results/current` is deleted after successful publication.

## Run

```bash
cd ~/measles_2026
git pull --ff-only
eval "$(conda shell.bash hook)"
conda activate measles-tcr-2026
python scripts/validate_repository.py
jid=$(sbatch --parsable scripts/run_analysis_refresh_slurm.sh)
echo "$jid"
tail -F "logs/refresh_slurm_${jid}.out"
```

For a final 1000-iteration rarefaction run:

```bash
export MEASLES_RAREFACTION_ITERATIONS=1000
jid=$(sbatch --parsable scripts/run_analysis_refresh_slurm.sh)
```
