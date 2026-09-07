# Measles TCR beta repertoire analysis

This repository contains the paired TRB repertoire workflow for the measles cohort.

## Cohort contract

- metadata: `config/analysis_metadata.csv`
- final library QC: at least 4400 productive TRB clonotypes
- children: younger than 18 years
- adults: 18 years and older
- controls: `CONTROL` only
- measles time points: T1, T2, T3
- accepted paired cohort: 22 T1-T2 patients, 6 children and 16 adults
- `new_106` and `new_111` are accepted T2 biological repeats for БА1_071 and are normalized separately before averaging

## Clone-level specification

Primary physical tracking unit:

```text
CDR3 nucleotide + TRBV + TRBJ
```

Functional and public analysis unit:

```text
CDR3 amino acid + TRBV + TRBJ
```

The final master table has one row per T1 nucleotide clonotype.

Primary T1 expanded threshold:

```text
0.10 percent of productive TRB reads
```

Final sensitivity thresholds:

```text
0.05 percent
0.10 percent
0.20 percent
0.50 percent
```

For every primary T1-expanded clone the workflow distinguishes:

```text
retained expanded
detected below the expansion threshold
not detected at T2
```

It also reports T1/T2 counts and fractions, absolute change, log2 T2 over T1 for nonzero T2, fold contraction, contribution to the original T1 expanded mass, sampling robustness, public/private status, control prevalence, available antigen annotation, motif cluster, T3 fate and candidate class.

A zero T2 value never receives an arbitrary pseudocount. A depth-based upper detection bound is retained.

## Branch 1: pair-specific repeated rarefaction

This is the primary sampling-controlled branch.

For each patient:

```text
pair depth = minimum productive-read depth among all accepted T1 and T2 libraries of that patient
```

Each contributing library is sampled without replacement to that pair-specific depth. The default is 500 iterations and at least 500 are required.

For БА1_071 the two T2 libraries are rarefied separately, converted to frequencies and only then averaged.

An internal T1 null produces empirical `P_detection_if_unchanged`. A high-frequency T1 clone that should nearly always be observed at the available depth but is absent from real T2 is flagged as a sampling-robust loss candidate.

This controls finite sequencing depth but cannot remove PCR amplification bias because these libraries have no UMI.

## Branch 2: no-downsampling Simpson analysis

A complementary branch uses the observed normalized repertoire without read downsampling.

Metrics follow PMID 37051229 / PMCID PMC10084475:

```text
Simpson D = sum(p_i squared)
Simpson clonality = sqrt(Simpson D)
normalized richness = unique productive clonotypes / total productive TRB reads
```

The paper also emphasizes top-expanded clonotypes above 1 percent. The workflow records the mass of such clones as article-style context, while the clone-fate analysis itself keeps the 0.10 percent primary threshold with 0.05, 0.20 and 0.50 percent sensitivity checks.

The same T1-to-T2 fate, direction, candidate, enrichment and specificity logic is repeated in this branch without downsampling.

## Biological interpretation used in the presentation

Primary T1-expanded clones are described as:

- increased: candidate developing or persisting measles response
- stable: stable circulating clone, potentially pre-existing memory or ongoing immune surveillance
- decreased: either normal contraction of a measles response or possible loss of pre-existing memory
- not detected: strongest disappearance category, interpreted together with the sampling null

The specificity fork is:

```text
measles-associated
known non-measles memory or other known specificity
unknown
```

The non-measles class is the most relevant for the immune-amnesia hypothesis.

TRBV, TRBJ and CDR3 3-mer features are tested among affected versus retained clones using patient-aware GEE adjusted for baseline T1 clone frequency.

T3 is used only where available. Identity QC uses abundance-aware repertoire similarity. ГА1_003 is repeated as a cross-batch sensitivity case because its T1 is legacy and T2 is new.

## Final figures

The finished run deliberately publishes a small final tree:

```text
figures/
  patients/
  rarefied/
  simpson_raw/
  validation/
  FIGURE_INDEX.txt
```

Only PNG figures are published. Temporary and duplicate exploratory figures are removed.

The patient folder retains the Venn-based patient panels. The rarefied and Simpson folders contain only figures that map directly to the presentation logic: shared-clone changes, clone fate, sampling or threshold robustness, candidate dumbbells, TRBV context, enrichment and specificity.

## Final tables

The finished run deliberately keeps only six presentation-level tables:

```text
results/final/01_clone_master.csv
results/final/02_candidates.csv
results/final/03_patient_summary.csv
results/final/04_threshold_sensitivity.csv
results/final/05_feature_enrichment.csv
results/final/06_validation.csv
```

Large reusable computational tables are moved to `results/cache`. MiXCR clone exports remain in `results/clonotypes`. Temporary `results/current` is deleted after successful publication.

The candidate table contains a PubMed-ready `literature_query` for each primary candidate. Exact antigen annotations are retained when available. The actual literature review is intentionally performed after the final candidate list is produced because top contracted clones can change after sampling correction and threshold sensitivity analysis.

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
