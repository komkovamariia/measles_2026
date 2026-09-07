# Measles TCR beta repertoire analysis

This repository contains the analysis workflow for bulk human TRB sequencing in the measles cohort

The authoritative metadata is `config/analysis_metadata.csv`

## Final cohort rules

- final library quality threshold: at least 4400 productive TRB clonotypes
- child: younger than 18 years
- adult: 18 years and older
- healthy controls: `CONTROL` only
- measles time points: T1, T2 and T3
- T0 is not used
- current accepted paired cohort: 22 T1-T2 patients, 6 children and 16 adults

T1 is early infection after infection and before the mature adaptive response is expected to be fully developed

T2 is immediately after clinical recovery

T3 is the follow-up point about one month after recovery

## Established analyses retained for reproducibility

The previous aaV analyses remain available

`aaV` means CDR3 beta amino-acid sequence together with TRBV gene

The established expanded threshold is 0.10 percent of productive TRB reads

Patient panels retain complete repertoire overlap, expanded repertoire overlap, T1 expanded-clone fate, shared expanded clones and TRBV changes

The old exploratory group figures are not all republished. In particular, old TRBV heatmaps are excluded from the final figure tree because the specification now tests gene enrichment from the clone-level candidate table instead of selecting genes visually

## Primary clone-level analysis

The central data product is

```text
results/current/paired_T1_clone_master_nt.csv
```

One row is one T1 clonotype from one of the 22 accepted T1-T2 patients

Physical tracking key

```text
CDR3 nucleotide + TRBV + TRBJ
```

Functional and public-clone key

```text
CDR3 amino acid + TRBV + TRBJ
```

For each T1 clone the table contains T1 and T2 counts and fractions, age, sequencing batch, expansion status, fate, absolute change, log2 T2 over T1 when T2 is nonzero, contraction, share of the original T1 expanded mass, sampling robustness, public/private annotation, healthy-control prevalence, available antigen annotation, sequence-neighbour cluster, T3 fate and final candidate class

The primary candidate pool contains all T1 clonotypes at or above 0.10 percent

Sensitivity thresholds are

```text
0.05 percent
0.10 percent
0.20 percent
0.50 percent
0.75 percent
1.00 percent
```

Primary T2 fate categories are

```text
retained expanded
detected below threshold
not detected
```

A clone such as 2.0 percent at T1 and 0.03 percent at T2 is therefore represented as a strong contraction rather than collapsed into a simple shared or not-shared label

Zero T2 values do not receive an arbitrary pseudocount. A 95 percent upper detection bound is reported instead

## Pair-specific repeated rarefaction and internal null

The authoritative depth-control analysis follows the longitudinal specification pair by pair

For patient `i`

```text
pair_depth_i = minimum productive-read depth among every accepted T1 and T2 library contributing to that patient
```

Every contributing sequencing library is sampled without replacement to that patient's `pair_depth_i`

The procedure is repeated at least 500 times, with 500 as the default and 1000 supported for the final sensitivity run

Observed-library rarefaction uses multivariate hypergeometric sampling without replacement

For each T1 clonotype an internal unchanged-repertoire null is also evaluated repeatedly. The empirical fraction of iterations in which the clonotype would be detected at the T2 sampling effort is stored as

```text
P_detection_if_unchanged
```

A T1-expanded clone that is absent in real T2 despite empirical `P_detection_if_unchanged >= 0.99` is flagged as sampling-robust not detected

The older analytical raw-depth detection probability is preserved as a separate checking column

This controls finite sequencing depth but does not remove PCR amplification bias because these bulk genomic DNA libraries do not contain UMI information

## БА1_071 biological repeat

`new_106` and `new_111` are accepted T2 biological-repeat libraries for patient БА1_071

They are not counted as two patient observations

During rarefaction each library is sampled separately. Frequencies are normalized within each library and only then averaged to form one T2 biological observation

## Affected-clone enrichment

After pair-specific sampling correction, the enrichment analysis compares

```text
sampling-robust not detected or strong contraction
```

against

```text
retained expanded
```

TRBV, TRBJ and CDR3 amino-acid 3-mer features are tested with patient-aware binomial GEE

Patient ID defines the within-patient dependence and baseline T1 clone frequency is included as a covariate

The final enrichment figure is generated only when at least one feature passes FDR 0.05. Otherwise the complete numerical table is retained without a decorative plot

## Sequence neighbourhoods

Candidate TCRs are also grouped into exploratory amino-acid sequence neighbourhoods within the same TRBV

The current implementation uses CDR3 amino-acid edit distance and is deliberately not labelled TCRdist

A cross-patient neighbourhood figure is generated only if at least two affected clones occur in a cluster spanning at least two patients

The table records possible clonotypic replacement when the exact clone is lost but a nearby T2 sequence remains

HLA compatibility is `not available` unless patient HLA typing is supplied

## T3 and identity validation

Available longitudinal validation is used without inventing missing points

```text
ЭГ1_042: T1 T2 T3
СА1_013: valid T2 T3 because T1 failed QC
ММ1_033: valid T1 T3 because T2 failed identity QC
```

For ЭГ1_042 the workflow explicitly records T1-expanded clones absent in T2 that reappear in T3

All-vs-all Morisita-Horn repertoire similarity is used for sample identity QC

ГА1_003 is summarized both in the complete paired cohort and after excluding this cross-batch T1 legacy to T2 new pair

## Final figure architecture

The final published tree is intentionally shallow

```text
figures/main
figures/patients
figures/downsampling_pairwise
figures/sensitivity_thresholds
figures/validation
figures/supplementary
figures/FIGURE_INDEX.txt
```

Core conclusions are limited to a small number of figures

```text
main/21_primary_fate_pairwise_rarefied.png
downsampling_pairwise/22_sampling_robust_clone_loss.png
sensitivity_thresholds/23_threshold_sensitivity_pairwise_rarefied.png
validation/26_validation_identity_T3_batch.png
```

Two additional main figures are conditional on the data

```text
main/24_lost_vs_retained_enrichment.png
main/25_sequence_neighbourhoods.png
```

Patient panels are kept separately

Supplementary contains only selected established contextual plots. Old TRBV heatmaps and duplicate exploratory plots are intentionally not republished

## Main tables

```text
paired_T1_clone_master_nt.csv
candidate_clones_compact.csv
pairwise_rarefaction_summary.csv
pairwise_threshold_sensitivity.csv
clone_sampling_null.csv
lost_vs_retained_feature_enrichment.csv
sequence_neighbour_cluster_summary.csv
extended_T3_clone_validation.csv
extended_T3_validation_summary.csv
repertoire_identity_morisita_matrix.csv
repertoire_identity_qc.csv
GA1_003_cross_batch_sensitivity.csv
```

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

To use 1000 repeated pair-specific rarefactions

```bash
export MEASLES_RAREFACTION_ITERATIONS=1000
```

Completed MiXCR clonotype exports are reused. The nucleotide analysis expects `nSeqCDR3`. If an existing TSV lacks it, re-export the existing `.clns` file rather than repeating FASTQ alignment and assembly

See `RUN_GUIDE.md` for the detailed execution and verification contract
