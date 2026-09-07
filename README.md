# Measles TCR beta repertoire analysis

This repository contains the final analysis workflow for bulk human TRB sequencing in the measles cohort

The authoritative repository metadata is `config/analysis_metadata.csv`

It is the analysis copy of the final reviewed metadata and includes library quality, final inclusion status, corrected biological patient identifiers, age groups and biological replicate handling

## Final cohort rules

- final library quality threshold: at least 4400 productive TRB clonotypes
- child: younger than 18 years
- adult: 18 years and older
- healthy controls: `CONTROL` only
- measles time points: T1, T2 and T3
- T0 is not used
- T4 is part of the intended clinical design at about six months but is not present among the current sequencing libraries

T1 is the early infection point after infection and before the mature adaptive response is expected to be fully developed

T2 is immediately after clinical recovery

T3 is the follow up point about one month after recovery

## Two complementary established analysis designs

### Cross sectional design

All accepted T1 observations are retained even when the patient has no T2 sample

Healthy controls are compared with all T1 observations and separately with all T2 observations within the same age group

Current biological observation counts after QC are

| Age group | Healthy controls | T1 | T2 |
| --- | ---: | ---: | ---: |
| Children younger than 18 years | 4 | 27 | 8 |
| Adults 18 years and older | 13 | 29 | 19 |

T1 and T2 are not treated as independent groups for a direct independent test because some observations belong to the same patient

### Paired longitudinal design

The primary immune amnesia analysis uses only patients with both accepted T1 and T2 observations

- children younger than 18 years: 6 paired patients
- adults 18 years and older: 16 paired patients

Clonality and Shannon diversity are calculated from the complete patient repertoire

Clone overlap is calculated only among expanded clonotypes

## Established expanded clonotypes

An expanded clonotype is defined as an `aaV` clone with at least 0.10 percent of productive TRB reads in a biological patient observation

`aaV` means CDR3 beta amino acid sequence together with TRBV gene

The expanded clone filter removes the low frequency background from overlap analysis. Bulk genomic DNA TRB sequencing does not directly phenotype naive, memory or effector T cells, therefore the workflow does not label every low frequency clone as naive

The established immune amnesia readouts are

- T1 to T2 change in clonality
- T1 to T2 change in Shannon diversity
- expanded T1 clones that remain expanded at T2
- expanded T1 clones detected at T2 below the expansion threshold
- expanded T1 clones not detected at T2
- expanded clone Jaccard similarity

New expanded clones at T2 are exported separately for exact specificity matching against the available reference database

## Additive nucleotide clone fate analysis

The established aaV workflow above is preserved. A new additive analysis uses all 22 accepted paired T1 and T2 patients and tracks nucleotide clonotypes directly

Primary physical tracking key

```text
CDR3 nucleotide sequence + TRBV + TRBJ
```

Functional and public clone key

```text
CDR3 amino acid sequence + TRBV + TRBJ
```

The primary candidate pool contains all T1 nucleotide clonotypes at or above 0.10 percent of productive TRB reads

Sensitivity is repeated at

```text
0.05 percent
0.10 percent
0.20 percent
0.50 percent
0.75 percent
1.00 percent
```

For every T1 clonotype the new master table reports T1 and T2 counts, normalized fractions, read equivalents, complete fate, absolute change, log2 T2 over T1 when T2 is nonzero, contraction, the fraction of original T1 expanded mass, public or private status, healthy control prevalence, available exact antigen annotation, T3 fate, batch information and candidate class

Zero T2 values do not receive an arbitrary pseudocount. The workflow reports a 95 percent upper detection bound based on the T2 productive read depth and an analytical `P_detection_if_unchanged`

Repeated depth control uses a requested depth of 50000 productive reads and 500 iterations by default. If the shallowest accepted paired library is below the requested depth the workflow automatically uses that lower observed depth rather than upsampling

The repeated sampling step is used for repertoire-level sensitivity summaries. Clone-level detection probability is calculated analytically from the binomial sampling model rather than estimated by unnecessary brute-force iterations

The sampling model controls finite read depth but does not remove PCR amplification bias because these libraries do not have UMI information

The nucleotide analysis also adds

- public and private amino-acid clone annotation
- prevalence among healthy controls
- simple amino-acid sequence-neighbour clusters
- a flag for possible clonotypic replacement when a nearby T2 sequence remains
- patient-aware TRBV, TRBJ and CDR3 motif enrichment among lost versus retained clones
- T3 follow-up of apparent T2 losses
- all-pairs Morisita-Horn repertoire identity fingerprinting
- a sensitivity summary with and without the cross-batch patient ГА1_003

The sequence-neighbour clusters are deliberately not labelled TCRdist clusters. They use CDR3 amino-acid edit distance within the same TRBV and are an exploratory motif analysis

HLA compatibility is reported as not available unless patient HLA typing is supplied

## БА1_071 T2 biological repeat

`new_106` and `new_111` are both accepted T2 libraries for patient БА1_071

`new_111` is the raw sample `111_БА1_071.2а`

The two libraries are not counted as two patients or two independent T2 observations

For patient level analysis each library is first normalized to clonotype frequencies. Missing clones receive frequency zero in the replicate where they are absent. The normalized frequencies are then averaged. The resulting averaged repertoire is one T2 observation for patient БА1_071

## Legacy cohort correction

UDP0097 through UDP0105 are measles samples but UDP is an adapter identifier rather than the biological patient identifier

The corrected biological patients are МР1_001, БА1_002, ГА1_003, ХИ1_004, ИМ1_005, ПИ1_006, СС1_007 and ШЕ1_008

UDP0106 through UDP0112 are healthy controls CGIMU_N1 through CGIMU_N7

UDP0108 remains excluded because its repertoire is unusually similar to UDP0107

The legacy sample ХИ1_004.2 is T2 for the same patient as ХИ1_004.1

The legacy sample ГА1_003.1 is T1 for the same patient as the new sample ГА1_003.2 at T2

## Main outputs

Group figures are published to

```text
figures/group
```

Patient figures are published to

```text
figures/patients
```

Final tables are stored in

```text
results/current
```

Established tables include

```text
cross_sectional_cohort_counts.csv
cross_sectional_statistics.csv
longitudinal_expanded_clone_summary.csv
expanded_clone_fates.csv
new_T2_expansions.csv
expanded_TRBV_usage.csv
patient_panel_index.csv
```

Additive nucleotide clone fate tables include

```text
paired_T1_clone_master_nt.csv
clone_fate_threshold_sensitivity.csv
rarefaction_summary.csv
sequence_neighbour_cluster_summary.csv
lost_vs_retained_feature_enrichment.csv
T3_clone_fate_validation.csv
repertoire_identity_morisita_matrix.csv
repertoire_identity_qc.csv
GA1_003_cross_batch_sensitivity.csv
clone_fate_v2_README.txt
```

Existing patient panels remain unchanged and continue to show complete repertoire overlap, expanded repertoire overlap, T1 expanded clone fate, shared expanded clones and TRBV changes

New compact group figures are added after the established figures rather than replacing them

## Refresh analysis without rerunning completed MiXCR exports

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

The refresh job reuses existing MiXCR clonotype exports. The nucleotide analysis expects the existing MiXCR clone export to contain `nSeqCDR3`. If an old export lacks nucleotide CDR3 it can be re-exported from the existing `.clns` file without repeating FASTQ alignment or assembly

Optional depth-control settings

```bash
export MEASLES_RAREFACTION_DEPTH=50000
export MEASLES_RAREFACTION_ITERATIONS=500
export MEASLES_RAREFACTION_SEED=1031
```

## Full run from raw FASTQ

```bash
bash run_all_from_scratch.sh
```

The full workflow uses MiXCR 4.7.0 with the project environment and the configured human DNA TCR multiplex preset

## Interpretation limits

Read fractions from multiplex genomic DNA TRB sequencing are relative abundance measures and can contain PCR amplification bias

The expanded clone thresholds are operational repertoire thresholds, not direct cell phenotype measurements

A sampling-robust loss means that simple finite read sampling is an unlikely explanation under the unchanged-frequency model. It does not by itself prove destruction of antigen-specific immune memory because trafficking, tissue redistribution, PCR bias and other biological or technical effects remain possible

Exact database matching has incomplete antigen coverage. A missing match must not be interpreted as absence of biological antigen specificity

See `config/SOURCE_NOTES.md` for cohort reconciliation decisions and analysis definitions
