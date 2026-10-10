# Cohort metadata and sample assignments

`config/analysis_metadata.csv` contains the reviewed library quality values, QC decisions, patient identifiers, age groups and biological replicate assignments.

Libraries must contain at least 4400 productive TRB clonotypes. Children are younger than 18 years; adults are 18 years and older. Healthy controls use `CONTROL`. Measles samples use T1, T2 and T3. T0 is absent from the study. T4 was planned at about six months, but the current dataset contains no T4 sequencing library.

## Clinical time points

T1 is early infection, before the mature adaptive response is expected to be fully developed. T2 is immediately after clinical recovery. T3 is about one month after recovery.

For measles libraries the sequencing suffix determines the time point:

| Suffix | Time point |
|---|---|
| `.1` | T1 |
| `.2` or `.2а` | T2 |
| `.3` | T3 |

Healthy controls keep `CONTROL` even when the raw sample label contains a numeric suffix.

## T2 replicates for БА1_071

`new_106`, labelled `106_БА1_071.2`, and `new_111`, labelled `111_БА1_071.2а`, are accepted T2 libraries for patient БА1_071. Both pass the 4400-clonotype threshold and contribute to one biological T2 observation.

Each library is normalized to frequencies summing to one before averaging. A clonotype missing from a replicate contributes zero for that replicate. This gives the libraries equal weight regardless of sequencing depth. The rarefaction branch samples each library separately before averaging its frequencies.

## Legacy sample assignments

UDP labels identify sequencing adapters. The patient assignments are:

| Adapter | Sample | Patient | Time point | Disease status |
|---|---|---|---|---|
| UDP0097 | МР1_001.1 | МР1_001 | T1 | measles |
| UDP0098 | БА1_002.1 | БА1_002 | T1 | measles |
| UDP0099 | ГА1_003.1 | ГА1_003 | T1 | measles |
| UDP0100 | ХИ1_004.1 | ХИ1_004 | T1 | measles |
| UDP0101 | ИМ1_005.1 | ИМ1_005 | T1 | measles |
| UDP0102 | ПИ1_006.1 | ПИ1_006 | T1 | measles |
| UDP0103 | СС1_007.1 | СС1_007 | T1 | measles |
| UDP0104 | ХИ1_004.2 | ХИ1_004 | T2 | measles |
| UDP0105 | ШЕ1_008.1 | ШЕ1_008 | T1 | measles |

UDP0106 through UDP0112 correspond to healthy adult controls CGIMU_N1 through CGIMU_N7. UDP0108 is excluded from biological inference because its repertoire is unusually similar to UDP0107.

ГА1_003.1 belongs to the legacy batch and ГА1_003.2 to the new batch. This patient is retained in the paired analysis and checked separately for batch sensitivity.

## Cross-sectional and paired analyses

Accepted T1 observations are compared with age-matched controls, including patients without T2. Accepted T2 observations are also summarized against controls. Some patients contribute to both time points, so a direct independent-group test of T1 against T2 would violate independence.

The paired analysis includes only patients with accepted T1 and T2 observations: 6 children and 16 adults. Clonality and Shannon diversity use the complete repertoire. Overlap uses expanded clonotypes.

## Expansion and clone identity

Expansion is defined as at least 0.10% of productive TRB reads within a patient observation. This operational threshold excludes the low-frequency repertoire from overlap calculations; bulk genomic DNA TRB sequencing does not identify naive, memory or effector phenotype. Threshold sensitivity is evaluated separately.

The early overlap analyses use `aaV`: CDR3 beta amino acid sequence and TRBV. Final physical clone tracking uses CDR3 nucleotide sequence, TRBV and TRBJ. Functional and public comparisons in the final analysis use CDR3 amino acid sequence, TRBV and TRBJ.

## Immune amnesia hypothesis

The hypothesis is that clonotypes representing established adaptive memory at T1 become depleted after measles. Readouts include changes in whole-repertoire clonality and Shannon diversity, expanded-set Jaccard and the fractions of T1-expanded clones retained above threshold, detected below threshold or absent at T2.

New T2 expansions are analyzed separately. If clonality rises after recovery, these clones are exported for exact specificity matching against the available reference. An absent match does not establish a lack of antigen specificity; reference coverage is reported explicitly.
