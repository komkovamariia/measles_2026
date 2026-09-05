# Source notes and authoritative cohort rules

## Current metadata source

The authoritative biological metadata is `config/measles_metadata_v2_age18_final.xlsx`

A plain CSV copy used by the analysis scripts is stored as `config/analysis_metadata.csv`

The current final rules are

- minimum library quality is 4400 productive TRB clonotypes
- children are younger than 18 years
- adults are 18 years and older
- healthy controls have clinical point `CONTROL` only
- measles cases use T1, T2 and T3
- T0 is not a study point
- T4 is part of the intended clinical design at about six months but no T4 sequencing library is present in the current dataset

## Biological meaning of the time points

T1 is the early infection point after measles infection and before the mature adaptive response is expected to be fully developed

T2 is the point immediately after clinical recovery

T3 is the follow up point about one month after recovery

The sequencing suffix is authoritative for the current data

- `.1` means T1
- `.2` and `.2а` mean T2
- `.3` means T3

Healthy controls are not assigned T1, T2 or T3 even when a raw sample label contains a numeric suffix

## Resolved biological repeat for БА1_071

`new_106` is `106_БА1_071.2`

`new_111` is `111_БА1_071.2а`

Both are accepted T2 libraries for the same patient and both pass the final 4400 clonotype quality threshold

They are not counted as two independent T2 observations

For patient level repertoire analyses each library is first normalized to clonotype frequencies that sum to one. The normalized frequencies are then averaged across the two T2 libraries. Missing clonotypes in one replicate contribute zero in that replicate. The resulting averaged T2 repertoire contributes one biological observation for patient БА1_071

This prevents the deeper library from receiving greater weight

## Legacy cohort reconciliation

UDP labels are sequencing adapter identifiers and are not biological patient identifiers

The corrected mappings are

- UDP0097 is МР1_001.1, patient МР1_001, T1, measles
- UDP0098 is БА1_002.1, patient БА1_002, T1, measles
- UDP0099 is ГА1_003.1, patient ГА1_003, T1, measles
- UDP0100 is ХИ1_004.1, patient ХИ1_004, T1, measles
- UDP0101 is ИМ1_005.1, patient ИМ1_005, T1, measles
- UDP0102 is ПИ1_006.1, patient ПИ1_006, T1, measles
- UDP0103 is СС1_007.1, patient СС1_007, T1, measles
- UDP0104 is ХИ1_004.2, patient ХИ1_004, T2, measles
- UDP0105 is ШЕ1_008.1, patient ШЕ1_008, T1, measles
- UDP0106 through UDP0112 are healthy adult controls CGIMU_N1 through CGIMU_N7

UDP0108 remains excluded from biological inference because its repertoire is unusually similar to UDP0107

The T1 sample ГА1_003.1 is from the legacy sequencing batch while ГА1_003.2 is from the new batch. This paired patient is retained but should be checked in a sensitivity analysis because the two clinical points come from different sequencing batches

## Main analysis questions

### Cross sectional analysis

All accepted T1 patient observations are compared with age matched healthy controls. T1 observations remain in the analysis even when the same patient has no T2 sample

All accepted T2 patient observations are also summarized against healthy controls

T1 and T2 are not treated as independent groups for a direct independent test because some observations belong to the same patient

### Paired longitudinal analysis

The primary immune amnesia analysis uses only patients with accepted T1 and T2 observations

The current paired sample sizes are

- children younger than 18 years: 6 patients
- adults 18 years and older: 16 patients

Clonality and Shannon diversity are calculated on the complete patient repertoire

Overlap is calculated only among expanded clonotypes

## Expanded clonotype definition

The working expanded clonotype threshold is 0.10 percent of productive TRB reads within a biological patient observation

The overlap feature is `aaV`, defined as CDR3 beta amino acid sequence together with TRBV gene

This threshold removes the low frequency repertoire from the overlap analysis. Bulk genomic DNA TRB sequencing does not directly identify naive, memory or effector phenotype, so the analysis does not claim that every excluded low frequency clonotype is biologically naive

The 0.10 percent threshold is consistent with published TCR repertoire studies that use 0.10 percent as a highly expanded clone threshold and with replicate work showing that TCR beta frequencies around this level are substantially more stable across replicate samples

References

https://pmc.ncbi.nlm.nih.gov/articles/PMC4998859/

https://pmc.ncbi.nlm.nih.gov/articles/PMC10374308/

## Immune amnesia interpretation

The main hypothesis is that clonotypes representing established adaptive immune memory at T1 are depleted after measles

The primary readouts are

- change in complete repertoire clonality from T1 to T2
- change in Shannon diversity from T1 to T2
- fraction of expanded T1 clonotypes that remain expanded at T2
- fraction detected at T2 but below the expansion threshold
- fraction no longer detected at T2
- expanded clone Jaccard similarity between T1 and T2

New T2 expansions are analyzed separately. If clonality rises after recovery, new expanded clonotypes are exported for exact specificity matching against the available reference database

An absent database match is not evidence that a clonotype has no antigen specificity. Reference coverage is reported explicitly
