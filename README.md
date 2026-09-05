# Measles TCR beta repertoire analysis

This repository contains the final analysis workflow for bulk human TRB sequencing in the measles cohort

The current biological metadata is stored in `config/measles_metadata_v2_age18_final.xlsx`

The analysis copy used by scripts is `config/analysis_metadata.csv`

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

## Two complementary analysis designs

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

## Expanded clonotypes

An expanded clonotype is defined as an `aaV` clone with at least 0.10 percent of productive TRB reads in a biological patient observation

`aaV` means CDR3 beta amino acid sequence together with TRBV gene

The expanded clone filter removes the low frequency background from overlap analysis. Bulk genomic DNA TRB sequencing does not directly phenotype naive, memory or effector T cells, therefore the workflow does not label every low frequency clone as naive

The primary immune amnesia readouts are

- T1 to T2 change in clonality
- T1 to T2 change in Shannon diversity
- expanded T1 clones that remain expanded at T2
- expanded T1 clones detected at T2 below the expansion threshold
- expanded T1 clones not detected at T2
- expanded clone Jaccard similarity

New expanded clones at T2 are exported separately for exact specificity matching against the available reference database

## БА1_071 T2 biological repeat

`new_106` and `new_111` are both accepted T2 libraries for patient БА1_071

`new_111` is the raw sample `111_БА1_071.2а`

The two libraries are not counted as two patients or two independent T2 observations

For patient level analysis each library is first normalized to clonotype frequencies. The normalized frequencies are then averaged. The resulting averaged repertoire is one T2 observation for patient БА1_071

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

Important tables include

```text
cross_sectional_cohort_counts.csv
cross_sectional_statistics.csv
longitudinal_expanded_clone_summary.csv
expanded_clone_fates.csv
new_T2_expansions.csv
expanded_TRBV_usage.csv
patient_panel_index.csv
```

Patient panels show complete repertoire diversity, overlap of expanded clonotypes only, fate of expanded clones from the earlier point and the largest new expansions at the later point

## Refresh analysis without rerunning MiXCR

```bash
cd ~/measles_2026
git pull --ff-only
python scripts/validate_repository.py
jid=$(sbatch --parsable scripts/run_analysis_refresh_slurm.sh)
echo "$jid"
tail -F "logs/refresh_slurm_${jid}.out"
```

The refresh job reuses existing MiXCR clonotype exports

## Full run from raw FASTQ

```bash
bash run_all_from_scratch.sh
```

The full workflow uses MiXCR 4.7.0 with the project environment and the configured human DNA TCR multiplex preset

## Interpretation limits

Read fractions from multiplex genomic DNA TRB sequencing are relative abundance measures and can contain PCR amplification bias

The expanded clone threshold is an operational repertoire threshold, not a direct cell phenotype measurement

Exact database matching has incomplete antigen coverage. A missing match must not be interpreted as absence of biological antigen specificity

See `config/SOURCE_NOTES.md` for all cohort reconciliation decisions and analysis definitions
