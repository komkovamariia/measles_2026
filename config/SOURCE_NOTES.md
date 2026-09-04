# Source rules used by the pipeline

## New cohort

`Корь_Влада.docx` is the source for the 96 sequencing-library records and their sequencing indices. The integer before the first underscore maps a FASTQ prefix to a sequencing record. The raw sequencing label is preserved in `raw_sample_id` for traceability.

The uploaded clinical workbook `Корь метадата.xlsx` is used for disease status and age metadata for sequenced new-cohort patients. The study owner subsequently clarified the authoritative longitudinal rule from the sequencing/sample list: **there is no T0 in this study** and the final suffix denotes the study time point directly:

- `.1` -> `T1`;
- `.2` and labelled variants such as `.2а` -> `T2`;
- `.3` -> `T3`.

The analysis-ready fields are stored in `config/sample_manifest.csv`:

- `disease_status`: `measles` / `healthy`;
- `is_control`: healthy-control flag;
- `age_group`: `adult` / `child`, with the current study rule child < 6 years and adult >= 6 years;
- `clinical_timepoint`: `T1`, `T2`, or `T3` for new sequencing libraries;
- `include`: whether the library is currently eligible for biological inference.

Patients present in the clinical workbook but absent from sequencing are not added to the repertoire analysis.

## Authoritative time-point rule: no T0

The previous interpretation that clinical workbook columns could remap the numeric suffix to another phase has been retired. For the new sequencing cohort, `analysis_timepoint` and `clinical_timepoint` must agree with the suffix in `raw_sample_id`.

Examples:

- `СА1_013.1/.2/.3` -> `T1/T2/T3`;
- `ММ1_033.1/.2/.3` -> `T1/T2/T3`;
- `ЭГ1_042.1/.2/.3` -> `T1/T2/T3`;
- `ПМ1_009.1/.2` -> `T1/T2`;
- `ГА1_003.2` -> `T2`;
- `ОЛ1_018.2` -> `T2`;
- `БА1_071.2а` -> `T2` for the time-point field.

`T0` is invalid and repository validation fails if any new sequencing library is assigned to `T0`.

## Reconciled patient identifiers

The following sequencing labels remain reconciled to patient identifiers in the clinical workbook while preserving the raw sequencing label. The longitudinal point itself is determined from the raw suffix:

- `21_ОЛ1_016.1` -> patient `ОЛ1_018`, `T1`; `75_ОЛ1_018.2` -> the same patient, `T2`;
- `22_СА1_017.1` -> patient `СА1_019`, `T1`;
- `23_МС1_018.1` -> patient `МС1_020`, `T1`;
- `57_НЛ1_043.2` -> patient `НЛ1_043`, `T2`; the workbook cell is written as `НЛ1_0432` without the dot.

## QC and unresolved sequencing/identity conflicts

The following libraries retain `include=False` and are excluded from biological group/longitudinal inference until resolved. Their time-point suffix is still recorded independently where possible:

- `new_24` / `24_ФО1_019.1`: patient identity is unresolved relative to the clinical workbook; time point is `T1` from the suffix;
- `new_39` / `39_ММ1_033.2`: strong sample-identity QC warning. Abundance-weighted TRB similarity is approximately 0.00004 for `.1` vs `.2` and 0.00016 for `.2` vs `.3`, while `.1` vs `.3` is approximately 0.57. This pattern is much more consistent with a sample identity/mislabelling problem than ordinary longitudinal turnover. The time points are nevertheless `T1`, `T2`, `T3` by suffix;
- `new_88` / `88_КА1_059.1`: a second sequencing library is again labelled `КА1_059.1`; it is `T1` by suffix but remains excluded pending clarification of whether it is a technical replicate or mislabelled sample;
- `new_111` / `111_БА1_071.2а`: time point is `T2` by suffix, but the library remains excluded by the current curated config pending clarification of the `2а` sample identity.

## СА1_013 T1 MiXCR export

`new_15` / `15_CA1_013.1` is unequivocally the `T1` sample of patient `СА1_013`. Its current MiXCR TSV was header-only and contained zero productive TRB clonotypes. It had been incorrectly accepted as a valid resume output because the previous validity check only inspected the header. The pipeline now removes invalid/header-only exports before resume and requires at least one positive productive TRB clonotype. Until that sample is successfully regenerated, T1 is shown in patient QC but excluded from biological overlap; valid `T2` and `T3` can still be compared.

## Antigen-specificity reference limitation

The current `IEDB_data.tsv` TCRMatch subset used by the project contains 763 unique Hepatitis B reference CDR3 sequences but no reference records labelled Measles, Poliovirus or Tetanus. Therefore zero exact matches for those three pathogens are **not biological evidence of absence**. The pipeline records reference coverage explicitly and reports unavailable pathogens as missing/unsupported rather than zero. Exact matching should be treated as a conservative annotation layer; a validated broader specificity/similarity reference is required before the measles and immune-amnesia hypotheses can be answered from antigen annotation.

## Legacy cohort

Old FASTQ files are stored under `/projects/mice_transplant_2025/measles/raw/pe150_res` and are named like `FT150008642_56_L01_UDP0112_1.fq.gz` and `_2.fq.gz`.

All legacy libraries are adults and are different people from the new cohort. Legacy identifiers remain independent `UDP0097` ... `UDP0112` patients and are never linked longitudinally to new biological codes.

- UDP0097-UDP0105: adult measles cases, retained as `legacy_acute` for sensitivity analyses;
- UDP0106-UDP0112: healthy adult controls;
- UDP0108: processed for QC but excluded from biological inference.

Primary case-control summaries use the new sequencing batch with one biological sample per patient to limit batch confounding. Adult all-batch comparisons are retained as sensitivity analyses.

## Downstream grouping contract

Biological analysis uses:

- `patient_id` for biological identity;
- `clinical_timepoint` for longitudinal phase, with new-cohort `.1/.2/.3 = T1/T2/T3` and no T0;
- `is_control` / `disease_status` for control-vs-measles comparisons;
- `age_group` for adult/child stratification using child < 6 years, adult >= 6 years;
- `include=False` to exclude unresolved/QC-failed libraries from biological inference while retaining traceability.
