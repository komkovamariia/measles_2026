# Source rules used by the pipeline

## New cohort

`Корь_Влада.docx` remains the source for the 96 sequencing-library records and their sequencing indices. The integer before the first underscore maps the FASTQ prefix to the sequencing-table record. The raw sequencing label is preserved in `raw_sample_id` for traceability.

The uploaded workbook `Корь метадата.xlsx` is now the authoritative source for **disease status (measles/healthy) and age group (adult/child)** for new-cohort patients that are represented in sequencing. The repository stores the sequenced subset of those patient-level annotations in `config/clinical_metadata_sequenced.csv`, and the analysis-ready values are propagated into `config/sample_manifest.csv` as `disease_status`, `is_control`, and `age_group`.

Patients present in the clinical workbook but absent from the sequencing manifest are not added to the repertoire analysis.

### Reconciled patient identifiers

Three sequencing labels were reconciled to patient identifiers in the uploaded clinical workbook while preserving the original raw label:

- `21_ОЛ1_016.1` -> patient `ОЛ1_018`, time point 1. The same patient has `75_ОЛ1_018.2`.
- `22_СА1_017.1` -> patient `СА1_019`, time point 1.
- `23_МС1_018.1` -> patient `МС1_020`, time point 1.

### Unresolved sequencing/clinical conflicts

The following libraries are still processed by MiXCR for QC/traceability but have `include=False`, so they are excluded from biological group statistics until their identity/time point is confirmed:

- `new_24` / `24_ФО1_019.1`: no exact `ФО1_019` patient exists in the uploaded clinical workbook; `ФО1_017.1` is already represented by another sequenced library.
- `new_88` / `88_КА1_059.1`: the clinical workbook contains `КА1_059.1` and `КА1_059.2`, while both `new_81` and `new_88` are labelled `.1` in the sequencing table.
- `new_111` / `111_БА1_071.2а`: the clinical workbook contains `БА1_071.2` and `БА1_071.3`, but `.2а` does not unambiguously identify which biological sample was sequenced.

The analysis therefore uses only unambiguous included libraries for age/status comparisons and longitudinal biological inference.

## Legacy cohort

The old FASTQ files are stored under `/projects/mice_transplant_2025/measles/raw/pe150_res` and are named like `FT150008642_56_L01_UDP0112_1.fq.gz` and `FT150008642_56_L01_UDP0112_2.fq.gz`.

The study owner clarified that all legacy libraries are adults and are different people from the patients in the new cohort. Therefore the legacy patient identifiers are kept as independent `UDP0097` ... `UDP0112` IDs and are never linked longitudinally to similarly numbered/new biological codes.

Disease-status rules:

- UDP0097-UDP0105: adult measles cases assigned to analysis time point 1;
- UDP0106-UDP0112: healthy adult controls;
- UDP0108: processed for the historical UDP0107/UDP0108 QC comparison but excluded from biological inference because of the previously identified similarity/depth problem.

The earlier tentative mapping of legacy UDP libraries to new patient IDs is retired and must not be restored.

## Grouping contract used by downstream analysis

`scripts/analysis_01_prepare.py` reads the resolved manifest created from `config/sample_manifest.csv`. Group separation therefore follows the current manifest fields:

- `is_control=True` / `disease_status=healthy` -> healthy controls;
- `disease_status=measles` -> measles cases;
- `age_group=adult` or `child` -> age-stratified analyses;
- `include=False` -> processed if FASTQ is available, but excluded from biological inference.
