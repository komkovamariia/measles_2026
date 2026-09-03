# Source rules used by the pipeline

## New cohort

`Корь_Влада.docx` is the source for the 96 new sequencing records. The integer before the first underscore maps the FASTQ prefix to the study-table record. The full biological code such as `СА1_013` is the patient ID. The final `.1/.2/.3` is the clinical time point.

Orange cells in the child/adult lists define healthy controls when that patient is present in the new sequencing table. The age lists are explicitly incomplete, therefore patients absent from both lists remain `unknown` and are excluded from age-stratified inference.

## Legacy cohort

The sequencing table maps UDP0097-UDP0112 to the old libraries. The explicit colleague rule is used for disease status:

- UDP0106-UDP0112: healthy CGIMU controls;
- UDP0097-UDP0105: measles cases assigned to analysis time point 1;
- UDP0108: processed for QC but excluded from biological inference because of the previously identified UDP0107 similarity/depth problem.

## Documented age-source inconsistency

The colleague prose says that children occur only in the new `fastq_files.tar` batch. The patient list in `Корь_Влада.docx` simultaneously labels several patient IDs that are present in the legacy table as children, including `ГА1_003`, `ХИ1_004`, `СС1_007` and `ШЕ1_008`. The new cohort also contains follow-up records such as `46_ГА1_003.2`, which links directly to the legacy `ГА1_003.1` patient.

The current manifest therefore keeps age attached to the biological patient ID from the child/adult table, so one person cannot change age group between old and new sequencing batches. This choice should be revisited only if the study owner provides a corrected age table or explicitly states that the legacy and new identical patient codes refer to different people.
