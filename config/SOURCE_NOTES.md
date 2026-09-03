# Source rules used by the pipeline

## New cohort

`Корь_Влада.docx` is the source for the 96 new sequencing records. The integer before the first underscore maps the FASTQ prefix to the study-table record. The full biological code such as `СА1_013` is the patient ID. The final `.1/.2/.3` is the clinical time point.

Orange cells in the child/adult lists define healthy controls when that patient is present in the new sequencing table. The age lists are incomplete, therefore patients absent from both lists remain `unknown` and are excluded from age-stratified inference.

## Legacy cohort

The old FASTQ files are stored under `/projects/mice_transplant_2025/measles/raw/pe150_res` and are named like `FT150008642_56_L01_UDP0112_1.fq.gz` and `FT150008642_56_L01_UDP0112_2.fq.gz`.

The study owner clarified that all legacy libraries are adults and are different people from the patients in the new cohort. Therefore the legacy patient identifiers are kept as independent `UDP0097` ... `UDP0112` IDs and are never linked longitudinally to similarly numbered/new biological codes.

Disease-status rules:

- UDP0097-UDP0105: adult measles cases assigned to analysis time point 1;
- UDP0106-UDP0112: healthy adult controls;
- UDP0108: processed for the historical UDP0107/UDP0108 QC comparison but excluded from biological inference because of the previously identified similarity/depth problem.

The earlier tentative mapping of legacy UDP libraries to new patient IDs is retired and must not be restored.
