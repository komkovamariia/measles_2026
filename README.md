# Measles TCR-beta repertoire analysis

A deliberately small Jupyter pipeline for the human measles TCR-beta project.

The workflow answers five practical questions:

1. Which TRB aaV clonotypes are retained, lost or newly detected between time points of the same patient?
2. Which TRBV genes disappear, appear, increase or decrease by their share of productive TRB reads?
3. Are there recurrent TRBV patterns across adults/children and healthy/acute measles groups?
4. Are measles- and other infection-associated CDR3β clonotypes differently represented during measles?
5. Does overall TCRβ repertoire diversity change during disease and recovery?

The implementation stays intentionally simple: MiXCR preprocessing, two notebooks, Venn overlaps, TRBV read fractions and simple non-parametric statistics.

## Patient identity and time points

For the new cohort, the first integer in a study-table record is the sequencing-library index. The biological patient ID is the complete code after that index and before the final time-point suffix.

```text
15_CA1_013.1
19_CA1_013.2
```

are the same patient `CA1_013` at time points 1 and 2. `22_CA1_017.1` is a different patient `CA1_017`.

The final `.1`, `.2`, `.3` means admission, clinical recovery/temperature normalization and later follow-up respectively.

FASTQ files are linked to the study table by the integer before the first underscore. For example all files from `10_S1_L001_R1_001.fastq.gz` through lane L004 and read mate R2 belong to sequencing index `10`.

## Input cohorts

New FASTQ directory:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

A complete new library must contain all four paired lanes L001-L004. The preprocessing step writes `results/fastq_audit.csv` and stops a full run if an included sample is incomplete.

Legacy FASTQ directory:

```text
/projects/mice_transplant_2025/measles/raw/pe150_res
```

The old files are paired as `_1.fq.gz` and `_2.fq.gz`, for example:

```text
FT150008642_56_L01_UDP0112_1.fq.gz
FT150008642_56_L01_UDP0112_2.fq.gz
```

The study owner clarified the old-cohort interpretation:

1. UDP0097-UDP0105 are adult measles cases and are assigned to analysis time point 1;
2. UDP0106-UDP0112 are healthy adult controls;
3. every legacy UDP library represents a different adult person;
4. legacy UDP identities are independent from all new-cohort patient codes and are never joined longitudinally to them;
5. UDP0108 is processed for the historical UDP0107/UDP0108 QC comparison but remains excluded from biological inference.

The complete analysis requires both cohorts. `run_all_from_scratch.sh` locates `pe150_res` automatically and stops if required input files are missing.

## One-shot run from scratch

```bash
git clone https://github.com/komkovamariia/measles_2026.git
cd measles_2026
nohup bash run_all_from_scratch.sh > logs/bootstrap_submit.log 2>&1 &
```

The bootstrap script creates or updates the pinned conda environment, validates the repository, caches the IEDB TCRMatch reference, locates both sequencing cohorts and submits the full Slurm job. No pilot is required for this command.

## Patient-level outputs

For every new-cohort patient with at least two time points:

```text
figures/patients/<PATIENT>_aav_venn.png/pdf
figures/patients/<PATIENT>_T1_to_T2_trbv_change.png/pdf
figures/patients/<PATIENT>_T2_to_T3_trbv_change.png/pdf

results/longitudinal/patients/<PATIENT>_T1_to_T2_trbv_changes.csv
results/longitudinal/patients/<PATIENT>_T2_to_T3_trbv_changes.csv
```

The Venn feature is the same simple `aaV = (CDR3aa, TRBV)` definition used by the mouse project. Separate CDR3-only and aaV overlap tables are also saved.

For every TRBV transition the table contains:

```text
fraction_old
fraction_new
delta_fraction
delta_percentage_points
presence_change = gained / lost / persistent / absent_both
direction = increased / decreased / unchanged
```

TRBV fraction is:

```text
reads assigned to TRBV / all productive TRB reads in the biological sample
```

## Group-level TRBV analysis

The project calculates:

1. adult acute versus adult healthy;
2. child acute versus child healthy;
3. adults versus children within acute measles;
4. adults versus children within healthy controls;
5. paired T1→T2 and T2→T3 changes separately in adults and children;
6. adult versus child differences in the patient-level longitudinal TRBV change.

Mann-Whitney U is used for independent groups and Wilcoxon signed-rank for paired longitudinal measurements. A p-value is omitted when a comparison has fewer than three independent observations per required group. Benjamini-Hochberg correction is applied within each TRBV comparison block.

Main outputs:

```text
results/longitudinal/trbv_group_statistics.csv
results/longitudinal/trbv_top_patterns.csv
figures/09_group_trbv_patterns.png/pdf
figures/10_group_trbv_gain_loss.png/pdf
```

## Other retained analyses

The project also keeps depth-normalized repertoire richness and Shannon diversity, top-10 clonotype dominance, TRBV/TRBJ usage, exact CDR3β matches to IEDB for measles, hepatitis B, poliovirus and tetanus, and longitudinal diversity trajectories.

Read downsampling is performed without replacement with seed `1031`.

## Project layout

```text
config/sample_manifest.csv                 source metadata
config/SOURCE_NOTES.md                     explicit cohort rules
notebooks/00_prepare_and_mixcr.ipynb       small preprocessing notebook
notebooks/01_measles_analysis.ipynb        patient and group analyses
scripts/preprocess_human_trb.py            FASTQ validation + MiXCR
scripts/analysis_01_prepare.py              biological sample assembly + diversity
scripts/analysis_05_longitudinal_venn.py   Venn + TRBV changes
scripts/run_pipeline.py                    notebook runner + output verification
scripts/validate_repository.py             static preflight checks
scripts/run_slurm.sh                       one 24-CPU Slurm job
run_all_from_scratch.sh                    environment + validation + submit
results/                                   CSV tables
figures/                                   PNG and PDF figures
audit_runs/                                executed notebooks
logs/                                      notebook and MiXCR logs
```

## Interpretation limits

The libraries are bulk TRB data from genomic-DNA multiplex sequencing. Read counts provide relative abundance and may contain multiplex-PCR bias, so TRBV changes are reported as read fractions and diversity is additionally assessed after common-depth downsampling.

The default MiXCR preset is `milab-human-dna-tcr-multiplex`. Keep it only if this is the confirmed wet-lab protocol; otherwise set `MEASLES_MIXCR_PRESET` to the appropriate protocol preset before submission.
