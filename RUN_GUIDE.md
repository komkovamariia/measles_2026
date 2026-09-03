# Run guide

## One unattended full run

On Aldan-3:

```bash
cd ~
if [ -d "$HOME/measles_2026/.git" ]; then
  cd "$HOME/measles_2026"
  git pull --ff-only
else
  git clone git@github.com:komkovamariia/measles_2026.git "$HOME/measles_2026"
  cd "$HOME/measles_2026"
fi

nohup bash run_all_from_scratch.sh > logs/bootstrap_submit.log 2>&1 &
echo "bootstrap PID: $!"
```

`run_all_from_scratch.sh` performs the steps in sequence:

1. creates or updates `measles-tcr-2026` from `environment.yml`;
2. checks MiXCR and runs `scripts/validate_repository.py`;
3. verifies the new FASTQ directory;
4. locates the old `pe150_res` cohort;
5. caches the IEDB TCRMatch reference before the compute job;
6. submits `scripts/run_slurm.sh`.

The full run stops if an included old or new library is missing. This prevents accidental analysis of a partial cohort.

## Data

New batch:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

Every new sequencing index must have exactly four paired lanes L001-L004. The first integer in the FASTQ name maps to `config/sample_manifest.csv`.

Legacy batch:

```text
/projects/mice_transplant_2025/measles/raw/pe150_res
```

Old files use names such as:

```text
FT150008642_56_L01_UDP0112_1.fq.gz
FT150008642_56_L01_UDP0112_2.fq.gz
```

The preprocessor supports this `_1.fq.gz` / `_2.fq.gz` pairing directly.

Legacy biological rules are fixed:

```text
UDP0097-UDP0105 = adult measles cases, analysis time point 1
UDP0106-UDP0112 = healthy adult controls
```

All 16 legacy UDP libraries are different adult people and are independent from all patients in the new cohort. They must never be joined longitudinally to new patient codes. UDP0108 is processed for QC but excluded from biological inference.

## Patient identity

For the new cohort, the complete biological code is the patient ID:

```text
15_СА1_013.1  -> patient СА1_013, time point 1
19_СА1_013.2  -> patient СА1_013, time point 2
22_СА1_017.1  -> patient СА1_017, time point 1
```

Therefore `СА1_013` and `СА1_017` are different people.

Legacy patient identity is simply the independent UDP code itself, for example `UDP0097` or `UDP0112`.

## Slurm

The submitted job requests 24 CPUs and 64 GB RAM. MiXCR receives 4 threads per library, so at most 6 independent libraries are processed simultaneously.

Monitor:

```bash
cat logs/bootstrap_submit.log
squeue -u "$USER"
tail -f logs/slurm_JOBID.out
```

After completion:

```bash
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
```

## Automatic validation

Before analysis, `scripts/validate_repository.py` checks the required files, notebook structure, Python syntax, 96 new records, 16 legacy records, patient identity, timepoint parsing, control mapping, corrected legacy adult case/control rules and pinned environment.

During preprocessing, `scripts/preprocess_human_trb.py` checks actual FASTQ availability. New samples require all eight FASTQ files corresponding to four paired lanes. Legacy samples require a valid paired old `_1/_2` or R1/R2 file pair. The audit is saved to:

```text
results/fastq_audit.csv
results/metadata_issues.csv
results/sample_manifest_resolved.csv
```

After the complete pipeline, `scripts/run_pipeline.py` verifies the expected result tables, figures 01-10 and at least one patient aaV Venn figure.

## Patient results

For each new-cohort patient with at least two time points:

```text
figures/patients/<PATIENT>_aav_venn.png/pdf
figures/patients/<PATIENT>_T1_to_T2_trbv_change.png/pdf
results/longitudinal/patients/<PATIENT>_T1_to_T2_trbv_changes.csv
```

When T3 exists, a T2→T3 TRBV table and figure are also created.

The Venn feature is `aaV = (CDR3aa, TRBV)`, matching the mouse-project Venn definition. Numerical CDR3-only and aaV overlaps are saved in:

```text
results/longitudinal/patient_cdr3_overlap.csv
results/longitudinal/patient_aav_overlap.csv
```

TRBV change is based on:

```text
TRBV read fraction = reads assigned to a TRBV / all productive TRB reads in that sample
```

## Group results

```text
results/longitudinal/trbv_group_statistics.csv
results/longitudinal/trbv_top_patterns.csv
figures/09_group_trbv_patterns.png/pdf
figures/10_group_trbv_gain_loss.png/pdf
```

The analysis includes adult/child case-control comparisons, within-status adult/child comparisons, paired T1→T2 and T2→T3 changes, and age-group differences in longitudinal change where sample sizes permit statistical testing.

Other retained results include diversity, top-clone dominance, V/J usage and exact IEDB matches for measles, hepatitis B, poliovirus and tetanus.

## Important external requirement

MiXCR requires a valid MiXCR license on the server. If the account has not previously activated MiXCR, activate the academic license once before leaving an unattended full run.
