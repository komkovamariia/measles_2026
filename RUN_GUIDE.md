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

`run_all_from_scratch.sh` updates the `measles-tcr-2026` environment, validates metadata and code, checks FASTQ, verifies MiXCR, and submits the Slurm job.

## Analysis-only rerun: use this after MiXCR has already completed

This is the recommended command after metadata/statistics/plot changes. It does **not** rerun MiXCR.

```bash
cd ~/measles_2026
git pull --ff-only

eval "$(conda shell.bash hook)"
conda activate measles-tcr-2026

python scripts/validate_repository.py

# Raw FASTQ depth is collected only if the table does not already exist.
if [ ! -s results/raw_fastq_depth.csv ]; then
  python scripts/collect_fastq_depth.py --workers 8
fi

python scripts/run_pipeline.py \
  --stage analysis \
  --resume \
  --cpus 24 \
  --threads-per-sample 4
```

The analysis notebook now runs:

```text
analysis_01_prepare.py
analysis_02_antigens.py
analysis_06_clinical_summary.py
analysis_07_age6_reporting.py
```

The project-wide age rule is:

```text
child = age < 6 years
adult = age >= 6 years
```

## Where to look after the analysis-only rerun

Main corrected statistics and QC:

```text
results/corrected/
figures/corrected/
figures/corrected/patients/
```

Additional age-6 reporting requested for interpretation:

```text
results/corrected_age6/
figures/corrected_age6/
```

Especially useful files:

```text
results/corrected_age6/age6_group_summary.csv
results/corrected_age6/timepoint_age_audit.csv
results/corrected_age6/TRBV_priority_scores.csv
results/corrected_age6/venn_T1_T2_vs_healthy_summary.csv

figures/corrected_age6/TRBV_usage_adult_annotated.png
figures/corrected_age6/TRBV_usage_child_annotated.png
figures/corrected_age6/TRBJ_usage_adult_annotated.png
figures/corrected_age6/TRBJ_usage_child_annotated.png
figures/corrected_age6/TRBV_priority_scores_adult.png
figures/corrected_age6/TRBV_priority_scores_child.png
figures/corrected_age6/venn_summary_T1_T2_vs_healthy.png
```

In annotated TRBV/TRBJ heatmaps each cell contains:

```text
patients with this gene / all patients in that group
total reads assigned to this gene across the group
```

The color represents mean percentage of productive TRB reads assigned to that gene.

## Data

New batch:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

Each new sequencing index has four paired lanes L001-L004. The first integer in the FASTQ name maps to `config/sample_manifest.csv`.

Legacy batch:

```text
/projects/mice_transplant_2025/measles/raw/pe150_res
```

Old files use names such as:

```text
FT150008642_56_L01_UDP0112_1.fq.gz
FT150008642_56_L01_UDP0112_2.fq.gz
```

Legacy biological rules remain fixed:

```text
UDP0097-UDP0105 = adult measles cases
UDP0106-UDP0112 = healthy adult controls
```

All 16 legacy UDP libraries are independent adult people. UDP0108 is processed for QC but excluded from biological inference.

## Clinical metadata and group assignment

`config/clinical_metadata_sequenced.csv` and `config/sample_manifest.csv` contain the curated sequenced cohort. `clinical_timepoint` is explicit and is not inferred from the raw `.1/.2/.3` suffix.

Important examples:

```text
СА1_013: .1=T1, .2=T2, .3=T3
ММ1_033: .1=T1, .2=T2, .3=T2
ЭГ1_042: .1=T0, .2=T1, .3=T2
```

`new_15 / СА1_013.1` currently has zero productive TRB after MiXCR, therefore it fails biological QC even though the clinical T1 sample exists.

`new_39 / ММ1_033.2` is excluded from biological inference because its repertoire identity is a strong outlier relative to `.1` and `.3`.

Current unresolved/excluded new libraries are:

```text
new_24
new_39
new_88
new_111
```

## Slurm full run

The full Slurm job requests 24 CPUs and 64 GB RAM. MiXCR receives 4 threads per library, so at most 6 libraries are processed simultaneously.

Monitor:

```bash
squeue -u "$USER"
tail -F logs/slurm_JOBID.out
```

After completion:

```bash
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
```

## Validation

Before analysis:

```bash
python scripts/validate_repository.py
```

The validator now checks the age rule explicitly: the included child set must correspond to patients younger than 6 years, while patients aged 6 years or older are adults.

After a full run `scripts/run_pipeline.py` verifies both the corrected outputs and the additional `corrected_age6` summary tables/figures.

## MiXCR license

MiXCR requires a valid license on the server. If needed:

```bash
conda activate measles-tcr-2026
mixcr activate-license
```
