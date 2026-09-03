# Run guide

## 1. Clone the project

```bash
cd ~
git clone https://github.com/komkovamariia/measles_2026.git
cd measles_2026
```

If the repository is already cloned:

```bash
cd ~/measles_2026
git pull
```

## 2. Create the environment once

```bash
conda env create -f environment.yml
conda activate measles-tcr-2026
```

If the environment already existed before `matplotlib-venn` was added:

```bash
conda env update -f environment.yml --prune
conda activate measles-tcr-2026
```

Check MiXCR:

```bash
mixcr --version
```

## 3. Data path

The current new FASTQ directory is already configured:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

Check several files:

```bash
find /projects/mice_transplant_2025/raw/20260817_161454/Fastq \
  -maxdepth 1 -name '*.fastq.gz' | head
```

The code identifies a library by the integer before the first underscore. All four lanes with this prefix belong to the same sequencing library.

The study-table naming is interpreted as follows:

```text
10_PM1_009.1
12_PM1_009.2
```

`10` and `12` are sequencing-library indices, `PM1` is one patient, and `.1/.2` are different clinical time points. The analysis therefore joins these records into one longitudinal patient.

For optional legacy UDP0097-UDP0112 data:

```bash
export MEASLES_LEGACY_FASTQ_DIR=/path/to/old/pe150_res
```

## 4. Pilot

```bash
python scripts/run_pipeline.py \
  --stage all \
  --pilot \
  --cpus 8 \
  --threads-per-sample 2
```

Check:

```bash
ls results/clonotypes | head
ls results/longitudinal
ls figures/patients | head
```

The pilot is mainly a technical check. Full patient Venn diagrams require the relevant repeated time points to have been processed.

## 5. Full Aldan-3 run

```bash
sbatch scripts/run_slurm.sh
```

The job requests 24 CPUs. MiXCR uses 4 threads per sample, so up to 6 libraries are processed at the same time.

Monitor:

```bash
squeue -u "$USER"
tail -f logs/slurm_JOBID.out
```

After completion:

```bash
sacct -j JOB_ID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
```

## 6. Restart after interruption

Run the same command again:

```bash
sbatch scripts/run_slurm.sh
```

Existing non-empty MiXCR TRB exports are skipped automatically.

## 7. Patient-level results

Each patient with at least two time points receives a Venn diagram:

```text
figures/patients/<PATIENT>_cdr3_venn.png
figures/patients/<PATIENT>_cdr3_venn.pdf
```

For two time points the plot is a 2-set Venn. For three time points it is a 3-set Venn. Sets contain exact CDR3β amino-acid clonotypes.

For each new time point the pipeline also creates a TRBV read-share change table and plot:

```text
results/longitudinal/patients/<PATIENT>_T1_to_T2_trbv_changes.csv
figures/patients/<PATIENT>_T1_to_T2_trbv_change.png/pdf
```

When a third time point exists:

```text
results/longitudinal/patients/<PATIENT>_T2_to_T3_trbv_changes.csv
figures/patients/<PATIENT>_T2_to_T3_trbv_change.png/pdf
```

Important columns:

```text
trbv
fraction_old
fraction_new
delta_percentage_points
presence_change
direction
```

`gained` means the TRBV had zero reads at the previous time point and non-zero reads at the new point. `lost` means the reverse. `persistent` means it is detected at both points.

## 8. Statistical group results

```text
results/longitudinal/trbv_group_statistics.csv
```

The statistical unit is a biological patient.

For adults and children separately:

```text
acute measles vs healthy control : Mann-Whitney U
T1 vs T2 within the same patient : paired Wilcoxon
T2 vs T3 within the same patient : paired Wilcoxon
```

The tested value is the fraction of all productive TRB reads assigned to each TRBV. Benjamini-Hochberg FDR correction is applied within each comparison.

Group figures:

```text
figures/09_group_trbv_patterns.png/pdf
figures/10_group_trbv_gain_loss.png/pdf
```

## 9. Other retained results

```text
results/diversity_metrics.csv
results/antigen_specificity_by_sample.csv
results/antigen_specific_clonotypes.csv
results/statistical_tests.csv
```

Existing study figures are retained:

```text
01_task1_measles_group_comparison.png/pdf
02_task1_oligoclonal_dominance.png/pdf
03_task2_other_infections_group_comparison.png/pdf
04_task2_other_infections_longitudinal.png/pdf
05_task3_diversity.png/pdf
06_task3_vj_usage.png/pdf
07_task4_richness_longitudinal.png/pdf
08_task4_shannon_longitudinal.png/pdf
```
