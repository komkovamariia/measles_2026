# Run guide

## 1. Clone the project

```bash
cd ~
git clone https://github.com/komkovamariia/measles_2026.git
cd measles_2026
```

## 2. Create the environment once

```bash
conda env create -f environment.yml
conda activate measles-tcr-2026
```

Check that MiXCR is available:

```bash
mixcr --version
```

## 3. Data paths

The new batch is already configured to use:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

Check it:

```bash
find /projects/mice_transplant_2025/raw/20260817_161454/Fastq -maxdepth 1 -name '*.fastq.gz' | head
```

For the old UDP0097-UDP0112 samples, first try:

```bash
find /projects/mice_transplant_2025 -type f -name '*UDP0107*R1*.fastq.gz' -print -quit
```

If a file is returned, set its directory once, for example:

```bash
export MEASLES_LEGACY_FASTQ_DIR=/path/to/old/pe150_res
```

The Slurm script performs the same UDP0107 search automatically. If the old data are not found, the new cohort still runs and legacy samples are marked as missing FASTQ in the metadata QC table.

## 4. Pilot

Run a small test before spending cluster time:

```bash
python scripts/run_pipeline.py \
  --stage all \
  --pilot \
  --cpus 8 \
  --threads-per-sample 2
```

The pilot attempts representative new libraries and the UDP0107/UDP0108 QC pair when legacy files are available.

Check:

```bash
ls results/clonotypes | head
ls figures
cat results/result_summary.txt
```

## 5. Full Aldan-3 run

```bash
sbatch scripts/run_slurm.sh
```

The job requests 24 CPUs. MiXCR uses 4 threads per sample, therefore at most 6 independent libraries are processed simultaneously.

Monitor:

```bash
squeue -u "$USER"
tail -f logs/slurm_JOBID.out
```

After completion:

```bash
sacct -j JOB_ID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
```

`64G` is an initial memory reservation. It can be reduced later after checking `MaxRSS`.

## 6. If the job stops

Run the same submission again:

```bash
sbatch scripts/run_slurm.sh
```

Completed MiXCR exports are detected and skipped, so the analysis continues from the remaining libraries.

## 7. Main outputs

```text
results/sample_manifest_resolved.csv
results/metadata_issues.csv
results/diversity_metrics.csv
results/antigen_specificity_by_sample.csv
results/antigen_specific_clonotypes.csv
results/longitudinal_metrics.csv
results/statistical_tests.csv
results/result_summary.txt
```

Main figures:

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
