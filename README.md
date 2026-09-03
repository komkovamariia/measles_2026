# Measles TCR-beta repertoire analysis

A deliberately small Jupyter pipeline for the human measles TCR-beta project.

The workflow has five practical questions:

1. Which exact CDR3β clonotypes are retained, lost or newly detected between time points of the same patient?
2. Which TRBV genes disappear, appear, increase or decrease by their share of productive TRB reads?
3. Are there recurrent TRBV changes across children/adults and healthy/acute measles groups?
4. Are measles- and other infection-associated CDR3β clonotypes differently represented during measles?
5. Does overall TCRβ repertoire diversity change during disease and recovery?

The implementation is intentionally simple: MiXCR preprocessing, two notebooks, exact Venn overlaps, TRBV read fractions and simple non-parametric statistics.

## Patient identity and time points

For the new FASTQ batch:

```text
10_S1_L001_R1_001.fastq.gz ... 10_S1_L004_R2_001.fastq.gz
```

The integer before the first underscore is the sequencing library index.

In the study table, records such as:

```text
10_PM1_009.1
12_PM1_009.2
```

belong to the same patient because the letter code `PM1` is the patient identity. The final `.1`, `.2`, `.3` is the clinical time point. The intermediate numeric code is retained as source metadata but is not used to split one patient into different people.

CGIMU controls are independent individuals and keep their complete IDs.

## Input data

The new FASTQ batch defaults to:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

Legacy UDP0097-UDP0112 FASTQ files are optional and are used when `MEASLES_LEGACY_FASTQ_DIR` is set.

## Install and run

```bash
git clone https://github.com/komkovamariia/measles_2026.git
cd measles_2026
conda env create -f environment.yml
conda activate measles-tcr-2026
```

Pilot:

```bash
python scripts/run_pipeline.py --stage all --pilot --cpus 8 --threads-per-sample 2
```

Full Aldan-3 run:

```bash
sbatch scripts/run_slurm.sh
```

## Main patient-level outputs

For every patient with at least two time points the pipeline creates:

```text
figures/patients/<PATIENT>_cdr3_venn.png/pdf
figures/patients/<PATIENT>_T1_to_T2_trbv_change.png/pdf
figures/patients/<PATIENT>_T2_to_T3_trbv_change.png/pdf   # when available

results/longitudinal/patients/<PATIENT>_T1_to_T2_trbv_changes.csv
results/longitudinal/patients/<PATIENT>_T2_to_T3_trbv_changes.csv
```

The Venn diagram compares exact amino-acid CDR3β clonotypes between time points.

For every TRBV and patient transition the table contains:

```text
fraction_old
fraction_new
delta_fraction
delta_percentage_points
presence_change = gained / lost / persistent / absent_both
direction = increased / decreased / unchanged
```

TRBV fraction is calculated as:

```text
reads assigned to TRBV / all productive TRB reads in the sample
```

An individual patient has one repertoire measurement per time point, so a patient-specific p-value is not statistically identifiable. Patient tables therefore report the observed effect. Statistical significance is tested across biological patients.

## Group-level TRBV analysis

Two simple comparisons are performed separately for adults and children:

1. acute measles versus healthy controls: Mann-Whitney U on TRBV read fractions;
2. paired time point 1 → 2 and 2 → 3: Wilcoxon signed-rank test on TRBV read fractions.

Benjamini-Hochberg correction is applied across TRBV genes within each comparison.

Main group outputs:

```text
results/longitudinal/trbv_group_statistics.csv
figures/09_group_trbv_patterns.png/pdf
figures/10_group_trbv_gain_loss.png/pdf
```

## Other retained analyses

The previous simple analyses remain in the project:

- normalized repertoire richness and Shannon diversity;
- top-10 clonotype dominance;
- TRBV/TRBJ usage heatmaps;
- exact CDR3β matches to IEDB for measles, hepatitis B, poliovirus and tetanus;
- longitudinal diversity trajectories.

## Project layout

```text
config/sample_manifest.csv            source metadata
notebooks/00_prepare_and_mixcr.ipynb  FASTQ discovery, QC and MiXCR
notebooks/01_measles_analysis.ipynb   patient and group analyses
scripts/run_pipeline.py               notebook runner
scripts/run_slurm.sh                  one 24-CPU Slurm job
scripts/analysis_05_longitudinal_venn.py  Venn + TRBV changes
results/                              CSV tables
figures/                              PNG and PDF figures
audit_runs/                           executed notebooks
logs/                                 notebook and MiXCR logs
```

## Interpretation limits

The libraries are bulk TRB data from genomic-DNA multiplex sequencing. Read counts provide a relative clonotype-abundance signal and are sensitive to multiplex-PCR amplification. For this reason TRBV changes are reported as fractions of productive TRB reads and clone diversity is also assessed after common-depth downsampling.
