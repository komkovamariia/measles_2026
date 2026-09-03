# Measles TCR-beta repertoire analysis

A deliberately small Jupyter pipeline for the human measles TCR-beta project.

It answers four questions:

1. Are measles-associated CDR3β clonotypes more frequent during acute measles?
2. Are clonotypes associated with hepatitis B, poliovirus and tetanus less represented during acute measles?
3. Is the TCRβ repertoire less diverse during measles?
4. How do repertoire metrics change from admission to clinical recovery and later follow-up?

The workflow is intentionally simple: MiXCR preprocessing, two notebooks, simple non-parametric statistics and exact CDR3β matches to the curated IEDB TCRMatch table.

## Get the project on Aldan-3

This project is stored in the private repository `komkovamariia/measles_2026`.

```bash
cd ~
git clone https://github.com/komkovamariia/measles_2026.git
cd measles_2026
```

## Input data

The new FASTQ batch defaults to the current server location:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

Therefore it does not need to be supplied manually unless the files are moved.

Legacy UDP0097-UDP0112 FASTQ files are used when `MEASLES_LEGACY_FASTQ_DIR` is set. The Slurm script also tries to locate UDP0107 automatically under `/projects/mice_transplant_2025`.

## Install once

```bash
conda env create -f environment.yml
conda activate measles-tcr-2026
```

## Recommended first run

Pilot first:

```bash
python scripts/run_pipeline.py --stage all --pilot --cpus 8 --threads-per-sample 2
```

Then submit the complete analysis:

```bash
sbatch scripts/run_slurm.sh
```

Monitor with:

```bash
squeue -u "$USER"
tail -f logs/slurm_JOBID.out
```

MiXCR libraries with an existing valid TRB export are skipped automatically.

## Project layout

```text
config/sample_manifest.csv            sample metadata
notebooks/00_prepare_and_mixcr.ipynb  FASTQ discovery, metadata QC and MiXCR
notebooks/01_measles_analysis.ipynb   all four biological tasks
scripts/run_pipeline.py               notebook runner
scripts/run_slurm.sh                  one simple 24-CPU Slurm job
results/                              CSV tables
figures/                              PNG and PDF figures
audit_runs/                           executed notebooks
logs/                                 run and MiXCR logs
```

## Main figures

Exactly two result figures are created for each requested biological task:

| Task | Figures |
|---|---|
| 1. Measles response | measles-specific CDR3β group comparison; top-10 clone dominance |
| 2. Other infections | HBV/polio/tetanus group comparison; longitudinal recovery |
| 3. Repertoire diversity | richness/Shannon comparison; V/J usage heatmaps |
| 4. Recovery dynamics | paired richness; paired Shannon diversity |

Adults and children are analysed separately when metadata permit it. If either comparison group contains fewer than three biological samples, the notebook shows the data descriptively and does not report a p-value.

## Interpretation limits

The libraries are bulk TRB data from genomic-DNA multiplex sequencing. Read counts therefore provide a relative clonotype-abundance signal and are sensitive to multiplex-PCR amplification.

The primary antigen analysis uses exact amino-acid CDR3β matches to experimentally curated IEDB records. Database matches are evidence for antigen association, while absence of a match cannot be interpreted as absence of immunity.

Observed repertoire richness depends on sequencing depth, so the primary richness and Shannon comparison is performed after reproducible common-depth downsampling with seed `1031`.
