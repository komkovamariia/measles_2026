# Measles TCR beta repertoire analysis

Analysis of TRB repertoire changes after measles, with paired clone tracking and separate outputs for children and adults.

## Data and analysis

Reviewed sample assignments and QC decisions are in [config/analysis_metadata.csv](config/analysis_metadata.csv). Libraries require at least 4400 productive TRB clonotypes. Children are younger than 18 years. The accepted paired cohort contains 22 patients: 6 children and 16 adults.

T1 is early infection, T2 follows clinical recovery and T3 is about one month after recovery. The T2 libraries `new_106` and `new_111` belong to patient БА1_071; their frequencies are normalized separately before averaging.

Physical clone tracking uses CDR3 nucleotide sequence with TRBV and TRBJ. Functional comparisons use CDR3 amino acid sequence with the same genes. The primary expansion threshold is 0.10% of productive TRB reads. Sensitivity checks cover 0.05%, 0.10% and 0.20% through 1.50% in 0.10 percentage-point steps.

The primary branch samples each library without replacement to the patient's minimum accepted T1 and T2 read depth. It uses at least 500 iterations with default seed 1031. The complementary branch uses observed frequencies and Simpson D without downsampling. Both track whether T1-expanded clones remain expanded, fall below threshold or are not detected at T2. Zero T2 counts receive a depth-based detection bound without pseudocounts.

Depth correction cannot remove PCR amplification bias: these libraries have no UMI. Antigen annotation uses exact CDR3-beta amino-acid matches against the [IEDB TCRMatch table](https://github.com/IEDB/TCRMatch/blob/master/data/IEDB_data.tsv). The reference downloads automatically; set `MEASLES_IEDB_FILE` to use a fixed local copy. HLA typing is unavailable and candidate literature review remains manual.

## Run

From the repository root, create the environment once:

```bash
conda env create -f environment.yml
```

On the Slurm cluster, reuse the existing MiXCR exports and regenerate final results:

```bash
eval "$(conda shell.bash hook)"
conda activate measles-tcr-2026
mkdir -p logs
python scripts/validate_repository.py
jid=$(sbatch --parsable scripts/run_analysis_refresh_slurm.sh)
tail -F "logs/refresh_slurm_${jid}.out"
```

For raw FASTQ input, set `MEASLES_NEW_FASTQ_DIR` and `MEASLES_LEGACY_FASTQ_DIR` as shown in [config/example.env](config/example.env), then run `python scripts/preprocess_human_trb.py` before submitting the analysis. For 1000 rarefaction iterations, run `export MEASLES_RAREFACTION_ITERATIONS=1000` before submission.

## Results

- `figures/`: final PNG figures, with separate trees for adults and children; `FIGURE_INDEX.txt` lists the published plots.
- `results/final/candidate_clones.csv`: candidate T1-expanded clones from accepted paired patients.
- `results/supplementary_tables/`: full clone master, patient summaries and supporting analyses.
- `results/clonotypes/`: MiXCR exports; `results/cache/`: reusable computational tables.

Repository validation checks metadata and analysis settings. Validating the full analysis requires the sequencing inputs.
