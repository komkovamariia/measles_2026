# Validation notes

The original minimal project was checked locally with notebook schema validation, Python compilation and a mocked MiXCR run.

The current version additionally includes patient-level CDR3β Venn diagrams, TRBV read-fraction change tables and group-level TRBV statistics. These additions are wired into `notebooks/01_measles_analysis.ipynb` and run automatically with `--stage analysis` or `--stage all`.

Real biological results still require execution on the Aldan-3 FASTQ data.

Recommended order:

1. `git pull`
2. update/activate the conda environment
3. run the pilot
4. inspect MiXCR exports and metadata mapping
5. submit the full Slurm run
6. inspect `results/longitudinal/` and `figures/patients/`
