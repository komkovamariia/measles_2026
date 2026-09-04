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

## Recommended analysis-only refresh after MiXCR

Use this launcher after metadata/statistics/plot changes. It does **not** rerun MiXCR. It explicitly runs the final QC-filtered reporting stage and then republishes the authoritative figures into the top-level `figures/` directory.

```bash
cd ~/measles_2026
git pull --ff-only

eval "$(conda shell.bash hook)"
conda activate measles-tcr-2026
python scripts/validate_repository.py

jid=$(sbatch --parsable scripts/run_analysis_refresh_slurm.sh)
echo "JOB ID: $jid"
tail -F "logs/refresh_slurm_${jid}.out"
```

The age rule is fixed project-wide:

```text
child = age < 6 years
adult = age >= 6 years
```

The final biological library QC is:

```text
>= 5,000 productive TRB clonotypes per library
```

Libraries below this threshold remain visible in QC metadata but are excluded from final biological group, overlap and longitudinal analyses.

## Final output locations

User-facing figures are published directly into:

```text
figures/
figures/patients/
```

The same authoritative files are retained in:

```text
figures/current/
figures/current/patients/
```

Final tables are in:

```text
results/current/
```

The marker file

```text
figures/CURRENT_OUTPUTS.txt
```

confirms that the top-level `figures/` directory was refreshed from the current final analysis.

Important final figures include:

```text
figures/00_library_QC_overview.png
figures/01_oligoclonality_case_control.png
figures/02_richness_case_control.png
figures/03_shannon_case_control.png
figures/04_T1_T2_oligoclonality.png
figures/05_T1_T2_richness.png
figures/06_T1_T2_shannon.png
figures/07_venn_and_weighted_similarity_summary.png
figures/08_TRBV_usage_adult.png
figures/08_TRBV_usage_child.png
figures/09_TRBJ_usage_adult.png
figures/09_TRBJ_usage_child.png
figures/10_specificity_reference_coverage.png
```

For `СА1_013`, the final patient figure is:

```text
figures/patients/СА1_013_repertoire_QC_and_overlap.png
```

It explicitly shows that clinical T1 exists but has 0 productive TRB clonotypes and is excluded, while T2 and T3 are compared biologically.

## Clinical metadata

Clinical phases are taken explicitly from curated metadata and are never inferred from the raw `.1/.2/.3` suffix.

Examples:

```text
СА1_013: .1=T1, .2=T2, .3=T3
ММ1_033: .1=T1, .2=T2, .3=T2
ЭГ1_042: .1=T0, .2=T1, .3=T2
```

`new_39 / ММ1_033.2` remains excluded from biological inference because repertoire identity QC is strongly suspicious.

## Data paths

New FASTQ:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

Legacy FASTQ:

```text
/projects/mice_transplant_2025/measles/raw/pe150_res
```

Legacy rules remain:

```text
UDP0097-UDP0105 = adult measles cases
UDP0106-UDP0112 = healthy adult controls
```

UDP0108 is processed for QC but excluded from biological inference.

## Monitoring

```bash
squeue -u "$USER"
```

After completion:

```bash
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
```

A successful refresh ends with:

```text
Refresh completed successfully.
Open: /home/mkomkova/measles_2026/figures
Tables: /home/mkomkova/measles_2026/results/current
```

## MiXCR license

MiXCR requires a valid license for full preprocessing. The analysis-only refresh above reuses existing clonotype TSV files and does not rerun MiXCR.
