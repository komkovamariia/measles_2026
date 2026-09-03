# Run guide

## Patient identity rule

Use the complete biological code after the sequencing index and before `.1/.2/.3`.

```text
15_CA1_013.1
19_CA1_013.2
```

are one patient `CA1_013`.

```text
22_CA1_017.1
```

is another patient `CA1_017`.

Thus `CA1_013` and `CA1_017` are never merged. The first integer is the sequencing-library index. The final integer after the dot is the time point.

## Run everything unattended from scratch

The repository contains `run_all_from_scratch.sh`. It creates the conda environment on the first run, updates it on later runs, and immediately submits the complete Slurm analysis. It does not perform a pilot.

From a server account where the repository has not yet been cloned, paste the following block once:

```bash
cd ~
if [ -d "$HOME/measles_2026/.git" ]; then
    cd "$HOME/measles_2026"
    git pull
else
    git clone https://github.com/komkovamariia/measles_2026.git "$HOME/measles_2026"
    cd "$HOME/measles_2026"
fi
mkdir -p logs
nohup bash run_all_from_scratch.sh > logs/bootstrap_submit.log 2>&1 &
```

After the `nohup` command you can disconnect from SSH. Environment creation/update continues in the background. When it finishes, the script submits the complete computation through Slurm, so the heavy MiXCR analysis runs on the cluster allocation.

The configured new-data path is:

```text
/projects/mice_transplant_2025/raw/20260817_161454/Fastq
```

The full Slurm job requests 24 CPUs. MiXCR receives 4 threads per library, allowing up to 6 independent libraries to be processed simultaneously.

If you later want to see what happened:

```bash
cd ~/measles_2026
cat logs/bootstrap_submit.log
squeue -u "$USER"
```

Once the job ID is visible, the Slurm output is:

```text
logs/slurm_<JOBID>.out
```

## Manual alternative

If the environment already exists and you only want to resubmit the complete calculation:

```bash
cd ~/measles_2026
git pull
conda activate measles-tcr-2026
sbatch scripts/run_slurm.sh
```

Existing non-empty MiXCR TRB exports are skipped automatically, so rerunning the Slurm job can continue an interrupted calculation.

## Main patient-level results

Each patient with at least two time points receives a Venn diagram:

```text
figures/patients/<PATIENT>_cdr3_venn.png
figures/patients/<PATIENT>_cdr3_venn.pdf
```

For two time points the plot is a 2-set Venn. For three time points it is a 3-set Venn. Sets contain exact CDR3β amino-acid clonotypes.

For each consecutive time-point transition the pipeline creates a TRBV read-share table and figure:

```text
results/longitudinal/patients/<PATIENT>_T1_to_T2_trbv_changes.csv
figures/patients/<PATIENT>_T1_to_T2_trbv_change.png/pdf
```

and, when available:

```text
results/longitudinal/patients/<PATIENT>_T2_to_T3_trbv_changes.csv
figures/patients/<PATIENT>_T2_to_T3_trbv_change.png/pdf
```

Important columns are:

```text
trbv
fraction_old
fraction_new
delta_percentage_points
presence_change
direction
```

`gained` means zero reads for that TRBV at the previous time point and non-zero reads at the new point. `lost` means the reverse. `persistent` means the TRBV is detected at both points.

## Group statistics

The statistical unit is a biological patient. Adults and children are analyzed separately.

```text
acute measles vs healthy control : Mann-Whitney U
T1 vs T2 within the same patient : paired Wilcoxon
T2 vs T3 within the same patient : paired Wilcoxon
```

The tested value is the fraction of all productive TRB reads assigned to each TRBV. Benjamini-Hochberg FDR correction is applied within each comparison.

Main outputs:

```text
results/longitudinal/trbv_group_statistics.csv
figures/09_group_trbv_patterns.png/pdf
figures/10_group_trbv_gain_loss.png/pdf
```

Other retained outputs include repertoire diversity, antigen-associated CDR3β matches, V/J usage and longitudinal diversity plots.
