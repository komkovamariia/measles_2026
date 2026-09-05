# Analysis run guide

## Current analysis contract

The final metadata file is `config/measles_metadata_v2_age18_final.xlsx`

The executable metadata table is `config/analysis_metadata.csv`

Quality threshold is 4400 productive TRB clonotypes per sequencing library

Children are younger than 18 years. Adults are 18 years and older

Controls use `CONTROL` only

The main case time points are T1, T2 and T3

`new_106` and `new_111` are accepted T2 biological repeats for БА1_071. They are averaged after within library frequency normalization and count as one patient observation

## Recommended refresh on the cluster

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

After completion

```bash
sacct -j "$jid" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

Expected state is `COMPLETED` with exit code `0:0`

## Analysis order

The refresh job runs

```text
analysis_01_prepare.py
analysis_02_antigens.py
analysis_06_clinical_summary.py
analysis_07_longitudinal_expanded.py
analysis_08_final_reporting.py
analysis_09_patient_panels.py
publish_current_figures.py
```

`analysis_01_prepare.py` applies the final metadata, library QC and replicate aggregation

`analysis_06_clinical_summary.py` implements the cross sectional design using all accepted T1 observations

`analysis_07_longitudinal_expanded.py` implements paired T1 and T2 analysis and expanded clone fate

`analysis_08_final_reporting.py` reports expanded TRBV patterns and specificity reference coverage

`analysis_09_patient_panels.py` creates patient level figures using expanded clone overlap only

## Expected cohort counts

Cross sectional observations

```text
Children: CONTROL 4, T1 27, T2 8
Adults: CONTROL 13, T1 29, T2 19
```

Paired T1 and T2 patients

```text
Children: 6
Adults: 16
```

Final unique patients after QC

```text
Children: 33
Adults: 45
Total: 78
```

## Main output folders

```text
figures/group
figures/patients
results/current
```

## Important verification

```bash
python - <<'PY'
import pandas as pd
m = pd.read_csv('results/current/sample_metadata_current.csv')
m['analysis_include_final'] = m['analysis_include_final'].astype(str).str.lower().isin(['true','1','yes'])
print(m[m['library_id'].isin(['new_106','new_111'])][['library_id','patient_id','clinical_timepoint','productive_trb_clonotypes','analysis_include_final']])
print('\nControls with wrong time point:', len(m[m['is_control'].astype(str).str.lower().isin(['true','1','yes']) & ~m['clinical_timepoint'].eq('CONTROL')]))
PY
```

Both `new_106` and `new_111` must be included at T2
