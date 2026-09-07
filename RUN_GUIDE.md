# Analysis run guide

## Current analysis contract

The authoritative repository metadata is `config/analysis_metadata.csv`

Quality threshold is 4400 productive TRB clonotypes per sequencing library

Children are younger than 18 years. Adults are 18 years and older

Controls use `CONTROL` only

The main case time points are T1, T2 and T3

`new_106` and `new_111` are accepted T2 biological repeats for БА1_071. They are averaged after within library frequency normalization and count as one patient observation

The established aaV analyses remain in place. The additive clone fate analysis tracks nucleotide clonotypes across the 22 accepted T1-T2 pairs

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
analysis_10_clone_fate_v2.py
analysis_11_clone_enrichment.py
publish_current_figures.py
```

`analysis_01_prepare.py` applies the final metadata, library QC and replicate aggregation

`analysis_06_clinical_summary.py` implements the cross sectional design using all accepted T1 observations

`analysis_07_longitudinal_expanded.py` implements the established paired T1 and T2 aaV analysis and expanded clone fate

`analysis_08_final_reporting.py` reports expanded TRBV patterns, patient averaged overlap and specificity reference coverage

`analysis_09_patient_panels.py` creates the established patient level panels

`analysis_10_clone_fate_v2.py` adds nucleotide CDR3 tracking, complete T1 clone fate, threshold sensitivity, analytical sampling detection probabilities, repeated equal-depth resampling, public/private annotation, T3 checks and the repertoire identity fingerprint

`analysis_11_clone_enrichment.py` tests TRBV, TRBJ and CDR3 amino-acid motif enrichment among lost versus retained T1-expanded clones using patient-aware GEE adjusted for baseline T1 clone frequency

## Nucleotide clone fate definitions

Primary physical tracking unit

```text
CDR3 nucleotide + TRBV + TRBJ
```

Functional and public clone unit

```text
CDR3 amino acid + TRBV + TRBJ
```

Primary T1 candidate threshold

```text
0.10 percent of productive TRB reads
```

Sensitivity thresholds

```text
0.05 0.10 0.20 0.50 0.75 1.00 percent
```

For a primary T1-expanded clone the T2 fate is one of

```text
retained expanded
detected below 0.10 percent
not detected
```

A zero T2 value is not assigned a pseudocount. The master table contains the 95 percent upper detection bound and `P_detection_if_unchanged`

## Depth-control settings

Defaults

```text
requested depth: 50000 productive reads
iterations: 500
seed: 1031
```

They can be changed before submitting the job

```bash
export MEASLES_RAREFACTION_DEPTH=50000
export MEASLES_RAREFACTION_ITERATIONS=500
export MEASLES_RAREFACTION_SEED=1031
```

The workflow never upsamples. If the shallowest accepted paired library is below the requested depth, that lower observed depth becomes the common depth

Clone-level detection probability is calculated analytically. Iterations are used only for repertoire-level depth-control sensitivity, not to estimate a probability that has a direct model solution

The sampling model does not correct PCR amplification bias because the libraries have no UMI

## Nucleotide export requirement

The existing MiXCR TSV files are expected to contain `nSeqCDR3`

Check one export with

```bash
python - <<'PY'
import pandas as pd
x = pd.read_csv('results/clonotypes/new_101.tsv', sep='\t', nrows=1)
print([c for c in x.columns if 'CDR3' in c])
PY
```

If `nSeqCDR3` is absent in an old export, re-export the existing `.clns` file with MiXCR. It is not necessary to repeat FASTQ alignment or assembly

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
Total: 22
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

Important new tables

```text
paired_T1_clone_master_nt.csv
clone_fate_threshold_sensitivity.csv
rarefaction_summary.csv
sequence_neighbour_cluster_summary.csv
lost_vs_retained_feature_enrichment.csv
T3_clone_fate_validation.csv
repertoire_identity_morisita_matrix.csv
repertoire_identity_qc.csv
GA1_003_cross_batch_sensitivity.csv
```

Important new figures

```text
14_clone_fate_v2_overview.png
15_lost_vs_retained_feature_enrichment.png
16_sequence_neighbour_clusters.png
17_repertoire_identity_fingerprint.png
18_T3_clone_fate_validation.png
```

Figures 15, 16 and 18 are generated only when the relevant data support a non-empty plot. Existing figures 01 through 13 and all patient panels are retained

## Important verification

```bash
python - <<'PY'
import pandas as pd
m = pd.read_csv('results/current/sample_metadata_current.csv')
m['analysis_include_final'] = m['analysis_include_final'].astype(str).str.lower().isin(['true','1','yes'])
print(m[m['library_id'].isin(['new_106','new_111'])][['library_id','patient_id','clinical_timepoint','productive_trb_clonotypes','analysis_include_final']])
print('\nControls with wrong time point:', len(m[m['is_control'].astype(str).str.lower().isin(['true','1','yes']) & ~m['clinical_timepoint'].eq('CONTROL')]))

x = pd.read_csv('results/current/paired_T1_clone_master_nt.csv')
print('\nPaired patients in nucleotide master table:', x.patient.nunique())
print('Primary T1 expanded candidates:', x.expanded_T1_0_10_percent.astype(str).str.lower().isin(['true','1','yes']).sum())
print('Sampling robust losses:', x.sampling_robust_not_detected.astype(str).str.lower().isin(['true','1','yes']).sum())
PY
```

Both `new_106` and `new_111` must be included at T2 and the nucleotide master table must contain 22 paired patients
