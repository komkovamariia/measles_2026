# Analysis run guide

## Current analysis contract

The authoritative metadata is `config/analysis_metadata.csv`

Quality threshold is 4400 productive TRB clonotypes per sequencing library

Children are younger than 18 years. Adults are 18 years and older

Controls use `CONTROL` only

Measles time points are T1, T2 and T3. T0 does not exist in this study

`new_106` and `new_111` are accepted T2 biological repeats for БА1_071. Each library is handled separately through depth control and normalization, then the normalized frequencies are averaged into one biological T2 observation

The established aaV analyses remain reproducible. The primary clone-fate analysis uses all 22 accepted T1-T2 pairs and tracks nucleotide clonotypes

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

```text
analysis_01_prepare.py
analysis_02_antigens.py
analysis_06_clinical_summary.py
analysis_07_longitudinal_expanded.py
analysis_08_final_reporting.py
analysis_09_patient_panels.py
analysis_10_clone_fate_v2.py
analysis_12_exact_rarefaction.py
analysis_11_clone_enrichment.py
analysis_13_compact_reporting.py
publish_current_figures.py
```

`analysis_10_clone_fate_v2.py` builds the nucleotide master table, public/private annotation, sequence-neighbour scaffolding, T3 scaffolding and repertoire identity fingerprint. Its older global resampling block is retained only for backward reproducibility and is reduced to one iteration by the refresh launcher

`analysis_12_exact_rarefaction.py` is the authoritative sampling analysis

`analysis_11_clone_enrichment.py` is run after sampling correction, so TRBV, TRBJ and CDR3 motif enrichment compares sampling-robust loss or strong contraction against retained expanded clones, using patient-aware GEE adjusted for baseline T1 clone frequency

`analysis_13_compact_reporting.py` produces the compact final candidate table, extended T3 validation, current GA1_003 sensitivity and only figures that answer a distinct question

## Clone definitions

Physical tracking unit

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

For every T1-expanded clone the T2 fate is

```text
retained expanded
detected below threshold
not detected
```

The master table contains T1 and T2 counts, fractions, absolute change, log2 T2 over T1 when T2 is nonzero, contraction, share of the original T1 expanded mass, public/private status, healthy-control prevalence, available exact antigen annotation, sequence-neighbour information, T3 fate, batch information and final candidate class

A zero T2 value never receives an arbitrary pseudocount. The table contains a 95 percent upper detection bound instead

## Pair-specific repeated rarefaction

This follows the specification literally

For patient `i`

```text
pair_depth_i = minimum productive-read depth among every accepted T1 and T2 library contributing to that patient
```

Each accepted T1 and T2 sequencing library is sampled without replacement to `pair_depth_i`

For БА1_071 the two T2 libraries are rarefied separately. Their frequencies are averaged only after rarefaction and within-library normalization

The procedure is repeated at least 500 times. Default

```text
iterations: 500
seed: 1031
```

Optional

```bash
export MEASLES_RAREFACTION_ITERATIONS=1000
export MEASLES_RAREFACTION_SEED=1031
```

There is no global 50000-read target in the authoritative analysis. Different patients may therefore have different pair depths. This is intentional because the central question is within-patient T1 to T2 clone fate

Observed libraries are rarefied with multivariate hypergeometric sampling without replacement

The internal null starts from repeatedly rarefied T1. For every clone, the workflow records how often that clone would be detected in an unchanged pseudo-T2 with the same sampling effort. The resulting empirical column is

```text
P_detection_if_unchanged
```

The old analytical raw-depth probability is retained separately as a check

A T1-expanded clone that is absent in real T2 and has empirical `P_detection_if_unchanged >= 0.99` is flagged as sampling-robust not detected

The sampling model controls finite sequencing depth. It does not remove PCR amplification bias because these libraries have no UMI

## T3 and identity validation

Available follow-up relations are used without inventing missing points

```text
ЭГ1_042: T1 T2 T3
СА1_013: valid T2 T3 because T1 failed QC
ММ1_033: valid T1 T3 because T2 failed identity QC
```

For ЭГ1_042 the workflow can explicitly identify a T1 clone absent at T2 that reappears at T3

Repertoire identity QC uses all-vs-all Morisita-Horn similarity. GA1_003 is additionally summarized with and without the cross-batch pair

## Nucleotide export requirement

Existing MiXCR TSV files must contain `nSeqCDR3`

```bash
python - <<'PY'
import pandas as pd
x = pd.read_csv('results/clonotypes/new_101.tsv', sep='\t', nrows=1)
print([c for c in x.columns if 'CDR3' in c])
PY
```

If `nSeqCDR3` is absent, re-export the existing `.clns` file. FASTQ alignment and assembly do not need to be repeated

## Expected cohort counts

```text
Cross sectional children: CONTROL 4, T1 27, T2 8
Cross sectional adults: CONTROL 13, T1 29, T2 19
Paired T1-T2: 6 children and 16 adults, total 22
Final unique patients: 33 children and 45 adults, total 78
```

## Final figure folders

The publisher intentionally keeps the tree shallow

```text
figures/main
figures/patients
figures/downsampling_pairwise
figures/sensitivity_thresholds
figures/validation
figures/supplementary
figures/FIGURE_INDEX.txt
```

`figures/main` contains the primary fate figure and only data-supported enrichment or sequence-neighbour figures

`figures/downsampling_pairwise` contains the empirical sampling-null conclusion figure

`figures/sensitivity_thresholds` contains one combined plot for all six expansion thresholds

`figures/validation` contains one combined identity, T3 and GA1_003 sensitivity figure

`figures/supplementary` contains only selected established context. Old TRBV heatmaps and duplicate exploratory plots are not republished

## Main tables

```text
paired_T1_clone_master_nt.csv
candidate_clones_compact.csv
pairwise_rarefaction_summary.csv
pairwise_threshold_sensitivity.csv
clone_sampling_null.csv
lost_vs_retained_feature_enrichment.csv
sequence_neighbour_cluster_summary.csv
extended_T3_clone_validation.csv
extended_T3_validation_summary.csv
repertoire_identity_morisita_matrix.csv
repertoire_identity_qc.csv
GA1_003_cross_batch_sensitivity.csv
```

## Quick verification after a completed job

```bash
python - <<'PY'
import pandas as pd
x = pd.read_csv('results/current/paired_T1_clone_master_nt.csv')
r = pd.read_csv('results/current/pairwise_rarefaction_summary.csv')
s = pd.read_csv('results/current/pairwise_threshold_sensitivity.csv')
print('Paired patients:', x.patient.nunique())
print('Rarefaction patients:', r.patient.nunique())
print('Iterations:', sorted(r.iterations.unique()))
print('Pair depth range:', int(r.pair_rarefaction_depth.min()), int(r.pair_rarefaction_depth.max()))
print('Thresholds percent:', sorted(s.threshold_percent.unique()))
print('Primary expanded candidates:', x.expanded_T1_0_10_percent.astype(str).str.lower().isin(['true','1','yes']).sum())
print('Sampling robust losses:', x.sampling_robust_not_detected.astype(str).str.lower().isin(['true','1','yes']).sum())
PY
```
