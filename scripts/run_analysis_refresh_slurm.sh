#!/usr/bin/env bash
#SBATCH --job-name=measles_refresh
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=medium
#SBATCH --constraint=hpc
#SBATCH --output=logs/refresh_slurm_%j.out

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT"
mkdir -p logs results figures audit_runs
export MPLBACKEND=Agg

if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "ERROR: conda was not found" >&2
  exit 1
fi
conda activate measles-tcr-2026
export MEASLES_REPO="$ROOT"
PAIR_ITERATIONS="${MEASLES_RAREFACTION_ITERATIONS:-500}"
if (( PAIR_ITERATIONS < 500 )); then
  echo "ERROR: MEASLES_RAREFACTION_ITERATIONS must be at least 500" >&2
  exit 1
fi
export MEASLES_RAREFACTION_SEED="${MEASLES_RAREFACTION_SEED:-1031}"

echo "Repository root: $ROOT"
echo "Git commit: $(git rev-parse HEAD)"
echo "Time points: T1 T2 T3, controls use CONTROL"
echo "Age rule: child younger than 18 years, adult 18 years and older"
echo "Final library QC: at least 4400 productive TRB clonotypes"
echo "Primary T1 expanded threshold: 0.10 percent of productive TRB reads"
echo "Sensitivity thresholds: 0.05 0.10 0.20 0.50 0.75 1.00 percent"
echo "Pair-specific repeated rarefaction: ${PAIR_ITERATIONS} iterations"
echo "Each pair uses its own minimum accepted T1 or T2 productive-read depth"
echo "Internal T1 sampling null provides empirical P_detection_if_unchanged"
echo "БА1_071 T2: new_106 and new_111 are rarefied separately, normalized, then averaged"
echo "Completed MiXCR exports are reused"

python scripts/validate_repository.py

if [[ ! -s results/clonotypes/new_111.tsv ]]; then
  echo "new_111 clonotype export is missing. Resuming preprocessing for missing exports"
  export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"
  export MEASLES_LEGACY_FASTQ_DIR="${MEASLES_LEGACY_FASTQ_DIR:-/projects/mice_transplant_2025/measles/raw/pe150_res}"
  export MEASLES_REQUIRE_ALL_INPUTS=1
  python scripts/preprocess_human_trb.py
else
  echo "Reuse existing new_111 clonotype export"
fi

rm -rf \
  results/current \
  figures/current \
  figures/group \
  figures/main \
  figures/patients \
  figures/downsampling_pairwise \
  figures/sensitivity_thresholds \
  figures/validation \
  figures/supplementary
mkdir -p results/current figures/current

# Established analyses are retained for reproducibility and selected supplementary context.
python scripts/analysis_01_prepare.py
python scripts/analysis_02_antigens.py
python scripts/analysis_06_clinical_summary.py
python scripts/analysis_07_longitudinal_expanded.py
python scripts/analysis_08_final_reporting.py
python scripts/analysis_09_patient_panels.py

# analysis_10 builds the nucleotide master table and sequence/T3/identity scaffolding.
# Its older global resampling block is retained in source for compatibility but reduced to
# one iteration here. The authoritative sampling analysis follows in analysis_12.
MEASLES_RAREFACTION_ITERATIONS=1 MEASLES_RAREFACTION_DEPTH=50000 \
  python scripts/analysis_10_clone_fate_v2.py

# Literal specification: each patient pair is repeatedly rarefied to its own common depth.
export MEASLES_RAREFACTION_ITERATIONS="$PAIR_ITERATIONS"
python scripts/analysis_12_exact_rarefaction.py

# Enrichment must use empirical sampling-aware candidate classes from analysis_12.
python scripts/analysis_11_clone_enrichment.py
python scripts/analysis_13_compact_reporting.py
python scripts/publish_current_figures.py

python - <<'PY'
from pathlib import Path
import pandas as pd

meta = pd.read_csv('results/current/sample_metadata_current.csv')
meta['analysis_include_final'] = meta['analysis_include_final'].astype(str).str.lower().isin(['true','1','yes'])
if meta['clinical_timepoint'].astype(str).eq('T0').any():
    raise SystemExit('T0 was reintroduced downstream')
controls = meta[meta['is_control'].astype(str).str.lower().isin(['true','1','yes'])]
if not controls['clinical_timepoint'].eq('CONTROL').all():
    raise SystemExit('A healthy control has a longitudinal time point')
rep = meta[meta['library_id'].isin(['new_106','new_111'])]
if len(rep) != 2 or not rep['analysis_include_final'].all() or not rep['clinical_timepoint'].eq('T2').all():
    raise SystemExit('БА1_071 T2 replicate handling is incorrect')

pairs = pd.read_csv('results/current/longitudinal_expanded_clone_summary.csv')
counts = pairs.groupby('age_group')['patient_id'].nunique().to_dict()
if counts.get('child', 0) != 6 or counts.get('adult', 0) != 16:
    raise SystemExit(f'Unexpected paired cohort counts: {counts}')

master = pd.read_csv('results/current/paired_T1_clone_master_nt.csv')
if master['patient'].nunique() != 22:
    raise SystemExit(f'Nucleotide clone master table has {master.patient.nunique()} paired patients instead of 22')
required_master = {
    'clone_id_nt','clone_id_aa','T1_fraction','T2_fraction','fate_0_10_percent',
    'P_detection_if_unchanged','P_detection_if_unchanged_empirical',
    'pair_rarefaction_depth','sampling_robust_not_detected','sampling_candidate_class',
    'public_private','control_prevalence','motif_cluster','T3_fate','candidate_class'
}
missing_cols = sorted(required_master - set(master.columns))
if missing_cols:
    raise SystemExit(f'Nucleotide clone master table is missing columns: {missing_cols}')

rare = pd.read_csv('results/current/pairwise_rarefaction_summary.csv')
if rare.patient.nunique() != 22:
    raise SystemExit(f'Pairwise rarefaction summary has {rare.patient.nunique()} patients instead of 22')
if (pd.to_numeric(rare.iterations, errors='coerce') < 500).any():
    raise SystemExit('Pairwise rarefaction used fewer than 500 iterations')
if (pd.to_numeric(rare.pair_rarefaction_depth, errors='coerce') <= 0).any():
    raise SystemExit('Pairwise rarefaction has a non-positive depth')

sens = pd.read_csv('results/current/pairwise_threshold_sensitivity.csv')
observed = {round(float(x), 8) for x in sens.threshold_percent.dropna().unique()}
expected = {0.05, 0.10, 0.20, 0.50, 0.75, 1.00}
if observed != expected:
    raise SystemExit(f'Unexpected threshold sensitivity set: {sorted(observed)}')

null = pd.read_csv('results/current/clone_sampling_null.csv')
if null.empty or null.patient.nunique() != 22:
    raise SystemExit('Empirical clone sampling null is incomplete')

required = [
    'figures/FIGURE_INDEX.txt',
    'figures/main/21_primary_fate_pairwise_rarefied.png',
    'figures/downsampling_pairwise/22_sampling_robust_clone_loss.png',
    'figures/sensitivity_thresholds/23_threshold_sensitivity_pairwise_rarefied.png',
    'figures/validation/26_validation_identity_T3_batch.png',
    'results/current/patient_panel_index.csv',
    'results/current/paired_T1_clone_master_nt.csv',
    'results/current/candidate_clones_compact.csv',
    'results/current/pairwise_rarefaction_summary.csv',
    'results/current/pairwise_threshold_sensitivity.csv',
    'results/current/clone_sampling_null.csv',
    'results/current/repertoire_identity_qc.csv',
    'results/current/GA1_003_cross_batch_sensitivity.csv',
    'results/current/extended_T3_clone_validation.csv',
    'results/current/extended_T3_validation_summary.csv',
]
for name in required:
    p = Path(name)
    if not p.exists() or p.stat().st_size == 0:
        raise SystemExit(f'Missing refreshed output: {name}')

patient_png = list(Path('figures/patients').glob('*.png'))
if len(patient_png) < 22:
    raise SystemExit(f'Expected at least 22 patient PNG panels, found {len(patient_png)}')
if Path('figures/group').exists():
    raise SystemExit('Obsolete figures/group directory was republished')
print('Final verification passed')
PY

echo "Refresh completed successfully"
echo "Main figures: $ROOT/figures/main"
echo "Patient figures: $ROOT/figures/patients"
echo "Pairwise downsampling: $ROOT/figures/downsampling_pairwise"
echo "Threshold sensitivity: $ROOT/figures/sensitivity_thresholds"
echo "Validation: $ROOT/figures/validation"
echo "Supplementary context: $ROOT/figures/supplementary"
echo "Tables: $ROOT/results/current"
