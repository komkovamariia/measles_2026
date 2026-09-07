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
export MEASLES_RAREFACTION_DEPTH="${MEASLES_RAREFACTION_DEPTH:-50000}"
export MEASLES_RAREFACTION_ITERATIONS="${MEASLES_RAREFACTION_ITERATIONS:-500}"
export MEASLES_RAREFACTION_SEED="${MEASLES_RAREFACTION_SEED:-1031}"

echo "Repository root: $ROOT"
echo "Git commit: $(git rev-parse HEAD)"
echo "Time points: T1 T2 T3, controls use CONTROL"
echo "Age rule: child younger than 18 years, adult 18 years and older"
echo "Final library QC: at least 4400 productive TRB clonotypes"
echo "Established expanded clone threshold: 0.10 percent of productive TRB reads"
echo "Additive nucleotide clone fate analysis keeps the established aaV analysis unchanged"
echo "Nucleotide candidate thresholds: 0.05 0.10 0.20 0.50 0.75 1.00 percent"
echo "Requested depth control: ${MEASLES_RAREFACTION_DEPTH} productive reads, ${MEASLES_RAREFACTION_ITERATIONS} iterations"
echo "БА1_071 T2: new_106 and new_111 are averaged after within library normalization"
echo "Completed MiXCR exports are reused"

python scripts/validate_repository.py

# The resolved biological repeat new_111 must have a MiXCR export.
# Existing exports are reused. If this one export is missing, preprocessing resumes
# and MiXCR skips every already valid library.
if [[ ! -s results/clonotypes/new_111.tsv ]]; then
  echo "new_111 clonotype export is missing. Resuming preprocessing for missing exports"
  export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"
  export MEASLES_LEGACY_FASTQ_DIR="${MEASLES_LEGACY_FASTQ_DIR:-/projects/mice_transplant_2025/measles/raw/pe150_res}"
  export MEASLES_REQUIRE_ALL_INPUTS=1
  python scripts/preprocess_human_trb.py
else
  echo "Reuse existing new_111 clonotype export"
fi

rm -rf results/current figures/current figures/group figures/patients
mkdir -p results/current figures/current

python scripts/analysis_01_prepare.py
python scripts/analysis_02_antigens.py
python scripts/analysis_06_clinical_summary.py
python scripts/analysis_07_longitudinal_expanded.py
python scripts/analysis_08_final_reporting.py
python scripts/analysis_09_patient_panels.py
python scripts/analysis_10_clone_fate_v2.py
python scripts/analysis_11_clone_enrichment.py
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
    'P_detection_if_unchanged','sampling_robust_not_detected','public_private',
    'control_prevalence','motif_cluster','T3_fate','candidate_class'
}
missing_cols = sorted(required_master - set(master.columns))
if missing_cols:
    raise SystemExit(f'Nucleotide clone master table is missing columns: {missing_cols}')

required = [
    'figures/FIGURE_INDEX.txt',
    'figures/group/01_cross_sectional_clonality_children.png',
    'figures/group/01_cross_sectional_clonality_adults.png',
    'figures/group/04_paired_clonality_children.png',
    'figures/group/04_paired_clonality_adults.png',
    'figures/group/06_expanded_clone_fate_children.png',
    'figures/group/06_expanded_clone_fate_adults.png',
    'figures/group/14_clone_fate_v2_overview.png',
    'figures/group/17_repertoire_identity_fingerprint.png',
    'results/current/cross_sectional_statistics.csv',
    'results/current/longitudinal_expanded_clone_summary.csv',
    'results/current/new_T2_expansions.csv',
    'results/current/patient_panel_index.csv',
    'results/current/paired_T1_clone_master_nt.csv',
    'results/current/clone_fate_threshold_sensitivity.csv',
    'results/current/rarefaction_summary.csv',
    'results/current/repertoire_identity_qc.csv',
    'results/current/GA1_003_cross_batch_sensitivity.csv',
]
for name in required:
    p = Path(name)
    if not p.exists() or p.stat().st_size == 0:
        raise SystemExit(f'Missing refreshed output: {name}')
print('Final verification passed')
PY

echo "Refresh completed successfully"
echo "Group figures: $ROOT/figures/group"
echo "Patient figures: $ROOT/figures/patients"
echo "Tables: $ROOT/results/current"
