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
mkdir -p logs results figures audit_runs results/cache
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
echo "Final library QC: at least 4400 productive TRB clonotypes"
echo "Age rule: child younger than 18 years, adult 18 years and older"
echo "Primary clone key: CDR3 nucleotide + TRBV + TRBJ"
echo "Functional/public clone key: CDR3 amino acid + TRBV + TRBJ"
echo "Primary T1 expanded threshold: 0.10 percent"
echo "Final sensitivity thresholds: 0.05 0.10 0.20 0.50 percent"
echo "Primary branch: pair-specific repeated rarefaction with ${PAIR_ITERATIONS} iterations"
echo "Complementary branch: no downsampling, article-style Simpson D and sqrt Simpson D"
echo "Completed MiXCR exports are reused"

python scripts/validate_repository.py

if [[ ! -s results/clonotypes/new_111.tsv ]]; then
  echo "new_111 clonotype export is missing. Resuming preprocessing for missing exports"
  export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"
  export MEASLES_LEGACY_FASTQ_DIR="${MEASLES_LEGACY_FASTQ_DIR:-/projects/mice_transplant_2025/measles/raw/pe150_res}"
  export MEASLES_REQUIRE_ALL_INPUTS=1
  python scripts/preprocess_human_trb.py
fi

# Restore reusable CSV caches as symlinks if the previous clean run moved them to results/cache.
for f in combined_trb_clonotypes.csv diversity_metrics.csv library_trb_clonotypes.csv library_trb_clonotypes_nt.csv; do
  if [[ ! -e "results/$f" && -s "results/cache/$f" ]]; then
    ln -s "cache/$f" "results/$f"
  fi
done

rm -rf \
  results/current results/final \
  figures/current figures/group figures/main figures/patients \
  figures/downsampling_pairwise figures/sensitivity_thresholds \
  figures/supplementary figures/rarefied figures/simpson_raw figures/validation
mkdir -p results/current figures/current

# Upstream established analyses create patient metadata, the existing Venn panels and
# the nucleotide master table scaffolding. Their exploratory figures are staging only.
python scripts/analysis_01_prepare.py
python scripts/analysis_02_antigens.py
python scripts/analysis_06_clinical_summary.py
python scripts/analysis_07_longitudinal_expanded.py
python scripts/analysis_08_final_reporting.py
python scripts/analysis_09_patient_panels.py

# Build the physical nucleotide clone master table. The legacy resampling block in this
# script is deliberately reduced because the literal pair-specific analysis follows next.
MEASLES_RAREFACTION_ITERATIONS=1 MEASLES_RAREFACTION_DEPTH=50000 \
  python scripts/analysis_10_clone_fate_v2.py

# Primary sampling-controlled branch.
export MEASLES_RAREFACTION_ITERATIONS="$PAIR_ITERATIONS"
python scripts/analysis_12_exact_rarefaction.py
python scripts/analysis_11_clone_enrichment.py

# Complementary branch requested from PMID 37051229: no downsampling, Simpson D,
# Simpson clonality = sqrt(Simpson D), normalized richness, then the same clone-fate logic.
python scripts/analysis_14_simpson_no_downsampling.py

# Verify the scientific contract before the final cleanup removes results/current.
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

master = pd.read_csv('results/current/paired_T1_clone_master_nt.csv')
if master.patient.nunique() != 22:
    raise SystemExit(f'Expected 22 paired patients, found {master.patient.nunique()}')
required = {
    'clone_id_nt','clone_id_aa','T1_fraction','T2_fraction','fate_0_10_percent',
    'P_detection_if_unchanged_empirical','pair_rarefaction_depth',
    'sampling_robust_not_detected','public_private','control_prevalence','motif_cluster','T3_fate'
}
missing = sorted(required - set(master.columns))
if missing:
    raise SystemExit(f'Master table is missing: {missing}')

rare = pd.read_csv('results/current/pairwise_rarefaction_summary.csv')
if rare.patient.nunique() != 22 or (pd.to_numeric(rare.iterations, errors='coerce') < 500).any():
    raise SystemExit('Pair-specific rarefaction is incomplete')

raw = pd.read_csv('results/current/simpson_raw_patient_metrics.csv')
if raw.patient.nunique() != 22 or set(raw.timepoint.astype(str)) != {'T1','T2'}:
    raise SystemExit('Simpson no-downsampling branch is incomplete')
for c in ['simpson_D','simpson_clonality','normalized_richness']:
    if c not in raw.columns:
        raise SystemExit(f'Simpson branch is missing {c}')

sens = pd.read_csv('results/current/simpson_raw_threshold_sensitivity.csv')
expected = {0.05,0.10,0.20,0.50}
observed = {round(float(x), 8) for x in sens.threshold_percent.dropna().unique()}
if observed != expected:
    raise SystemExit(f'Unexpected final raw sensitivity thresholds: {sorted(observed)}')

print('Pre-publication scientific contract passed')
PY

# Create only the figures and tables that map to the requested presentation logic.
python scripts/analysis_15_tz_final_reporting.py

# Final filesystem contract. results/current and figures/current must no longer exist.
python - <<'PY'
from pathlib import Path
import pandas as pd

if Path('results/current').exists():
    raise SystemExit('Temporary results/current was not cleaned')
if Path('figures/current').exists():
    raise SystemExit('Temporary figures/current was not cleaned')

required_tables = [
    'results/final/01_clone_master.csv',
    'results/final/02_candidates.csv',
    'results/final/03_patient_summary.csv',
    'results/final/04_threshold_sensitivity.csv',
    'results/final/05_feature_enrichment.csv',
    'results/final/06_validation.csv',
]
for f in required_tables:
    p=Path(f)
    if not p.exists() or p.stat().st_size == 0:
        raise SystemExit(f'Missing final table: {f}')

master = pd.read_csv('results/final/01_clone_master.csv')
if master.patient.nunique() != 22:
    raise SystemExit('Final clone master does not contain 22 patients')

required_figs = [
    'figures/FIGURE_INDEX.txt',
    'figures/rarefied/01_shared_and_fate.png',
    'figures/rarefied/02_sampling_and_sensitivity.png',
    'figures/rarefied/03_candidate_dumbbell_and_trbv.png',
    'figures/rarefied/04_enrichment_and_specificity.png',
    'figures/simpson_raw/01_simpson_metrics_and_fate.png',
    'figures/simpson_raw/02_shared_and_sensitivity.png',
    'figures/simpson_raw/03_candidate_dumbbell_and_trbv.png',
    'figures/simpson_raw/04_enrichment_and_specificity.png',
    'figures/validation/02_branch_concordance.png',
]
for f in required_figs:
    p=Path(f)
    if not p.exists() or p.stat().st_size == 0:
        raise SystemExit(f'Missing final figure: {f}')

n_patient = len(list(Path('figures/patients').glob('*.png')))
if n_patient < 22:
    raise SystemExit(f'Expected at least 22 patient Venn panels, found {n_patient}')

print('Final curated output contract passed')
PY

echo "Refresh completed successfully"
echo "Patient Venn panels: $ROOT/figures/patients"
echo "Primary rarefied branch: $ROOT/figures/rarefied"
echo "No-downsampling Simpson branch: $ROOT/figures/simpson_raw"
echo "Validation: $ROOT/figures/validation"
echo "Final tables only: $ROOT/results/final"
echo "Reusable computational caches: $ROOT/results/cache"
