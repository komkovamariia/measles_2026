#!/usr/bin/env bash
#SBATCH --job-name=measles_refresh
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
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
  echo "ERROR: conda was not found." >&2
  exit 1
fi
conda activate measles-tcr-2026

export MEASLES_REPO="$ROOT"
export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"
export MEASLES_LEGACY_FASTQ_DIR="${MEASLES_LEGACY_FASTQ_DIR:-/projects/mice_transplant_2025/measles/raw/pe150_res}"

echo "Repository root: $ROOT"
echo "Git commit: $(git rev-parse HEAD)"
echo "Age split: child < 6 years; adult >= 6 years"
echo "Final biological QC: >= 5000 productive TRB clonotypes/library"
echo "MiXCR will NOT be rerun."

echo "Validating repository..."
python scripts/validate_repository.py

if [[ ! -s results/raw_fastq_depth.csv ]]; then
  echo "Collecting raw FASTQ depth..."
  python scripts/collect_fastq_depth.py --workers "${SLURM_CPUS_PER_TASK:-8}"
else
  echo "Reuse results/raw_fastq_depth.csv"
fi

echo "Removing stale derived outputs..."
rm -rf \
  results/corrected \
  results/corrected_age6 \
  results/current \
  figures/corrected \
  figures/corrected_age6 \
  figures/current \
  figures/group \
  figures/patients
find figures -maxdepth 1 -type f \( -name '*.png' -o -name '*.pdf' -o -name 'CURRENT_OUTPUTS.txt' -o -name 'FIGURE_INDEX.txt' \) -delete
mkdir -p results figures

echo "Running analysis_01_prepare.py..."
python scripts/analysis_01_prepare.py

echo "Running analysis_02_antigens.py..."
python scripts/analysis_02_antigens.py

echo "Running analysis_06_clinical_summary.py..."
python scripts/analysis_06_clinical_summary.py

echo "Running analysis_07_age6_reporting.py..."
python scripts/analysis_07_age6_reporting.py

echo "Running final reporting with read-only NumPy compatibility fix..."
python scripts/analysis_08_final_reporting_safe.py

echo "Running analysis_09_patient_panels.py..."
python scripts/analysis_09_patient_panels.py

echo "Publishing final figures into figures/group and figures/patients..."
python scripts/publish_current_figures.py

for f in \
  figures/FIGURE_INDEX.txt \
  figures/group/00_library_QC_overview.png \
  figures/group/01_oligoclonality_case_control.png \
  figures/group/05_T1_T2_richness.png \
  figures/group/07_venn_and_weighted_similarity_summary.png \
  figures/group/08_TRBV_usage_adult.png \
  figures/group/09_TRBJ_usage_adult.png \
  "figures/patients/СА1_013_repertoire_QC_and_overlap.png" \
  results/current/sample_metadata_current.csv \
  results/current/patient_panel_index.csv \
  results/current/patients/СА1_013_shared_clonotypes.csv
do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: expected refreshed output is missing: $f" >&2
    exit 1
  fi
done

echo "Removing intermediate figure directories..."
rm -rf figures/corrected figures/corrected_age6 figures/current

echo
echo "===== FINAL FIGURE TREE ====="
find figures -maxdepth 2 -type f -printf '%TY-%Tm-%Td %TH:%TM:%TS  %p\n' | sort

echo
echo "Refresh completed successfully."
echo "Group figures:   $ROOT/figures/group"
echo "Patient figures: $ROOT/figures/patients"
echo "Tables:          $ROOT/results/current"
