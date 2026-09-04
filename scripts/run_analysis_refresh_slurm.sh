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

export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"
export MEASLES_LEGACY_FASTQ_DIR="${MEASLES_LEGACY_FASTQ_DIR:-/projects/mice_transplant_2025/measles/raw/pe150_res}"

echo "Repository root: $ROOT"
echo "Git commit: $(git rev-parse HEAD)"
echo "Age split: child < 6 years; adult >= 6 years"
echo "Final biological QC: >= 5000 productive TRB clonotypes/library"
echo "MiXCR will NOT be rerun."

python scripts/validate_repository.py

if [[ ! -s results/raw_fastq_depth.csv ]]; then
  python scripts/collect_fastq_depth.py --workers "${SLURM_CPUS_PER_TASK:-8}"
else
  echo "Reuse results/raw_fastq_depth.csv"
fi

# Rebuild notebook-derived biological tables and corrected metadata.
python scripts/run_pipeline.py \
  --stage analysis \
  --cpus "${SLURM_CPUS_PER_TASK:-8}" \
  --threads-per-sample 4

# Explicit final stage: this was missing from the earlier launcher path.
echo "Running final QC-filtered reporting..."
python scripts/analysis_08_final_reporting.py

# Publish authoritative final outputs into top-level figures/ for easy viewing.
echo "Publishing final figures to top-level figures/..."
python scripts/publish_current_figures.py

# Fail loudly if the user-facing outputs were not actually refreshed.
for f in \
  figures/CURRENT_OUTPUTS.txt \
  figures/00_library_QC_overview.png \
  figures/01_oligoclonality_case_control.png \
  figures/05_T1_T2_richness.png \
  figures/07_venn_and_weighted_similarity_summary.png \
  figures/08_TRBV_usage_adult.png \
  figures/09_TRBJ_usage_adult.png \
  "figures/patients/СА1_013_repertoire_QC_and_overlap.png" \
  results/current/sample_metadata_current.csv \
  results/current/patient_sample_QC_current.csv
do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: expected refreshed output is missing: $f" >&2
    exit 1
  fi
done

echo
echo "===== USER-FACING FIGURES ====="
find figures -maxdepth 1 -type f -printf '%TY-%Tm-%Td %TH:%TM:%TS  %p\n' | sort

echo
echo "===== CA1_013 ====="
ls -lh "figures/patients/СА1_013_repertoire_QC_and_overlap.png"

echo
echo "Refresh completed successfully."
echo "Open: $ROOT/figures"
echo "Tables: $ROOT/results/current"
