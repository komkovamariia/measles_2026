#!/usr/bin/env bash
#SBATCH --job-name=measles_analysis
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --constraint=hpc
#SBATCH --output=logs/analysis_slurm_%j.out

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
echo "Environment: ${CONDA_DEFAULT_ENV:-unknown}"
echo "Python: $(command -v python)"
echo "Age split: child < 6 years; adult >= 6 years"

python scripts/validate_repository.py

if [[ ! -s results/raw_fastq_depth.csv ]]; then
  python scripts/collect_fastq_depth.py --workers "${SLURM_CPUS_PER_TASK:-8}"
else
  echo "Reuse results/raw_fastq_depth.csv"
fi

python scripts/run_pipeline.py \
  --stage analysis \
  --resume \
  --cpus "${SLURM_CPUS_PER_TASK:-8}" \
  --threads-per-sample 4

echo "Analysis-only run completed successfully."
