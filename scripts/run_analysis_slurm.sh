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
mkdir -p logs

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
export MPLBACKEND=Agg

python scripts/validate_repository.py
rm -rf results/current figures/current figures/group figures/patients
mkdir -p results/current figures/current
python scripts/run_pipeline.py --stage analysis --cpus "${SLURM_CPUS_PER_TASK:-8}" --threads-per-sample 4

echo "Analysis completed"
echo "Figures: $ROOT/figures/group and $ROOT/figures/patients"
echo "Tables: $ROOT/results/current"
