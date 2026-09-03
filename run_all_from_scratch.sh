#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

mkdir -p logs

# Make conda available in non-interactive shells.
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

ENV_NAME="measles-tcr-2026"

# First run: create the environment. Later runs: update it to the repository definition.
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Updating existing conda environment: $ENV_NAME"
  conda env update -n "$ENV_NAME" -f environment.yml --prune
else
  echo "Creating conda environment: $ENV_NAME"
  conda env create -f environment.yml
fi

conda activate "$ENV_NAME"

# Current human measles FASTQ batch.
export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"

# Submit the complete analysis directly. No pilot run is performed.
JOB_ID="$(sbatch --parsable scripts/run_slurm.sh)"

echo "Submitted full measles analysis as Slurm job: $JOB_ID"
echo "Bootstrap log: $ROOT/logs/bootstrap_submit.log"
echo "Slurm log after the job starts: $ROOT/logs/slurm_${JOB_ID}.out"
