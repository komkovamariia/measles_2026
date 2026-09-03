#!/usr/bin/env bash
#SBATCH --job-name=measles_tcr
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --partition=medium
#SBATCH --constraint=hpc
#SBATCH --output=logs/slurm_%j.out

set -euo pipefail

# Slurm executes a copied job script from its spool directory, so $0 does not
# point to the repository. SLURM_SUBMIT_DIR is the directory from which sbatch
# was called and is therefore the correct repository root for this project.
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
  echo "ERROR: conda was not found in the Slurm shell." >&2
  exit 1
fi
conda activate measles-tcr-2026

export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"
export MEASLES_REQUIRE_ALL_INPUTS=1

if [[ -z "${MEASLES_LEGACY_FASTQ_DIR:-}" ]]; then
  legacy_named="$(find /projects/mice_transplant_2025 -type d -name 'pe150_res' -print -quit 2>/dev/null || true)"
  if [[ -n "$legacy_named" ]]; then
    export MEASLES_LEGACY_FASTQ_DIR="$legacy_named"
  else
    legacy_probe="$(find /projects/mice_transplant_2025 -type f -name '*UDP0107*R1*.fastq.gz' -print -quit 2>/dev/null || true)"
    if [[ -n "$legacy_probe" ]]; then
      export MEASLES_LEGACY_FASTQ_DIR="$(dirname "$legacy_probe")"
    fi
  fi
fi

if [[ -z "${MEASLES_LEGACY_FASTQ_DIR:-}" || ! -d "$MEASLES_LEGACY_FASTQ_DIR" ]]; then
  echo "ERROR: legacy UDP0097-UDP0112 FASTQ directory was not found." >&2
  echo "Set MEASLES_LEGACY_FASTQ_DIR to the old pe150_res directory and resubmit." >&2
  exit 1
fi

echo "Repository root: $ROOT"
echo "New FASTQ:    $MEASLES_NEW_FASTQ_DIR"
echo "Legacy FASTQ: $MEASLES_LEGACY_FASTQ_DIR"
mixcr -v > results/slurm_mixcr_version.txt 2>&1
head -n 1 results/slurm_mixcr_version.txt
python scripts/validate_repository.py

python scripts/run_pipeline.py \
  --stage all \
  --resume \
  --cpus "${SLURM_CPUS_PER_TASK:-24}" \
  --threads-per-sample 4
