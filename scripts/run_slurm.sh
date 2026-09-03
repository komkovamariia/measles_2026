#!/usr/bin/env bash
#SBATCH --job-name=measles_tcr
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --partition=medium
#SBATCH --constraint=hpc
#SBATCH --output=logs/slurm_%j.out

set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p logs

eval "$(conda shell.bash hook)"
conda activate measles-tcr-2026

# Current 2026 FASTQ batch. Override only if the data are moved.
export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"

# Legacy UDP0097-UDP0112 data are optional for starting the new batch.
# If the variable was not supplied, try to locate UDP0107 under the study project.
if [[ -z "${MEASLES_LEGACY_FASTQ_DIR:-}" ]]; then
  legacy_probe="$(find /projects/mice_transplant_2025 -type f -name '*UDP0107*R1*.fastq.gz' -print -quit 2>/dev/null || true)"
  if [[ -n "$legacy_probe" ]]; then
    export MEASLES_LEGACY_FASTQ_DIR="$(dirname "$legacy_probe")"
    echo "Auto-detected legacy FASTQ directory: $MEASLES_LEGACY_FASTQ_DIR"
  else
    echo "Legacy UDP FASTQ directory was not found automatically. New batch will still run."
  fi
fi

python scripts/run_pipeline.py \
  --stage all \
  --resume \
  --cpus "${SLURM_CPUS_PER_TASK:-24}" \
  --threads-per-sample 4
