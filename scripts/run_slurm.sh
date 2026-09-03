#!/usr/bin/env bash
#SBATCH --job-name=measles_tcr
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --partition=medium
#SBATCH --constraint=hpc
#SBATCH --output=logs/slurm_%j.out

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT"
mkdir -p logs results figures audit_runs work
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
export MEASLES_LEGACY_FASTQ_DIR="${MEASLES_LEGACY_FASTQ_DIR:-/projects/mice_transplant_2025/measles/raw/pe150_res}"
export MEASLES_REQUIRE_ALL_INPUTS=1
export MEASLES_MIXCR_MEMORY_GB="${MEASLES_MIXCR_MEMORY_GB:-8}"

# The working mouse project did not use a conda MiXCR package. Its
# presenting_repseq.ipynb called ~/soft/mixcr/mixcr explicitly. Reuse that
# executable first so we inherit the already working installation/license.
if [[ -z "${MEASLES_MIXCR_BIN:-}" ]]; then
  if [[ -x "$HOME/soft/mixcr/mixcr" ]]; then
    export MEASLES_MIXCR_BIN="$HOME/soft/mixcr/mixcr"
  elif command -v mixcr >/dev/null 2>&1; then
    export MEASLES_MIXCR_BIN="$(command -v mixcr)"
  else
    echo "ERROR: MiXCR was not found. Expected $HOME/soft/mixcr/mixcr." >&2
    exit 1
  fi
fi

if [[ ! -d "$MEASLES_NEW_FASTQ_DIR" ]]; then
  echo "ERROR: new FASTQ directory not found: $MEASLES_NEW_FASTQ_DIR" >&2
  exit 1
fi
if [[ ! -d "$MEASLES_LEGACY_FASTQ_DIR" ]]; then
  echo "ERROR: legacy FASTQ directory not found: $MEASLES_LEGACY_FASTQ_DIR" >&2
  exit 1
fi

echo "Repository root: $ROOT"
echo "New FASTQ:    $MEASLES_NEW_FASTQ_DIR"
echo "Legacy FASTQ: $MEASLES_LEGACY_FASTQ_DIR"
echo "MiXCR binary: $MEASLES_MIXCR_BIN"
"$MEASLES_MIXCR_BIN" -v > results/slurm_mixcr_version.txt 2>&1
head -n 1 results/slurm_mixcr_version.txt

python scripts/validate_repository.py

# Fail fast on one tiny real subset before launching six MiXCR processes in
# parallel. This catches license, preset and basic input problems immediately.
python scripts/mixcr_smoke_test.py

python scripts/run_pipeline.py \
  --stage all \
  --resume \
  --cpus "${SLURM_CPUS_PER_TASK:-24}" \
  --threads-per-sample 4
