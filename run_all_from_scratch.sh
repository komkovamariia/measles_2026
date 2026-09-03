#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p logs results figures audit_runs resources work

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
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Updating existing conda environment: $ENV_NAME"
  conda env update -n "$ENV_NAME" -f environment.yml --prune
else
  echo "Creating conda environment: $ENV_NAME"
  conda env create -f environment.yml
fi
conda activate "$ENV_NAME"

export MEASLES_NEW_FASTQ_DIR="${MEASLES_NEW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/20260817_161454/Fastq}"
export MEASLES_LEGACY_FASTQ_DIR="${MEASLES_LEGACY_FASTQ_DIR:-/projects/mice_transplant_2025/measles/raw/pe150_res}"

if [[ ! -d "$MEASLES_NEW_FASTQ_DIR" ]]; then
  echo "ERROR: new FASTQ directory does not exist: $MEASLES_NEW_FASTQ_DIR" >&2
  exit 1
fi
if [[ ! -d "$MEASLES_LEGACY_FASTQ_DIR" ]]; then
  echo "ERROR: legacy FASTQ directory does not exist: $MEASLES_LEGACY_FASTQ_DIR" >&2
  exit 1
fi

# Match the working mice project: prefer the user's existing standalone MiXCR.
if [[ -z "${MEASLES_MIXCR_BIN:-}" ]]; then
  if [[ -x "$HOME/soft/mixcr/mixcr" ]]; then
    export MEASLES_MIXCR_BIN="$HOME/soft/mixcr/mixcr"
  elif command -v mixcr >/dev/null 2>&1; then
    export MEASLES_MIXCR_BIN="$(command -v mixcr)"
  else
    echo "ERROR: MiXCR not found. Expected $HOME/soft/mixcr/mixcr." >&2
    exit 1
  fi
fi

echo "New FASTQ:    $MEASLES_NEW_FASTQ_DIR"
echo "Legacy FASTQ: $MEASLES_LEGACY_FASTQ_DIR"
echo "MiXCR binary: $MEASLES_MIXCR_BIN"
"$MEASLES_MIXCR_BIN" -v > results/bootstrap_mixcr_version.txt 2>&1
head -n 1 results/bootstrap_mixcr_version.txt

python scripts/validate_repository.py

python - <<'PY'
from pathlib import Path
import requests
p = Path("resources/IEDB_data.tsv")
if not p.exists() or p.stat().st_size < 100:
    url = "https://raw.githubusercontent.com/IEDB/TCRMatch/master/data/IEDB_data.tsv"
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    p.write_bytes(r.content)
print(f"IEDB reference ready: {p} ({p.stat().st_size} bytes)")
PY

export MEASLES_REQUIRE_ALL_INPUTS=1
JOB_ID="$(sbatch --parsable scripts/run_slurm.sh)"

echo "Submitted full measles analysis as Slurm job: $JOB_ID"
echo "Bootstrap log: $ROOT/logs/bootstrap_submit.log"
echo "Slurm log after the job starts: $ROOT/logs/slurm_${JOB_ID}.out"
