#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p logs results figures audit_runs resources

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
if [[ ! -d "$MEASLES_NEW_FASTQ_DIR" ]]; then
  echo "ERROR: new FASTQ directory does not exist: $MEASLES_NEW_FASTQ_DIR" >&2
  exit 1
fi

# Locate the required old cohort before submission.
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
  echo "ERROR: old UDP0097-UDP0112 / pe150_res FASTQ directory was not found." >&2
  exit 1
fi

echo "New FASTQ:    $MEASLES_NEW_FASTQ_DIR"
echo "Legacy FASTQ: $MEASLES_LEGACY_FASTQ_DIR"

mixcr -v | head -n 1
python scripts/validate_repository.py

# Cache the small IEDB TCRMatch reference on the login node so the compute job
# does not depend on outbound web access for antigen annotation.
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

# Exported variables are inherited by sbatch in the normal Slurm configuration.
export MEASLES_REQUIRE_ALL_INPUTS=1
JOB_ID="$(sbatch --parsable scripts/run_slurm.sh)"

echo "Submitted full measles analysis as Slurm job: $JOB_ID"
echo "Bootstrap log: $ROOT/logs/bootstrap_submit.log"
echo "Slurm log after the job starts: $ROOT/logs/slurm_${JOB_ID}.out"
