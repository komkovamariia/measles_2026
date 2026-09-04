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
echo "Final QC threshold: >= 5000 productive TRB clonotypes/library"
echo "This run does NOT rerun MiXCR."

python scripts/validate_repository.py

if [[ ! -s results/raw_fastq_depth.csv ]]; then
  python scripts/collect_fastq_depth.py --workers "${SLURM_CPUS_PER_TASK:-8}"
else
  echo "Reuse results/raw_fastq_depth.csv"
fi

# Force a true downstream refresh. Base MiXCR/clonotype outputs are preserved,
# but every derived corrected/current table and figure is deleted first.
echo "Removing stale downstream outputs before regeneration..."
rm -rf \
  results/corrected \
  results/corrected_age6 \
  results/current \
  figures/corrected \
  figures/corrected_age6 \
  figures/current
rm -f audit_runs/01_measles_analysis.executed.ipynb
mkdir -p results figures audit_runs

python scripts/run_pipeline.py \
  --stage analysis \
  --cpus "${SLURM_CPUS_PER_TASK:-8}" \
  --threads-per-sample 4

# Hard verification: the run is considered successful only if the new current
# outputs and the corrected CA1_013 patient figure were recreated.
for f in \
  results/current/sample_metadata_current.csv \
  results/current/cross_sectional_statistics_current.csv \
  results/current/T1_T2_paired_statistics_current.csv \
  results/current/venn_T1_T2_vs_healthy_current.csv \
  figures/current/00_library_QC_overview.png \
  figures/current/01_oligoclonality_case_control.png \
  figures/current/05_T1_T2_richness.png \
  figures/current/07_venn_and_weighted_similarity_summary.png \
  figures/current/08_TRBV_usage_adult.png \
  figures/current/09_TRBJ_usage_adult.png \
  "figures/current/patients/СА1_013_repertoire_QC_and_overlap.png"
do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: expected fresh output is missing: $f" >&2
    exit 1
  fi
done

cp figures/current/GENERATED_CURRENT.txt figures/CURRENT_OUTPUTS.txt

echo
echo "===== FRESH CURRENT FIGURES ====="
find figures/current -type f -printf '%TY-%Tm-%Td %TH:%TM:%TS  %p\n' | sort | tail -30

echo
echo "===== CURRENT RESULTS ====="
find results/current -maxdepth 1 -type f -printf '%TY-%Tm-%Td %TH:%TM:%TS  %p\n' | sort

echo
echo "Analysis-only run completed successfully."
echo "Authoritative figures: $ROOT/figures/current"
echo "Authoritative tables:  $ROOT/results/current"
echo "CA1_013 QC figure:     $ROOT/figures/current/patients/СА1_013_repertoire_QC_and_overlap.png"
