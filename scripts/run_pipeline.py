#!/usr/bin/env python3
"""Notebook runner for the measles TRB project"""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
STAGES = {
    "prepare": ROOT / "notebooks" / "00_prepare_and_mixcr.ipynb",
    "analysis": ROOT / "notebooks" / "01_measles_analysis.ipynb",
}


def execute_notebook(path: Path, extra_env: dict[str, str]) -> None:
    output = ROOT / "audit_runs" / f"{path.stem}.executed.ipynb"
    log_path = ROOT / "logs" / f"{path.stem}.log"
    output.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    old_env = {k: os.environ.get(k) for k in extra_env}
    old_cwd = Path.cwd()
    os.environ.update(extra_env)
    os.chdir(ROOT)
    started = time.time()
    try:
        notebook = nbformat.read(path, as_version=4)
        client = NotebookClient(
            notebook,
            timeout=None,
            kernel_name="python3",
            resources={"metadata": {"path": str(ROOT)}},
            allow_errors=False,
        )
        client.execute()
        nbformat.write(notebook, output)
        with log_path.open("a", encoding="utf-8") as log:
            log.write(f"DONE {path.name} in {(time.time() - started) / 60:.1f} min\n")
    except Exception:
        if "notebook" in locals():
            nbformat.write(notebook, output)
        raise
    finally:
        os.chdir(old_cwd)
        for key, value in old_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def require_files(paths, label):
    missing = []
    for rel in paths:
        p = ROOT / rel
        if not p.exists() or (p.is_file() and p.stat().st_size == 0):
            missing.append(rel)
    if missing:
        raise RuntimeError(f"{label} verification failed: {missing}")


def verify_prepare_outputs():
    require_files(
        [
            "results/sample_manifest_resolved.csv",
            "results/fastq_audit.csv",
            "results/sample_qc.csv",
            "results/mixcr_run_status.csv",
        ],
        "Preprocessing",
    )
    if not list((ROOT / "results" / "clonotypes").glob("*.tsv")):
        raise RuntimeError("No MiXCR clonotype TSV files were produced")


def verify_analysis_outputs():
    require_files(
        [
            "results/library_trb_clonotypes.csv",
            "results/combined_trb_clonotypes.csv",
            "results/diversity_metrics.csv",
            "results/current/sample_metadata_current.csv",
            "results/current/cross_sectional_statistics.csv",
            "results/current/longitudinal_expanded_clone_summary.csv",
            "results/current/expanded_clone_fates.csv",
            "results/current/new_T2_expansions.csv",
            "results/current/expanded_TRBV_usage.csv",
            "results/current/patient_panel_index.csv",
            "figures/FIGURE_INDEX.txt",
        ],
        "Analysis",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage", choices=["prepare", "analysis", "all"], default="all"
    )
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument(
        "--cpus", type=int, default=int(os.getenv("SLURM_CPUS_PER_TASK", "1"))
    )
    parser.add_argument("--threads-per-sample", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    env = {
        "MEASLES_REPO": str(ROOT),
        "MEASLES_CPUS": str(args.cpus),
        "MEASLES_THREADS_PER_SAMPLE": str(args.threads_per_sample),
        "MEASLES_PILOT": "1" if args.pilot else "0",
    }
    if args.stage in {"prepare", "all"}:
        execute_notebook(STAGES["prepare"], env)
        verify_prepare_outputs()
    if args.stage in {"analysis", "all"}:
        if args.stage == "analysis":
            verify_prepare_outputs()
        execute_notebook(STAGES["analysis"], env)
        verify_analysis_outputs()


if __name__ == "__main__":
    main()
