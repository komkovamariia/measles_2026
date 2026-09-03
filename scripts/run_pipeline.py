#!/usr/bin/env python3
"""Minimal notebook runner for the measles TCR project."""

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
    os.environ.update(extra_env)
    started = time.time()

    with log_path.open("a", encoding="utf-8") as log:
        def say(message: str) -> None:
            line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} | {message}"
            print(line, flush=True)
            log.write(line + "\n")
            log.flush()

        say(f"START {path.name}")
        notebook = nbformat.read(path, as_version=4)
        client = NotebookClient(
            notebook,
            timeout=None,
            kernel_name="python3",
            resources={"metadata": {"path": str(ROOT)}},
            allow_errors=False,
        )
        try:
            client.execute()
            nbformat.write(notebook, output)
            say(f"DONE {path.name} in {(time.time() - started) / 60:.1f} min")
        except Exception:
            nbformat.write(notebook, output)
            say(f"FAILED {path.name}; partial executed notebook saved to {output}")
            raise
        finally:
            for key, value in old_env.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


def verify_full_outputs() -> None:
    required = [
        "results/sample_manifest_resolved.csv",
        "results/fastq_audit.csv",
        "results/sample_qc.csv",
        "results/combined_trb_clonotypes.csv",
        "results/diversity_metrics.csv",
        "results/antigen_specificity_by_sample.csv",
        "results/longitudinal/patient_aav_overlap.csv",
        "results/longitudinal/trbv_patient_transition_summary.csv",
        "results/longitudinal/trbv_group_statistics.csv",
        "results/longitudinal/trbv_top_patterns.csv",
        "results/statistical_tests.csv",
        "results/result_summary.txt",
    ]
    required += [f"figures/{i:02d}_" for i in range(1, 11)]

    missing = []
    for item in required:
        if item.endswith("_"):
            prefix = Path(item).name
            parent = ROOT / Path(item).parent
            if not any(p.name.startswith(prefix) and p.suffix == ".png" for p in parent.glob("*.png")):
                missing.append(item + "*.png")
        else:
            path = ROOT / item
            if not path.exists() or (path.is_file() and path.stat().st_size == 0):
                missing.append(item)

    patient_venns = list((ROOT / "figures" / "patients").glob("*_aav_venn.png"))
    if not patient_venns:
        missing.append("figures/patients/*_aav_venn.png")

    if missing:
        raise RuntimeError(f"Output verification failed; missing: {missing}")
    print(f"Output verification passed. Patient Venn figures: {len(patient_venns)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["prepare", "analysis", "all"], default="all")
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--cpus", type=int, default=int(os.getenv("SLURM_CPUS_PER_TASK", "1")))
    parser.add_argument("--threads-per-sample", type=int, default=4)
    parser.add_argument("--resume", action="store_true", help="Valid MiXCR exports are skipped automatically.")
    args = parser.parse_args()

    if args.cpus < 1 or args.threads_per_sample < 1:
        raise SystemExit("CPU values must be positive.")

    env = {
        "MEASLES_REPO": str(ROOT),
        "MEASLES_CPUS": str(args.cpus),
        "MEASLES_THREADS_PER_SAMPLE": str(args.threads_per_sample),
        "MEASLES_PILOT": "1" if args.pilot else "0",
    }

    selected = ["prepare", "analysis"] if args.stage == "all" else [args.stage]
    for stage in selected:
        execute_notebook(STAGES[stage], env)

    if args.stage == "all" and not args.pilot:
        verify_full_outputs()


if __name__ == "__main__":
    main()
