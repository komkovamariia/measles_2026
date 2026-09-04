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
    old_cwd = Path.cwd()
    os.environ.update(extra_env)
    os.chdir(ROOT)
    started = time.time()

    with log_path.open("a", encoding="utf-8") as log:
        def say(message: str) -> None:
            line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} | {message}"
            print(line, flush=True)
            log.write(line + "\n")
            log.flush()

        say(f"START {path.name}")
        say(f"ROOT {ROOT}")
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
            os.chdir(old_cwd)
            for key, value in old_env.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


def require_files(paths: list[str], label: str) -> None:
    missing = []
    for rel in paths:
        path = ROOT / rel
        if not path.exists() or (path.is_file() and path.stat().st_size == 0):
            missing.append(rel)
    if missing:
        raise RuntimeError(f"{label} verification failed; missing: {missing}")


def verify_prepare_outputs() -> None:
    require_files([
        "results/sample_manifest_resolved.csv",
        "results/fastq_audit.csv",
        "results/sample_qc.csv",
        "results/mixcr_run_status.csv",
    ], "Preprocessing")

    clonotypes = list((ROOT / "results" / "clonotypes").glob("*.tsv"))
    if not clonotypes:
        raise RuntimeError("Preprocessing verification failed; no MiXCR clonotype TSV files were produced.")
    print(f"Preprocessing verification passed. Clonotype TSV files: {len(clonotypes)}")


def verify_full_outputs() -> None:
    required = [
        "results/sample_manifest_resolved.csv",
        "results/fastq_audit.csv",
        "results/sample_qc.csv",
        "results/raw_fastq_depth.csv",
        "results/combined_trb_clonotypes.csv",
        "results/diversity_metrics.csv",
        "results/antigen_reference_used.csv",
        "results/corrected/sample_metadata_enriched_with_analysis.csv",
        "results/corrected/patient_summary.csv",
        "results/corrected/cross_sectional_statistics.csv",
        "results/corrected/longitudinal_statistics.csv",
        "results/corrected/T1_T2_repertoire_similarity.csv",
        "results/corrected/TRBV_T1_T2_paired_statistics.csv",
        "results/corrected/corrected_vj_usage_by_sample.csv",
        "results/corrected/specificity_reference_coverage.csv",
        "results/corrected/MM1_033_identity_QC.csv",
        "results/corrected_age6/sample_metadata_age6.csv",
        "results/corrected_age6/age6_group_summary.csv",
        "results/corrected_age6/TRBV_priority_scores.csv",
        "results/corrected_age6/venn_T1_T2_vs_healthy_summary.csv",
        "results/corrected_age6/timepoint_age_audit.csv",
    ]
    require_files(required, "Corrected full output")

    corrected_figures = ROOT / "figures" / "corrected"
    required_prefixes = [f"{i:02d}_" for i in range(0, 16)]
    missing = []
    for prefix in required_prefixes:
        if not any(p.name.startswith(prefix) and p.suffix == ".png" for p in corrected_figures.glob("*.png")):
            missing.append(f"figures/corrected/{prefix}*.png")

    patient_venns = list((corrected_figures / "patients").glob("*_clinical_aav_venn.png"))
    if not patient_venns:
        missing.append("figures/corrected/patients/*_clinical_aav_venn.png")

    age6_figures = ROOT / "figures" / "corrected_age6"
    for name in [
        "TRBV_usage_adult_annotated.png",
        "TRBJ_usage_adult_annotated.png",
        "TRBV_usage_child_annotated.png",
        "TRBJ_usage_child_annotated.png",
        "TRBV_priority_scores_adult.png",
        "TRBV_priority_scores_child.png",
        "venn_summary_T1_T2_vs_healthy.png",
    ]:
        p = age6_figures / name
        if not p.exists() or p.stat().st_size == 0:
            missing.append(f"figures/corrected_age6/{name}")

    if missing:
        raise RuntimeError(f"Corrected full output verification failed; missing: {missing}")
    print(f"Corrected full output verification passed. Clinical patient Venn figures: {len(patient_venns)}")


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

    if args.stage in {"prepare", "all"}:
        execute_notebook(STAGES["prepare"], env)
        verify_prepare_outputs()

    if args.stage in {"analysis", "all"}:
        if args.stage == "analysis":
            verify_prepare_outputs()
        execute_notebook(STAGES["analysis"], env)

    if args.stage == "all" and not args.pilot:
        verify_full_outputs()


if __name__ == "__main__":
    main()
