#!/usr/bin/env python3
"""Count raw FASTQ read pairs per sequencing library without rerunning MiXCR.

The preprocessing stage creates work/inputs/<library>_Lxxx_R1/R2.fastq.gz symlinks.
For paired-end data, the number of R1 records equals the number of read pairs, so
only R1 files are decompressed for counting. raw_reads_total is reported as two
reads per pair (R1 + R2). This is an input-depth QC metric and is distinct from
productive_trb_reads after MiXCR.
"""
from __future__ import annotations

import argparse
import gzip
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "work" / "inputs"
RESULTS = ROOT / "results"


def fastq_records(path: Path) -> int:
    lines = 0
    with gzip.open(path, "rb") as fh:
        for _ in fh:
            lines += 1
    if lines % 4:
        raise RuntimeError(f"FASTQ line count is not divisible by 4: {path} ({lines} lines)")
    return lines // 4


def one_library(library_id: str) -> dict:
    r1 = sorted(INPUTS.glob(f"{library_id}_L*_R1.fastq.gz"))
    r2 = sorted(INPUTS.glob(f"{library_id}_L*_R2.fastq.gz"))
    if not r1 or len(r1) != len(r2):
        return {
            "library_id": library_id,
            "r1_files": len(r1),
            "r2_files": len(r2),
            "raw_read_pairs": None,
            "raw_reads_total": None,
            "compressed_bytes_all": None,
            "status": "missing_or_unpaired_fastq",
        }
    pairs = sum(fastq_records(p) for p in r1)
    compressed = sum(p.resolve().stat().st_size for p in r1 + r2)
    return {
        "library_id": library_id,
        "r1_files": len(r1),
        "r2_files": len(r2),
        "raw_read_pairs": int(pairs),
        "raw_reads_total": int(2 * pairs),
        "compressed_bytes_all": int(compressed),
        "status": "ok",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    manifest = pd.read_csv(RESULTS / "sample_manifest_resolved.csv")
    ids = manifest["library_id"].astype(str).tolist()
    rows = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(one_library, lib): lib for lib in ids}
        for fut in as_completed(futures):
            row = fut.result()
            rows.append(row)
            print(row["library_id"], row["status"], row["raw_read_pairs"])
    out = pd.DataFrame(rows).sort_values("library_id")
    path = RESULTS / "raw_fastq_depth.csv"
    out.to_csv(path, index=False)
    print(f"Wrote {path}")
    print("Libraries with raw read counts:", int(out.status.eq("ok").sum()), "/", len(out))


if __name__ == "__main__":
    main()
