#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config" / "sample_manifest.csv"
CLONOTYPES = ROOT / "results" / "clonotypes"


def positive_data_row(path: Path) -> bool:
    """Require a real clonotype row with a positive count, not only a TSV header."""
    if not path.exists() or path.stat().st_size < 100:
        return False
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if not header:
                return False
            count_col = next(
                (c for c in ["cloneCount", "readCount", "count"] if c in header), None
            )
            aa_col = next(
                (c for c in ["aaSeqCDR3", "cdr3aa", "cdr3_aa"] if c in header), None
            )
            if count_col is None or aa_col is None:
                return False
            ci, ai = header.index(count_col), header.index(aa_col)
            for row in reader:
                if len(row) <= max(ci, ai):
                    continue
                try:
                    count = float(row[ci])
                except (TypeError, ValueError):
                    continue
                if count > 0 and str(row[ai]).strip():
                    return True
    except Exception:
        return False
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--remove-invalid", action="store_true")
    args = parser.parse_args()

    manifest = pd.read_csv(MANIFEST)
    included = manifest[
        manifest["include"].astype(str).str.lower().isin(["true", "1", "yes"])
    ].copy()

    invalid = []
    for lib in included["library_id"].astype(str):
        path = CLONOTYPES / f"{lib}.tsv"
        if path.exists() and not positive_data_row(path):
            invalid.append(lib)
            if args.remove_invalid:
                path.unlink()
                print(f"Removed invalid/header-only MiXCR export: {path}")

    if args.remove_invalid:
        if invalid:
            print(
                "Invalid exports removed so MiXCR resume will regenerate them:",
                ", ".join(invalid),
            )
        else:
            print("No invalid existing MiXCR exports found.")
        return

    missing_or_invalid = []
    for lib in included["library_id"].astype(str):
        path = CLONOTYPES / f"{lib}.tsv"
        if not positive_data_row(path):
            missing_or_invalid.append(lib)

    if missing_or_invalid:
        raise SystemExit(
            "CLONOTYPE EXPORT VALIDATION FAILED: included libraries without a positive productive TRB export: "
            + ", ".join(missing_or_invalid)
        )
    print(
        f"Clonotype export validation passed: {len(included)} included libraries have positive TRB exports."
    )


if __name__ == "__main__":
    main()
