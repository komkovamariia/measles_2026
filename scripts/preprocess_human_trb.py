from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import os
import re
import shutil
import subprocess

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
MANIFEST = ROOT / "config" / "sample_manifest.csv"
KNOWN_ISSUES = ROOT / "config" / "metadata_known_issues.csv"
NEW_DIR = Path(os.environ.get(
    "MEASLES_NEW_FASTQ_DIR",
    "/projects/mice_transplant_2025/raw/20260817_161454/Fastq",
)).expanduser()
legacy_value = os.environ.get("MEASLES_LEGACY_FASTQ_DIR", "").strip()
LEGACY_DIR = Path(legacy_value).expanduser() if legacy_value else None
PRESET = os.environ.get("MEASLES_MIXCR_PRESET", "milab-human-dna-tcr-multiplex")
CPUS = int(os.environ.get("MEASLES_CPUS", os.environ.get("SLURM_CPUS_PER_TASK", "1")))
THREADS_PER_SAMPLE = int(os.environ.get("MEASLES_THREADS_PER_SAMPLE", "4"))
MIXCR_MEMORY_GB = int(os.environ.get("MEASLES_MIXCR_MEMORY_GB", "8"))
PILOT = os.environ.get("MEASLES_PILOT", "0") == "1"
REQUIRE_ALL = os.environ.get("MEASLES_REQUIRE_ALL_INPUTS", "1") == "1" and not PILOT
EXPECTED_NEW_LANES = {"001", "002", "003", "004"}

RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
LOGS = ROOT / "logs"
INPUTS = ROOT / "work" / "inputs"
MIXCR_DIR = RESULTS / "mixcr"
CLONOTYPES = RESULTS / "clonotypes"
for p in [RESULTS, FIGURES, LOGS, INPUTS, MIXCR_DIR, CLONOTYPES]:
    p.mkdir(parents=True, exist_ok=True)


def resolve_mixcr() -> str:
    """Prefer the exact external MiXCR location used by presenting_repseq.ipynb."""
    configured = os.environ.get("MEASLES_MIXCR_BIN", "").strip()
    candidates: list[Path] = []
    if configured:
        candidates.append(Path(configured).expanduser())
    candidates.append(Path.home() / "soft" / "mixcr" / "mixcr")
    path_mixcr = shutil.which("mixcr")
    if path_mixcr:
        candidates.append(Path(path_mixcr))
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    raise RuntimeError(
        "MiXCR executable was not found. The mouse project used ~/soft/mixcr/mixcr. "
        "Set MEASLES_MIXCR_BIN explicitly if the executable is elsewhere."
    )


MIXCR_BIN = resolve_mixcr()

manifest = pd.read_csv(MANIFEST)
for col in ["include", "is_control"]:
    manifest[col] = manifest[col].astype(str).str.lower().isin(["true", "1", "yes"])
manifest["analysis_timepoint"] = pd.to_numeric(manifest["analysis_timepoint"], errors="coerce")

issues = pd.read_csv(KNOWN_ISSUES)
auto_issues: list[dict] = []

dup_mask = manifest.duplicated(["patient_id", "analysis_timepoint", "batch"], keep=False)
for _, r in manifest[dup_mask].iterrows():
    auto_issues.append({
        "library_id": r.library_id,
        "raw_sample_id": r.raw_sample_id,
        "issue": "duplicate_patient_timepoint",
        "note": "Multiple libraries map to the same patient/timepoint; retained for library QC and aggregated for biological analysis.",
    })
for _, r in manifest[manifest.age_group.eq("unknown")].iterrows():
    auto_issues.append({
        "library_id": r.library_id,
        "raw_sample_id": r.raw_sample_id,
        "issue": "unknown_age",
        "note": "Patient is absent from the current child/adult lists; retained and excluded from age-stratified inference.",
    })

NEW_RE = re.compile(r"^(?P<idx>\d+)_S\d+_L(?P<lane>\d{3})_R(?P<read>[12])_001\.fastq\.gz$")
expected_new = set(manifest.loc[manifest.batch.eq("new"), "sequencing_index"].astype(str))
observed_by_index: dict[str, dict[tuple[str, str], Path]] = {}
unexpected_files: list[str] = []
if NEW_DIR.exists():
    for path in NEW_DIR.glob("*.fastq.gz"):
        m = NEW_RE.match(path.name)
        if not m:
            unexpected_files.append(path.name)
            continue
        idx, lane, read = m.group("idx"), m.group("lane"), m.group("read")
        observed_by_index.setdefault(idx, {})[(lane, read)] = path
        if idx not in expected_new:
            unexpected_files.append(path.name)
else:
    auto_issues.append({
        "library_id": "NEW_BATCH",
        "raw_sample_id": str(NEW_DIR),
        "issue": "missing_new_fastq_directory",
        "note": "Configured new FASTQ directory does not exist.",
    })


def discover_new(index) -> tuple[list[tuple[Path, Path]], list[str]]:
    idx = str(index)
    found = observed_by_index.get(idx, {})
    pairs: list[tuple[Path, Path]] = []
    missing: list[str] = []
    for lane in sorted(EXPECTED_NEW_LANES):
        p1 = found.get((lane, "1"))
        p2 = found.get((lane, "2"))
        if p1 is None:
            missing.append(f"L{lane}_R1")
        if p2 is None:
            missing.append(f"L{lane}_R2")
        if p1 is not None and p2 is not None:
            pairs.append((p1, p2))
    return pairs, missing


def discover_legacy(key) -> list[tuple[Path, Path]]:
    """Find files such as FT..._L01_UDP0112_1.fq.gz and _2.fq.gz."""
    if LEGACY_DIR is None or not LEGACY_DIR.exists():
        return []
    r1_candidates: set[Path] = set()
    for pattern in [
        f"*{key}*_1.fq.gz",
        f"*{key}*_1.fastq.gz",
        f"*{key}*_R1*.fq.gz",
        f"*{key}*_R1*.fastq.gz",
    ]:
        r1_candidates.update(LEGACY_DIR.rglob(pattern))

    pairs: list[tuple[Path, Path]] = []
    for p1 in sorted(r1_candidates):
        s = str(p1)
        possible: list[str] = []
        if s.endswith("_1.fq.gz"):
            possible.append(s[:-8] + "_2.fq.gz")
        if s.endswith("_1.fastq.gz"):
            possible.append(s[:-11] + "_2.fastq.gz")
        possible.extend([s.replace("_R1_", "_R2_"), s.replace("R1", "R2")])
        p2 = next((Path(x) for x in possible if x != s and Path(x).exists()), None)
        if p2 is not None:
            pairs.append((p1, p2))
    return pairs


def link_pairs(library_id: str, pairs: list[tuple[Path, Path]]) -> None:
    for old in INPUTS.glob(f"{library_id}_L*_R*.fastq.gz"):
        old.unlink()
    for lane_i, (p1, p2) in enumerate(pairs, start=1):
        for read, src in [(1, p1), (2, p2)]:
            dst = INPUTS / f"{library_id}_L{lane_i:03d}_R{read}.fastq.gz"
            dst.symlink_to(src.resolve())


resolved_rows = []
fastq_audit_rows = []
for _, r in manifest.iterrows():
    if r.batch == "new":
        pairs, missing = discover_new(r.sequencing_index)
        ready = len(pairs) == 4 and not missing
        missing_text = ";".join(missing)
    else:
        pairs = discover_legacy(r.fastq_key)
        ready = len(pairs) >= 1
        missing_text = "" if ready else "legacy_read1_read2_pair"

    fastq_audit_rows.append({
        "library_id": r.library_id,
        "sequencing_index": r.sequencing_index,
        "batch": r.batch,
        "n_complete_pairs": len(pairs),
        "missing": missing_text,
        "fastq_ready": ready,
    })

    if not ready:
        auto_issues.append({
            "library_id": r.library_id,
            "raw_sample_id": r.raw_sample_id,
            "issue": "incomplete_new_fastq" if r.batch == "new" else "missing_legacy_fastq",
            "note": (
                f"Expected four paired lanes L001-L004; complete_pairs={len(pairs)}; missing={missing_text}"
                if r.batch == "new"
                else "No complete old _1/_2 or R1/R2 pair was found in MEASLES_LEGACY_FASTQ_DIR."
            ),
        })
    else:
        link_pairs(r.library_id, pairs)

    resolved_rows.append({**r.to_dict(), "fastq_ready": ready, "n_fastq_pairs": len(pairs)})

for name in unexpected_files:
    auto_issues.append({
        "library_id": "UNEXPECTED_FASTQ",
        "raw_sample_id": name,
        "issue": "unexpected_fastq",
        "note": "FASTQ is outside the current 96-record study manifest or has an unexpected filename.",
    })

resolved = pd.DataFrame(resolved_rows)
fastq_audit = pd.DataFrame(fastq_audit_rows)
all_issues = pd.concat([issues, pd.DataFrame(auto_issues)], ignore_index=True).drop_duplicates()
resolved.to_csv(RESULTS / "sample_manifest_resolved.csv", index=False)
fastq_audit.to_csv(RESULTS / "fastq_audit.csv", index=False)
all_issues.to_csv(RESULTS / "metadata_issues.csv", index=False)

print("Libraries in manifest:", len(resolved))
print("New FASTQ-ready:", int(resolved.loc[resolved.batch.eq("new"), "fastq_ready"].sum()), "/", int(resolved.batch.eq("new").sum()))
print("Legacy FASTQ-ready:", int(resolved.loc[resolved.batch.eq("legacy"), "fastq_ready"].sum()), "/", int(resolved.batch.eq("legacy").sum()))
print("Unexpected FASTQ files:", len(unexpected_files))

if REQUIRE_ALL:
    missing_required = resolved[resolved["include"] & ~resolved["fastq_ready"]]
    if not missing_required.empty:
        ids = ", ".join(missing_required.library_id.astype(str).tolist())
        raise RuntimeError(
            "Full analysis requires every included old and new sample from the manifest. "
            f"Missing/incomplete libraries: {ids}. See results/fastq_audit.csv and results/metadata_issues.csv."
        )

version = subprocess.run([MIXCR_BIN, "-v"], capture_output=True, text=True, check=True)
(RESULTS / "mixcr_version.txt").write_text(version.stdout + version.stderr, encoding="utf-8")
print("MiXCR binary:", MIXCR_BIN)
print((version.stdout or version.stderr).splitlines()[0])
print("MiXCR preset:", PRESET)

pilot_ids = ["new_10", "new_11", "new_12", "legacy_UDP0107", "legacy_UDP0108"]
run_table = resolved[resolved.fastq_ready].copy()
if PILOT:
    run_table = run_table[run_table.library_id.isin(pilot_ids)].copy()


def exported_table_is_valid(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 100:
        return False
    try:
        header = pd.read_csv(path, sep="\t", nrows=2)
    except Exception:
        return False
    has_aa = any(c in header.columns for c in ["aaSeqCDR3", "cdr3aa", "cdr3_aa"])
    has_count = any(c in header.columns for c in ["cloneCount", "readCount", "count"])
    return has_aa and has_count


def log_tail(path: Path, n: int = 60) -> str:
    if not path.exists():
        return "<log file was not created>"
    return "\n".join(path.read_text(encoding="utf-8", errors="replace").splitlines()[-n:])


def run_mixcr(row):
    lib = row.library_id
    final_tsv = CLONOTYPES / f"{lib}.tsv"
    if exported_table_is_valid(final_tsv):
        return lib, "skipped_valid_export"

    r1_pattern = str(INPUTS / f"{lib}_L{{{{n}}}}_R1.fastq.gz")
    r2_pattern = str(INPUTS / f"{lib}_L{{{{n}}}}_R2.fastq.gz")
    prefix = MIXCR_DIR / lib
    log_file = LOGS / f"mixcr_{lib}.log"

    cmd = [
        MIXCR_BIN,
        f"-Xmx{MIXCR_MEMORY_GB}g",
        "analyze",
        PRESET,
        "-f",
        "--threads", str(THREADS_PER_SAMPLE),
        r1_pattern,
        r2_pattern,
        str(prefix),
    ]
    with log_file.open("w", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(cmd) + "\n")
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{lib}: MiXCR analyze failed with exit code {proc.returncode}. "
            f"Last lines of {log_file}:\n{log_tail(log_file)}"
        )

    candidates = sorted(MIXCR_DIR.glob(f"{lib}*.clns"))
    if not candidates:
        raise RuntimeError(f"MiXCR finished for {lib}, but no .clns file was produced.")
    clns = candidates[0]

    export_cmd = [
        MIXCR_BIN,
        f"-Xmx{MIXCR_MEMORY_GB}g",
        "exportClones",
        "--dont-split-files",
        "--chains", "TRB",
        "--filter-out-of-frames",
        "--filter-stops",
        str(clns),
        str(final_tsv),
    ]
    with log_file.open("a", encoding="utf-8") as log:
        log.write("\nCOMMAND: " + " ".join(export_cmd) + "\n")
        proc = subprocess.run(export_cmd, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{lib}: MiXCR exportClones failed with exit code {proc.returncode}. "
            f"Last lines of {log_file}:\n{log_tail(log_file)}"
        )
    if not exported_table_is_valid(final_tsv):
        raise RuntimeError(f"MiXCR export for {lib} is missing required CDR3/count columns.")
    return lib, "done"


max_workers = max(1, CPUS // max(1, THREADS_PER_SAMPLE))
print(
    f"CPUs={CPUS}; threads/sample={THREADS_PER_SAMPLE}; "
    f"MiXCR JVM memory/sample={MIXCR_MEMORY_GB}G; concurrent MiXCR samples={max_workers}"
)
print("Libraries selected for MiXCR:", len(run_table))
statuses = {}
with ThreadPoolExecutor(max_workers=max_workers) as pool:
    futures = {pool.submit(run_mixcr, row): row.library_id for row in run_table.itertuples(index=False)}
    for future in as_completed(futures):
        lib = futures[future]
        try:
            _, status = future.result()
            statuses[lib] = status
            print(lib, status)
        except Exception as exc:
            statuses[lib] = f"FAILED: {exc}"
            print(statuses[lib])
            pd.DataFrame([
                {"library_id": k, "mixcr_status": v} for k, v in statuses.items()
            ]).to_csv(RESULTS / "mixcr_run_status.csv", index=False)
            raise

pd.DataFrame([
    {"library_id": k, "mixcr_status": v} for k, v in statuses.items()
]).to_csv(RESULTS / "mixcr_run_status.csv", index=False)


def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


qc_rows = []
for _, r in resolved.iterrows():
    path = CLONOTYPES / f"{r.library_id}.tsv"
    if not exported_table_is_valid(path):
        continue
    df = pd.read_csv(path, sep="\t")
    count_col = pick_col(df, ["cloneCount", "readCount", "count"])
    frac_col = pick_col(df, ["cloneFraction", "readFraction", "fraction"])
    aa_col = pick_col(df, ["aaSeqCDR3", "cdr3aa", "cdr3_aa"])
    counts = pd.to_numeric(df[count_col], errors="coerce").fillna(0)
    total = float(counts.sum())
    if frac_col:
        frac = pd.to_numeric(df[frac_col], errors="coerce").fillna(0)
        largest = float(frac.max()) if len(frac) else np.nan
    else:
        largest = float(counts.max() / total) if total > 0 and len(counts) else np.nan
    qc_rows.append({
        "library_id": r.library_id,
        "patient_id": r.patient_id,
        "analysis_timepoint": r.analysis_timepoint,
        "batch": r.batch,
        "productive_trb_reads": int(round(total)),
        "productive_trb_clonotypes": int(df[aa_col].notna().sum()),
        "largest_clone_fraction": largest,
    })

qc = pd.DataFrame(qc_rows)
qc.to_csv(RESULTS / "sample_qc.csv", index=False)


def exact_cdr3_set(lib):
    path = CLONOTYPES / f"{lib}.tsv"
    if not exported_table_is_valid(path):
        return set()
    df = pd.read_csv(path, sep="\t")
    aa = pick_col(df, ["aaSeqCDR3", "cdr3aa", "cdr3_aa"])
    return set(df[aa].dropna().astype(str))


if {"legacy_UDP0107", "legacy_UDP0108"}.issubset(set(qc.library_id)):
    s7, s8 = exact_cdr3_set("legacy_UDP0107"), exact_cdr3_set("legacy_UDP0108")
    union = s7 | s8
    udp = qc[qc.library_id.isin(["legacy_UDP0107", "legacy_UDP0108"])].copy()
    udp["cdr3_jaccard_0107_0108"] = len(s7 & s8) / len(union) if union else np.nan
    udp.to_csv(RESULTS / "udp0107_udp0108_qc.csv", index=False)
    print("UDP0107/UDP0108 exact CDR3β Jaccard:", udp["cdr3_jaccard_0107_0108"].iloc[0])

print("Preprocessing completed.")
