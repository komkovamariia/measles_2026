#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
NEW_DIR = Path(os.environ.get(
    "MEASLES_NEW_FASTQ_DIR",
    "/projects/mice_transplant_2025/raw/20260817_161454/Fastq",
)).expanduser()
PRESET = os.environ.get("MEASLES_MIXCR_PRESET", "milab-human-dna-tcr-multiplex")
LOG = ROOT / "logs" / "mixcr_smoke.log"
OUT_DIR = ROOT / "work" / "mixcr_smoke"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)


def resolve_mixcr() -> str:
    configured = os.environ.get("MEASLES_MIXCR_BIN", "").strip()
    candidates = []
    if configured:
        candidates.append(Path(configured).expanduser())
    path_mixcr = shutil.which("mixcr")
    if path_mixcr:
        candidates.append(Path(path_mixcr))
    candidates.append(Path.home() / "soft" / "mixcr" / "mixcr")
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    raise RuntimeError(
        "MiXCR executable was not found after activating the project environment. "
        "Set MEASLES_MIXCR_BIN explicitly only if a different installation is required."
    )


mixcr = resolve_mixcr()
rx = re.compile(r"^(?P<idx>\d+)_S(?P<s>\d+)_L001_R1_001\.fastq\.gz$")
first = None
for candidate in sorted(NEW_DIR.glob("*_L001_R1_001.fastq.gz")):
    if rx.match(candidate.name):
        first = candidate
        break
if first is None:
    raise RuntimeError(f"No valid L001 R1 FASTQ found in {NEW_DIR}")

m = rx.match(first.name)
idx = m.group("idx")
sample_s = m.group("s")

# Reproduce the full-run multi-lane input contract using four lightweight symlinks.
for old in OUT_DIR.glob("smoke_L*_R*.fastq.gz"):
    old.unlink()
for lane in ["001", "002", "003", "004"]:
    for read in ["1", "2"]:
        src = NEW_DIR / f"{idx}_S{sample_s}_L{lane}_R{read}_001.fastq.gz"
        if not src.exists():
            raise RuntimeError(f"Smoke-test lane is missing: {src}")
        dst = OUT_DIR / f"smoke_L{lane}_R{read}.fastq.gz"
        dst.symlink_to(src.resolve())

for old in OUT_DIR.glob("smoke*.clns"):
    old.unlink()
for old in OUT_DIR.glob("smoke*.vdjca"):
    old.unlink()

r1_pattern = str(OUT_DIR / "smoke_L{{n}}_R1.fastq.gz")
r2_pattern = str(OUT_DIR / "smoke_L{{n}}_R2.fastq.gz")
prefix = OUT_DIR / "smoke"

# Keep the command shape used by the historical presenting_repseq helper:
# executable, JVM memory, analyze, preset, -f, input pair, output prefix.
cmd = [
    mixcr,
    "-Xmx4g",
    "analyze",
    PRESET,
    "-f",
    "--threads", "2",
    "--limit-input", "10000",
    r1_pattern,
    r2_pattern,
    str(prefix),
]

version = subprocess.run([mixcr, "-v"], capture_output=True, text=True, check=True)
print("MiXCR binary:", mixcr)
print((version.stdout or version.stderr).splitlines()[0])
print("Smoke sequencing index:", idx)
print("Smoke R1 pattern:", r1_pattern)
print("Smoke R2 pattern:", r2_pattern)
print("Smoke preset:", PRESET)

with LOG.open("w", encoding="utf-8") as log:
    log.write("COMMAND: " + " ".join(cmd) + "\n")
    proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)

if proc.returncode != 0:
    tail = LOG.read_text(encoding="utf-8", errors="replace").splitlines()[-80:]
    print("\nMiXCR smoke test failed. Last log lines:\n" + "\n".join(tail))
    raise RuntimeError(
        f"MiXCR smoke test failed with exit code {proc.returncode}. "
        f"See {LOG}. Full cohort processing was not started."
    )

if not list(OUT_DIR.glob("smoke*.clns")):
    raise RuntimeError(
        f"MiXCR smoke test returned success but no .clns was produced in {OUT_DIR}. See {LOG}."
    )

print("MiXCR smoke test passed.")
