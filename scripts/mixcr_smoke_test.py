#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import os
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
    candidates.append(Path.home() / "soft" / "mixcr" / "mixcr")
    path_mixcr = shutil.which("mixcr")
    if path_mixcr:
        candidates.append(Path(path_mixcr))
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    raise RuntimeError(
        "MiXCR executable was not found. Expected the proven mouse-project binary at "
        "~/soft/mixcr/mixcr or set MEASLES_MIXCR_BIN explicitly."
    )


mixcr = resolve_mixcr()
r1_candidates = sorted(NEW_DIR.glob("*_L001_R1_001.fastq.gz"))
if not r1_candidates:
    raise RuntimeError(f"No L001 R1 FASTQ found in {NEW_DIR}")
r1 = r1_candidates[0]
r2 = Path(str(r1).replace("_R1_", "_R2_"))
if not r2.exists():
    raise RuntimeError(f"Paired R2 FASTQ is missing for smoke test: {r2}")

prefix = OUT_DIR / "smoke"
for old in OUT_DIR.glob("smoke*"):
    if old.is_file() or old.is_symlink():
        old.unlink()

# Follow the working presenting_repseq pattern from the mouse project:
# explicit MiXCR executable, JVM memory before `analyze`, preset, -f, R1, R2, prefix.
cmd = [
    mixcr,
    "-Xmx4g",
    "analyze",
    PRESET,
    "-f",
    "--threads", "2",
    "--limit-input", "10000",
    str(r1),
    str(r2),
    str(prefix),
]

version = subprocess.run([mixcr, "-v"], capture_output=True, text=True, check=True)
print("MiXCR binary:", mixcr)
print((version.stdout or version.stderr).splitlines()[0])
print("Smoke input R1:", r1)
print("Smoke input R2:", r2)
print("Smoke preset:", PRESET)

with LOG.open("w", encoding="utf-8") as log:
    log.write("COMMAND: " + " ".join(cmd) + "\n")
    proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)

if proc.returncode != 0:
    tail = LOG.read_text(encoding="utf-8", errors="replace").splitlines()[-60:]
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
