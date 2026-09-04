#!/usr/bin/env python3
"""Run analysis_08_final_reporting.py with a runtime compatibility fix.

Pandas/NumPy can return a read-only ndarray from Series.to_numpy().  The
original Morisita-Horn helper used in-place normalization (a /= a.sum()),
which raises ``ValueError: output array is read-only`` on the current server
environment.  This wrapper applies the minimal safe source transformation in
memory and then executes the authoritative script with ``__file__`` set to the
original script path.

No input files are modified and no MiXCR step is rerun.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "scripts" / "analysis_08_final_reporting.py"

source = TARGET.read_text(encoding="utf-8")
old = "a/=a.sum();b/=b.sum();den=(a*a).sum()+(b*b).sum();return 2*(a*b).sum()/den if den else np.nan"
new = "a=a/a.sum();b=b/b.sum();den=(a*a).sum()+(b*b).sum();return 2*(a*b).sum()/den if den else np.nan"

if old not in source:
    raise RuntimeError(
        "Expected Morisita-Horn normalization code was not found in "
        f"{TARGET}. Refusing to run an unverified transformation."
    )

source = source.replace(old, new, 1)
namespace = {
    "__file__": str(TARGET),
    "__name__": "__main__",
    "__package__": None,
}
exec(compile(source, str(TARGET), "exec"), namespace)
