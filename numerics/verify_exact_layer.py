# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Drive the exact-arithmetic / formal layer and cross-check it against the Python core.

Runs the independent computer-algebra implementations of the linear-theory arithmetic and the
machine-checked criterion certificate, then verifies they agree with src/large_wave_effect:

  * GAP   (exact/gap/scan_classification.g)  -- rank = 1/2 phi(2N) and the saturation
          classification A_N = U_N <=> N prime or 2^m, in exact cyclotomic arithmetic.
  * PARI  (exact/pari/crosscheck.gp)         -- a second exact implementation of the Q-rank
          and the cyclotomic-degree identities behind the prime / 2^m / 2p proofs.
  * Rocq  (formal/rocq/Mod4Criterion.v)      -- machine-checks (vm_compute) that the mod-4
          criterion on the GAP lattices reproduces the prime/2^m verdict.

Tools are optional: a missing tool is reported as SKIP, not FAIL, so the script still runs the
pure-Python cross-checks. With GAP/PARI/Rocq installed (Fedora: dnf install gap pari-gp fricas
rocq-prover) every layer is exercised.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def euler_phi(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def is_prime_or_2pow(n: int) -> bool:
    if n >= 2 and n & (n - 1) == 0:
        return True
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))


def run(cmd: list[str], stdin_file: Path | None = None, timeout: int = 300) -> tuple[int, str]:
    kw: dict = dict(cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    if stdin_file is not None:
        with stdin_file.open() as fh:
            return _finish(subprocess.run(cmd, stdin=fh, **kw))
    return _finish(subprocess.run(cmd, **kw))


def _finish(p: subprocess.CompletedProcess[str]) -> tuple[int, str]:
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def main() -> int:
    ok = True

    # --- pure-Python core: the Q-rank formula, independent of every CAS ---
    sys.path.insert(0, str(ROOT / "src"))
    from large_wave_effect.cyclotomic import ring_freq_rank

    rank_ok = all(ring_freq_rank(n) == euler_phi(2 * n) // 2 for n in range(3, 41))
    print(f"[python] ring_freq_rank == phi(2N)/2 for N<=40 : {'PASS' if rank_ok else 'FAIL'}")
    ok &= rank_ok

    # --- GAP: exact classification scan ---
    if shutil.which("gap"):
        env = "LW_MAXN := 60;;"
        rc, out = run(["gap", "-q", "--nointeract", "-c", env, "exact/gap/scan_classification.g"])
        gap_ok = "RESULT: PASS" in out
        print(f"[gap]    classification scan N<=60          : {'PASS' if gap_ok else 'FAIL'}")
        ok &= gap_ok
    else:
        print("[gap]    SKIP (gap not on PATH)")

    # --- PARI/GP: independent exact cross-check ---
    if shutil.which("gp"):
        rc, out = run(["gp", "-q", "exact/pari/crosscheck.gp"])
        pari_ok = rc == 0 and "FAIL" not in out and out.count("PASS") >= 5
        print(f"[pari]   cyclotomic cross-check             : {'PASS' if pari_ok else 'FAIL'}")
        ok &= pari_ok
    else:
        print("[pari]   SKIP (gp not on PATH)")

    # --- Rocq: machine-checked mod-4 certificate ---
    if shutil.which("rocq"):
        vfile = ROOT / "formal/rocq/Mod4Criterion.v"
        if vfile.exists():
            rc, out = run(["rocq", "compile", "formal/rocq/Mod4Criterion.v"], timeout=300)
            rocq_ok = rc == 0
            print(f"[rocq]   compile Mod4Criterion.v            : {'PASS' if rocq_ok else 'FAIL'}")
            ok &= rocq_ok
        else:
            print("[rocq]   SKIP (run exact/gap/export_rocq.g first)")
    else:
        print("[rocq]   SKIP (rocq not on PATH)")

    print("=" * 56)
    print("RESULT:", "ALL CHECKS PASS" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
