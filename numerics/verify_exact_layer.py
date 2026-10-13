# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Drive the exact-arithmetic / formal layer and cross-check it against the Python core.

Runs the independent computer-algebra implementations of the linear-theory arithmetic and the
machine-checked criterion certificate, then verifies they agree with src/large_wave_effect:

  * GAP   (exact/gap/scan_classification.g)  -- rank = 1/2 phi(2N) and the saturation
          classification A_N = U_N <=> N prime or 2^m, in exact cyclotomic arithmetic.
  * GAP   (exact/gap/schrodinger_saturation.g) -- rank_Q{lambda_r} = 1/2 phi(N) and the prime-power
          relation lattice of the Schrodinger frequencies (second implementation of the PARI check).
  * PARI  (exact/pari/crosscheck.gp)         -- a second exact implementation of the Q-rank
          and the cyclotomic-degree identities behind the prime / 2^m / 2p proofs.
  * PARI  (exact/pari/bn_prime_flat.gp)      -- the exact inputs of B_p = sqrt(p): Q-independence
          of the Schrodinger eigenvalues and the flat Gauss-sum kernel, p <= 53.
  * PARI  (exact/pari/bn_prime_square.gp)    -- B_N = sqrt(N) for N = 4, 9, 25, 49: flat points of the
          orbit closure that satisfy every eigenvalue relation, checked exactly.
  * PARI  (exact/pari/defect_sos.gp)         -- exact SOS certificates of the defect U_N - A_N for
          N = 9, 15, 21, with mutated certificates that must be rejected.
  * PARI  (exact/pari/defect_lower_bound.gp) -- Proposition deficit-lower: the class relations of every
          odd prime p | N, checked in Q(zeta_2N), and the rational lower bounds for all
          non-saturating N <= 240 (triplets for 3 | N, the N = 2p sandwich).
  * PARI  (exact/pari/kjump_relations.gp)    -- exact Q-relation lattices of the k-jump circulant
          frequencies for N <= 126 (relations of C_N({1,3}) at N = 28, 56, 84, 112) and their mod-4 test.
  * PARI  (exact/pari/kjump_nnn_primes.gp)   -- the exact identities behind A_p = U_p for C_p({1,2}),
          every prime 7 <= p <= 211, and a Legendre-character certificate of the conclusion.
  * PARI  (exact/pari/bn_saturation.gp)      -- the Schrodinger saturation B_N = sqrt(N): rank of the
          lambda_r, the Gauss chirp and the even Frank point as exact flat points, the Golay search.
  * PARI  (exact/pari/support3_relations.gp) -- Lemma cleancomb: every support-3 relation
          omega_a + omega_b = omega_c for 3 | N <= 330, each decided modulo Phi_2N.
  * PARI  (exact/pari/kjump_kummer.gp)       -- Kummer certificate: saturation of C_p({1,2}) for primes
          7 <= p <= 97 and Q-independence of the NNN-ring frequencies at N = 15.
  * Julia (exact/julia/c15_flat_certificate.jl) -- Nullstellensatz certificate D y^4 = sum N_j m_j f_j
          for the N = 15 flatness system, by P-adic lifting, checked over Z: H_15 has no flat point.
  * PARI  (exact/pari/c15_flat_check.gp)     -- rebuilds the N = 15 flatness polynomials and re-checks
          that certificate by exact polynomial arithmetic.
  * Python (numerics/bn_flat_certificate_3q.py) -- equivariant Nullstellensatz certificates 1 = sum g_j A_j
          over the autocorrelations of the y-parametrized torus, N = 15 and 21, solved and checked over Q
          with python-flint: H_15 and H_21 have no flat point.
  * PARI  (exact/pari/c3q_flat_check.gp)     -- rebuilds the y-parametrization and the autocorrelations
          from their definitions and re-checks both certificates exactly in Z[L].
  * Rocq  (formal/rocq/Mod4Criterion.v)      -- machine-checks (vm_compute) that the mod-4
          criterion on the GAP lattices reproduces the prime/2^m verdict; Classification.v and
          StaggeredCurrent.v are compiled as well.

Tools are optional: a missing tool is reported as SKIP, not FAIL, so the script still runs the
pure-Python cross-checks. With GAP/PARI/Julia/Rocq installed (Fedora: dnf install gap pari-gp fricas
rocq-prover; Julia from juliaup) every layer is exercised.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def euler_phi(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def is_prime_or_2pow(n: int) -> bool:
    if n >= 2 and n & (n - 1) == 0:
        return True
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))


def run(
    cmd: list[str],
    stdin_file: Path | None = None,
    timeout: int = 300,
    env: dict[str, str] | None = None,
) -> tuple[int, str]:
    kw: dict = dict(cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    if env is not None:
        kw["env"] = {**os.environ, **env}
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
        gap_ok = (
            rc == 0 and "RESULT: PASS" in out
        )  # require clean exit, like the PARI/Rocq branches
        print(f"[gap]    classification scan N<=60          : {'PASS' if gap_ok else 'FAIL'}")
        ok &= gap_ok
        rc, out = run(["gap", "-q", "-A", "--quitonbreak", "exact/gap/schrodinger_saturation.g"])
        sgap_ok = rc == 0 and "RESULT: PASS" in out
        print(f"[gap]    Schrodinger rank + p^a lattice     : {'PASS' if sgap_ok else 'FAIL'}")
        ok &= sgap_ok
    else:
        print("[gap]    SKIP (gap not on PATH)")

    # --- PARI/GP: independent exact cross-check ---
    if shutil.which("gp"):
        rc, out = run(["gp", "-q", "exact/pari/crosscheck.gp"])
        pari_ok = rc == 0 and "FAIL" not in out and out.count("PASS") >= 5
        print(f"[pari]   cyclotomic cross-check             : {'PASS' if pari_ok else 'FAIL'}")
        ok &= pari_ok
        rc, out = run(
            ["gp", "-q", "-f", "exact/pari/bn_prime_flat.gp"], stdin_file=Path("/dev/null")
        )
        bp_ok = rc == 0 and "***" not in out and "RESULT: PASS" in out
        print(f"[pari]   B_p = sqrt(p) exact inputs         : {'PASS' if bp_ok else 'FAIL'}")
        ok &= bp_ok
        rc, out = run(
            ["gp", "-q", "-f", "exact/pari/bn_prime_square.gp"], stdin_file=Path("/dev/null")
        )
        bq_ok = rc == 0 and "***" not in out and "RESULT: PASS" in out
        print(f"[pari]   B_(p^2) = p for p <= 7             : {'PASS' if bq_ok else 'FAIL'}")
        ok &= bq_ok
        rc, out = run(["gp", "-q", "-f", "exact/pari/defect_sos.gp"], stdin_file=Path("/dev/null"))
        ds_ok = rc == 0 and "***" not in out and "RESULT: PASS" in out
        print(f"[pari]   defect SOS certificates N=9,15,21  : {'PASS' if ds_ok else 'FAIL'}")
        ok &= ds_ok
        rc, out = run(
            ["gp", "-q", "-f", "exact/pari/defect_lower_bound.gp"], stdin_file=Path("/dev/null")
        )
        dl_ok = rc == 0 and "***" not in out and "RESULT: PASS" in out
        print(f"[pari]   deficit lower bounds, N<=240       : {'PASS' if dl_ok else 'FAIL'}")
        ok &= dl_ok
        rc, out = run(
            ["gp", "-q", "-f", "exact/pari/kjump_relations.gp"], stdin_file=Path("/dev/null")
        )
        kr_ok = rc == 0 and "***" not in out and "RESULT: PASS" in out
        print(f"[pari]   k-jump relation lattices N<=126    : {'PASS' if kr_ok else 'FAIL'}")
        ok &= kr_ok
        rc, out = run(
            ["gp", "-q", "-f", "exact/pari/kjump_nnn_primes.gp"], stdin_file=Path("/dev/null")
        )
        kp_ok = rc == 0 and "***" not in out and "RESULT: PASS" in out
        print(f"[pari]   C_p({{1,2}}) saturation p<=211       : {'PASS' if kp_ok else 'FAIL'}")
        ok &= kp_ok
        # gp exits 0 after an error, so the *** marker is checked as well as the verdict line
        rc, out = run(["gp", "-q", "-f", "exact/pari/bn_saturation.gp"], timeout=600)
        sat_ok = rc == 0 and "RESULT: PASS" in out and "***" not in out
        print(f"[pari]   Schrodinger saturation B_N         : {'PASS' if sat_ok else 'FAIL'}")
        ok &= sat_ok
        rc, out = run(["gp", "-q", "-f", "exact/pari/support3_relations.gp"], timeout=600)
        s3_ok = rc == 0 and "RESULT: PASS" in out and "***" not in out
        print(f"[pari]   support-3 relations, 3|N<=330      : {'PASS' if s3_ok else 'FAIL'}")
        ok &= s3_ok
        rc, out = run(["gp", "-q", "-f", "exact/pari/kjump_kummer.gp"], timeout=600)
        kj_ok = rc == 0 and "RESULT: PASS" in out and "***" not in out
        print(f"[pari]   Kummer certificate, C_p({{1,2}})    : {'PASS' if kj_ok else 'FAIL'}")
        ok &= kj_ok
    else:
        print("[pari]   SKIP (gp not on PATH)")

    # --- Julia: Nullstellensatz certificate for N = 15, re-checked in PARI ---
    if shutil.which("julia"):
        with tempfile.TemporaryDirectory() as tmp:
            cert = Path(tmp) / "c15_certificate.gp"
            jl = ["julia", "--startup-file=no", "exact/julia/c15_flat_certificate.jl", str(cert)]
            rc, out = run(jl, timeout=900)
            jl_ok = rc == 0 and "RESULT: PASS" in out
            print(f"[julia]  N=15 certificate, no flat point    : {'PASS' if jl_ok else 'FAIL'}")
            ok &= jl_ok
            if jl_ok and shutil.which("gp"):
                gp = ["gp", "-q", "-f", "exact/pari/c15_flat_check.gp"]
                rc, out = run(gp, timeout=600, env={"C15_CERT": str(cert)})
                chk_ok = rc == 0 and "RESULT: PASS" in out and "***" not in out
                print(
                    f"[pari]   N=15 certificate re-checked        : {'PASS' if chk_ok else 'FAIL'}"
                )
                ok &= chk_ok
    else:
        print("[julia]  SKIP (julia not on PATH)")

    # --- python-flint: N = 15 and N = 21 certificates in the y-parametrization, re-checked in PARI ---
    if shutil.which("uv"):
        with tempfile.TemporaryDirectory() as tmp:
            for q in (5, 7):
                cert = Path(tmp) / f"c{3 * q}_y_certificate.txt"
                py = ["uv", "run", "--script", "numerics/bn_flat_certificate_3q.py", str(q)]
                rc, out = run([*py, "--out", str(cert)], timeout=1200)
                py_ok = rc == 0 and "RESULT: PASS" in out
                print(
                    f"[flint]  N={3 * q} y-certificate, no flat point  : {'PASS' if py_ok else 'FAIL'}"
                )
                ok &= py_ok
                if py_ok and shutil.which("gp"):
                    gp = ["gp", "-q", "-f", "exact/pari/c3q_flat_check.gp"]
                    rc, out = run(
                        gp, stdin_file=Path("/dev/null"), timeout=600, env={"C3Q_CERT": str(cert)}
                    )
                    chk_ok = rc == 0 and "RESULT: PASS" in out and "***" not in out
                    print(
                        f"[pari]   N={3 * q} y-certificate re-checked      : "
                        f"{'PASS' if chk_ok else 'FAIL'}"
                    )
                    ok &= chk_ok
    else:
        print("[flint]  SKIP (uv not on PATH)")

    # --- Rocq: machine-checked certificates ---
    if shutil.which("rocq"):
        for vname in (
            "formal/rocq/Mod4Criterion.v",
            "formal/rocq/Classification.v",
            "formal/rocq/StaggeredCurrent.v",
        ):
            if (ROOT / vname).exists():
                rc, out = run(["rocq", "compile", vname], timeout=300)
                rocq_ok = rc == 0
                print(f"[rocq]   compile {Path(vname).name:<22}: {'PASS' if rocq_ok else 'FAIL'}")
                ok &= rocq_ok
            else:
                print(f"[rocq]   SKIP {vname} (regenerate first)")
    else:
        print("[rocq]   SKIP (rocq not on PATH)")

    print("=" * 56)
    print("RESULT:", "ALL CHECKS PASS" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
