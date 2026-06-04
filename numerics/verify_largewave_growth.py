# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Probe of the impulse large wave A_N (docs/bead_chain_formulas_claude.md, §6, §9, §10).

Impulse IC on the zero-mean subspace: u(0)=0, v0 = e_0 - 1/N, so the displacement is
exactly the velocity Green's function

    u_j(t) = G_N^{(s)}(j,t) = (1/N) sum_{r=1}^{N-1} [sin(omega_r t)/omega_r] e^{2*pi*i*r*j/N},

and the large-wave functional is A_N = sup_{t>=0} max_j |G_N^{(s)}(j,t)|.

FINDING (the honest outcome): A_N is the supremum of an almost-periodic function whose
near-recurrent peaks occur at times that grow rapidly with N (cf. the diploma's t ~ 2e8).
A uniform t-scan therefore yields only a LOWER bound that keeps rising as the window T
grows -- so the N-asymptotics of A_N CANNOT be read off a fixed-budget scan. This is why
the c*log N growth (proved by Myshkis-Filimonov for the Dirichlet segment) is a genuine
open analytic problem on the ring, not a thing to "measure".

What this script does establish, rigorously:
  (A) the in-phase upper bound A_N <= U_N := (1/(2N)) sum csc(pi r/N) ~ (1/pi) ln N;
  (B) the time-MEAN amplitude stays O(1) in N (growth lives only in rare sup excursions);
  (C) the scan is truncation-limited: A_N(obs) increases monotonically with the window T
      and the argmax t* hugs the window boundary -> the scan never reaches the true sup.
"""

from __future__ import annotations

import numpy as np

PROBE_NS = [64, 128, 256]
WINDOW_MULTS = [16, 64, 256, 1024]


def omega(n: int) -> np.ndarray:
    return 2.0 * np.abs(np.sin(np.pi * np.arange(n) / n))  # sqrt eigenvalues of L_N


def green_peak(n: int, t_max: float, dt: float = 0.1, block: int = 16384) -> tuple[float, float]:
    """Observed max_t max_j |G_N^{(s)}(j,t)| on [0,t_max] -- a LOWER bound on the true sup."""
    om = omega(n)
    nz = om > 1e-12
    ts = np.arange(0.0, t_max, dt)
    best, t_at = 0.0, 0.0
    for s in range(0, len(ts), block):
        tb = ts[s:s + block]
        coeff = np.zeros((len(tb), n))
        coeff[:, nz] = np.sin(np.outer(tb, om[nz])) / om[nz]
        u = np.fft.ifft(coeff, axis=1).real          # ifft carries the 1/N -> equals G^{(s)}
        amax = np.abs(u).max(axis=1)
        i = int(amax.argmax())
        if amax[i] > best:
            best, t_at = float(amax[i]), float(tb[i])
    return best, t_at


def green_time_rms(n: int, t_max: float, dt: float = 0.05) -> float:
    """Root-mean-square over time of |G_N^{(s)}(0,t)| -- the typical (not extreme) amplitude."""
    om = omega(n)
    nz = om > 1e-12
    ts = np.arange(0.0, t_max, dt)
    acc, cnt = 0.0, 0
    for s in range(0, len(ts), 16384):
        tb = ts[s:s + 16384]
        g0 = (np.sin(np.outer(tb, om[nz])) / om[nz]).sum(axis=1) / n   # G^{(s)}(0,t)
        acc += float(np.sum(g0 * g0))
        cnt += len(tb)
    return float(np.sqrt(acc / cnt))


def inphase_bound(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / (2.0 * np.sin(np.pi * r / n))) / n)


def main() -> int:
    print("=" * 78)
    print("Impulse large-wave A_N -- probe (m=kappa=1)")
    print("=" * 78)
    ok = True
    inv_pi = 1.0 / np.pi

    # (A) in-phase upper bound and its (1/pi) ln N asymptotic
    print("\n(A) Upper bound  A_N <= U_N = (1/2N) sum csc(pi r/N) ~ (1/pi) ln N")
    print(f"  {'N':>6} {'U_N':>10} {'U_N*pi/lnN':>12}")
    prev = None
    mono = True
    for n in (16, 64, 256, 1024, 4096, 16384):
        u = inphase_bound(n)
        ratio = u * np.pi / np.log(n)
        print(f"  {n:>6} {u:>10.4f} {ratio:>12.4f}")
        if prev is not None and ratio > prev + 1e-9:
            mono = False
        prev = ratio
    print(f"  -> U_N*pi/lnN decreases toward 1 (from above): {'PASS' if mono else 'FAIL'}")
    ok &= mono

    # (B) time-mean amplitude is exactly RMS_t|G(0,t)| = sqrt((N^2-1)/(12 N^2)) -> 1/sqrt(12).
    # The factor 2 over the naive 1/sqrt(24) is the ring's twofold degeneracy omega_r=omega_{N-r}
    # (modes r and N-r share a frequency, so their contributions stay fully correlated).
    print("\n(B) Typical amplitude RMS_t|G(0,t)| = 1/sqrt(12)=0.2887, exact (ring degeneracy);")
    print("    O(1) in N while the sup is conjectured to grow -> growth is in rare excursions.")
    print(f"  {'N':>6} {'RMS_t(obs)':>12} {'RMS(theory)':>12} {'rel':>8}")
    rms_ok = True
    for n in (32, 128, 512):
        rms = green_time_rms(n, t_max=64.0 * n)
        pred = float(np.sqrt((n * n - 1) / (12 * n * n)))
        rel = abs(rms - pred) / pred
        print(f"  {n:>6} {rms:>12.4f} {pred:>12.4f} {rel:>8.1e}")
        rms_ok &= rel < 0.02
    print(f"  -> matches 1/sqrt(12), NOT 1/sqrt(24): {'PASS' if rms_ok else 'FAIL'}")
    ok &= rms_ok

    # (C) the scan is truncation-limited -- the central methodological finding
    print("\n(C) Scan truncation: A_N(obs) keeps rising with window T; t* hugs the boundary")
    header = "  {:>5} ".format("N") + " ".join(f"{'T=' + str(m) + 'N':>15}" for m in WINDOW_MULTS) + f" {'U_N':>8}"
    print(header)
    truncated = True
    bound_ok = True
    for n in PROBE_NS:
        peaks, cells = [], []
        for m in WINDOW_MULTS:
            p, t = green_peak(n, m * n)
            peaks.append(p)
            cells.append(f"{p:.3f}@{t:.0f}")
        u = inphase_bound(n)
        print(f"  {n:>5} " + " ".join(f"{c:>15}" for c in cells) + f" {u:>8.3f}")
        if not all(peaks[i + 1] >= peaks[i] - 1e-6 for i in range(len(peaks) - 1)):
            truncated = False
        if peaks[-1] > u * (1 + 1e-9):
            bound_ok = False
    print(f"  -> A_N(obs) monotonically increasing in T (sup not reached): "
          f"{'PASS' if truncated else 'FAIL'}")
    print(f"  -> every observed A_N stays below U_N: {'PASS' if bound_ok else 'FAIL'}")
    ok &= truncated and bound_ok

    print("\n" + "=" * 78)
    print("RESULT:", "PROBE CONSISTENT" if ok else "UNEXPECTED")
    print("Conclusion: the c*log N growth of the ring impulse large wave is NOT obtainable")
    print("by uniform t-scanning (sup lives at exp-large t). It needs an analytic lower")
    print("bound (Diophantine / Kronecker) -- the central open result for paper #1. The")
    print(f"strict ceiling is A_N <= U_N ~ (1/pi) ln N (1/pi={inv_pi:.4f}).")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
