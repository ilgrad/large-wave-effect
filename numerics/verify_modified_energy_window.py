# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Modified-energy probe for the DNLS self-trapping window gamma*P in [0.43, 4) (paper, eq:stagger).

The staggered momentum S = Re sum_n conj(b_n) b_{n+2} = P<cos 2k> drifts only through the nonlinearity,
    dS/dt = mu * F[b],   F[b] = -sum_n (|b_n|^2 - |b_{n+2}|^2) Im(conj(b_n) b_{n+2}),   mu = gamma*P.
Lemma lem:stagger-leading (verify_staggered_current.{py,mac}, formal/rocq/StaggeredCurrent.v) proves
F[b_lin] = 0 on the LINEAR flow, so the O(mu) drift cancels and the first nonzero contribution is O(mu^2).
This probe confirms the consequence on the full NONLINEAR flow and quantifies the route:
  (1) mu-scaling of sup_t|Delta S| has exponent p ~ 2 (not 1)  => leading secular drift cancels;
  (2) the linear-flow secular rate Fbar ~ 0  (independent check of the lemma);
  (3) sup_t|Delta S| is N-uniform (drift is a bounded transient, not secular accumulation);
  (4) inf_t <v^2> = 2 - 2 S/P stays > 0 across the window (the dispersive virial margin).

DNLS i b' = -(L b) + mu|b|^2 b, single-site seed ||b||^2 = 1.  Norm-conserving Strang split-step,
batched over mu (one batched FFT sweep).  Runs on CPU by default; for the GPU path use
  uv run --with cupy-cuda12x --with nvidia-cufft-cu12 --with nvidia-nvjitlink-cu12 \
         python numerics/verify_modified_energy_window.py
"""
from __future__ import annotations

import time
from types import ModuleType

import numpy as np

try:
    # cupy-cuda12x finds elementwise kernels via the driver, but cuFFT is a separate shared lib; on a
    # CUDA-13 host the CUDA-12 cuFFT (soname .11) lives in the nvidia-*-cu12 wheels and is not on the loader
    # path, so preload it (with its nvJitLink dependency) via RTLD_GLOBAL before importing cupy.
    import ctypes as _ct
    import os as _os

    import nvidia.cufft as _cf
    import nvidia.nvjitlink as _nj

    for _m, _so in [(_nj, "libnvJitLink.so.12"), (_cf, "libcufft.so.11")]:
        _ct.CDLL(_os.path.join(next(iter(_m.__path__)), "lib", _so), mode=_ct.RTLD_GLOBAL)
    import cupy as _cp

    XP: ModuleType = _cp if _cp.cuda.runtime.getDeviceCount() > 0 else np
except Exception:
    XP = np

GPU = XP is not np
EPS = (0.2, 0.1, 0.05)
DT = 0.01
CHUNK = 2_000_000
NS = (1024, 2048, 4096) if GPU else (512, 1024, 2048)
CHECKPOINTS = (1e3, 1e4, 1e5)


def strang(b, phase_half, mu_col, dt, xp):
    b = xp.fft.ifft(xp.fft.fft(b, axis=-1) * phase_half, axis=-1)
    b = b * xp.exp(-1j * mu_col * (xp.abs(b) ** 2) * dt)
    return xp.fft.ifft(xp.fft.fft(b, axis=-1) * phase_half, axis=-1)


def run_sweep(n: int, mus: np.ndarray, t_final: float, xp: ModuleType) -> dict:
    m = len(mus)
    mu_col = xp.asarray(mus).reshape(m, 1)
    lam = 2.0 - 2.0 * xp.cos(2.0 * np.pi * xp.arange(n) / n)
    phase_half = xp.exp(1j * lam * (DT / 2))
    b = xp.zeros((m, n), dtype=xp.complex128)
    b[:, 0] = 1.0  # single-site seed, ||b||^2 = 1
    n_steps = round(t_final / DT)
    stride = max(1, n_steps // 2000)  # S varies on an O(1) t-scale; ~2000 samples resolve the running max
    maxS = np.zeros(m)
    minv2 = np.full(m, 2.0)
    for it in range(n_steps + 1):
        if it % stride == 0:
            p = xp.sum(xp.abs(b) ** 2, axis=-1)
            s = xp.real(xp.sum(xp.conj(b) * xp.roll(b, -2, axis=-1), axis=-1))
            s_np = np.asarray(s.get() if GPU else s)
            v2_np = np.asarray((2.0 - 2.0 * s / p).get() if GPU else 2.0 - 2.0 * s / p)
            maxS = np.maximum(maxS, np.abs(s_np))
            minv2 = np.minimum(minv2, v2_np)
        if it < n_steps:
            b = strang(b, phase_half, mu_col, DT, xp)
    return {"maxS": maxS, "minv2": minv2}


def linear_secular(n: int, t_final: float, xp: ModuleType) -> float:
    """Time-average of F along the PURE LINEAR flow from the seed (the leading secular rate; lemma => 0)."""
    lam = 2.0 - 2.0 * xp.cos(2.0 * np.pi * xp.arange(n) / n)
    phase = xp.exp(1j * lam * DT)
    b = xp.zeros((1, n), dtype=xp.complex128)
    b[:, 0] = 1.0
    acc = 0.0
    n_steps = round(t_final / DT)
    for _ in range(n_steps):
        dens = xp.abs(b) ** 2
        cur = xp.imag(xp.conj(b) * xp.roll(b, -2, axis=-1))
        f = -xp.sum((dens - xp.roll(dens, -2, axis=-1)) * cur)
        acc += float(f) * DT
        b = xp.fft.ifft(xp.fft.fft(b, axis=-1) * phase, axis=-1)
    return acc / t_final


def main() -> int:
    print("=" * 80)
    print(f"DNLS modified-energy window probe   backend={'GPU' if GPU else 'CPU'}   N={NS}")
    print("=" * 80)
    mus = np.round(np.linspace(0.30, 4.40, 42), 3)
    win = (mus >= 0.43) & (mus <= 1.5)  # perturbative part of the window, for the mu-scaling fit
    winf = (mus >= 0.43) & (mus < 4.0)
    t_start = time.time()
    maxS_by_n, results = {}, []
    for n in NS:
        r = run_sweep(n, mus, n / 5.0, XP)
        fbar = linear_secular(n, min(n / 5.0, 300.0), XP)
        p_exp = float(np.polyfit(np.log(mus[win]), np.log(r["maxS"][win]), 1)[0])
        maxS_by_n[n] = r["maxS"]
        results.append((n, p_exp, fbar, r["maxS"][winf], r["minv2"][winf]))
        print(f"N={n:>5}: mu-scaling p={p_exp:.2f}  Fbar_lin={fbar:.2e}  "
              f"inf<v2>={r['minv2'][winf].min():.3f}  sup|dS|@1.5={r['maxS'][win][-1]:.3e}")
    # N-uniformity: spread of sup|dS| at a fixed mu across N
    j15 = int(np.argmin(abs(mus - 1.5)))
    vals = np.array([maxS_by_n[n][j15] for n in NS])
    n_unif = float(np.ptp(vals) / np.mean(vals))
    p_all = np.mean([row[1] for row in results])
    fbar_max = max(abs(row[2]) for row in results)
    v2_min = min(row[4].min() for row in results)
    print("-" * 80)
    print(f"mean mu-scaling exponent p = {p_all:.2f}  (p~2 => leading secular drift cancels; lemma)")
    print(f"max |Fbar_lin| = {fbar_max:.1e}  (~0 => no leading secular obstruction; lemma)")
    print(f"N-uniformity of sup|dS|@1.5 = {n_unif:.1e} (relative spread across N)")
    print(f"min <v^2> over window [0.43,4) = {v2_min:.3f}  (>0 => dispersive margin)")
    print(f"elapsed {time.time()-t_start:.1f}s")
    ok = (1.7 < p_all < 2.4) and (fbar_max < 1e-9) and (n_unif < 0.05) and (v2_min > 0.3)
    print("=" * 80)
    print("RESULT:", "MODIFIED-ENERGY ROUTE SUPPORTED (drift O(mu^2), N-uniform)" if ok else "CHECK FAILED")
    print("=" * 80)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
