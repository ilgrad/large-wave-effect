# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Self-trapping transition of the single-site DNLS seed: gamma*P = 4 is energy-sharp, not dynamically sharp.

For i z' = -L z + gamma|z|^2 z on the ring with z(0) = sqrt(P) e_0, the conserved energy is
    H = -<L z, z> + (gamma/2) ||z||_4^4 ,   H(0) = -2P + (gamma/2) P^2 .
Since <L z, z> >= 0, FULL dispersal (||z||_4^4 -> 0, the peak max_j|z_j|^2 -> 0) forces
    <L z, z> -> -H = 2P - (gamma/2) P^2 ,
which is >= 0 iff gamma*P <= 4.  Here 4 = 2<L e_0, e_0> is the infinite-lattice bandwidth
(||L_N|| = 4 cos^2(pi/2N) < 4 for odd N).  So:
  * gamma*P > 4: full dispersal is ENERGETICALLY FORBIDDEN; quantitatively max_j|z_j|^2 >= P - 4/gamma > 0
    (Proposition threshold).
  * gamma*P <= 4: full dispersal is energetically allowed -- but energy does not decide whether it happens.

This script checks, by split-step integration of the single-site seed (P=1):
  [1] N=512, T=200: H is conserved (drift < 1e-3); the bound inf_t max_j|z_j|^2 >= P - 4/gamma holds for
      every gamma*P > 4; and at the dispersal minimum <L z, z> = 2P - (gamma/2) P^2 to within 0.01 P;
  [2] N=4096, T=800 (front |j| <= 2T inside N/2, no wrap-around): the onset of trapping below 4.  The power
      within ten sites of the seed decays like the linear 1/t at gamma*P = 3.4, more slowly at 3.5, and
      settles near a third at 3.6 and near half at 3.7;
  [3] --long, N=16384, T=3200: at 3.4 the 1/t decay persists; at 3.5 a weak remnant (2-7% within ten
      sites) survives and rebounds.
So dispersal is established numerically only up to gamma*P ~ 3.4, while it is proved only for
gamma*P < 0.43 (Proposition dispersal); the long runs of Johansson-Hornquist-Riklund place the dynamical
threshold near 3.45.  Compare verify_dispersal_weak.py and verify_selftrapping_threshold.py.
"""

from __future__ import annotations

import sys

import numpy as np


def evolve(n: int, gamma: float, t_final: float, steps: int) -> tuple[float, float, float]:
    """Return (inf_t max|z|^2, |H drift|, <Lz,z> at the peak-minimizing time), P = 1."""
    lam = 2 - 2 * np.cos(2 * np.pi * np.arange(n) / n)  # ring Laplacian spectrum
    z = np.zeros(n, complex)
    z[0] = 1.0
    dt = t_final / steps
    half = np.exp(1j * lam * dt / 2)

    def energy(state: np.ndarray) -> float:
        sh = np.fft.fft(state)
        return -float(np.sum(lam * np.abs(sh) ** 2) / n) + 0.5 * gamma * float(
            np.sum(np.abs(state) ** 4)
        )

    def dirichlet(state: np.ndarray) -> float:
        sh = np.fft.fft(state)
        return float(np.sum(lam * np.abs(sh) ** 2) / n)

    h0 = energy(z)
    peak_min, luu_at_min = 1.0, dirichlet(z)
    for _ in range(steps):
        z = np.fft.ifft(half * np.fft.fft(z))
        z = np.exp(-1j * gamma * np.abs(z) ** 2 * dt) * z
        z = np.fft.ifft(half * np.fft.fft(z))
        peak = float((np.abs(z) ** 2).max())
        if peak < peak_min:
            peak_min, luu_at_min = peak, dirichlet(z)
    return peak_min, abs(energy(z) - h0), luu_at_min


def local_power(
    n: int, gamma: float, dt: float, times: tuple[float, ...], radius: int = 10
) -> list[float]:
    """Power within |j| <= radius of the seed at each of `times` (P = 1, Strang splitting)."""
    lam = 4 * np.sin(np.pi * np.arange(n) / n) ** 2
    e_full, e_half = np.exp(1j * lam * dt), np.exp(1j * lam * dt / 2)
    u0 = np.zeros(n, complex)
    u0[0] = 1.0
    w = np.fft.ifft(e_half * np.fft.fft(u0))  # half a linear step ahead of u
    marks = {round(t / dt): t for t in times}
    out = []
    for k in range(1, max(marks) + 1):
        w = w * np.exp(-1j * gamma * dt * (w.real**2 + w.imag**2))
        wh = np.fft.fft(w)
        w = np.fft.ifft(e_full * wh)
        if k in marks:
            p2 = np.abs(np.fft.ifft(e_half * wh)) ** 2
            out.append(float(p2[: radius + 1].sum() + p2[-radius:].sum()))
    return out


def main() -> int:
    long_run = "--long" in sys.argv[1:]
    ok = True

    n, t_final, steps = 512, 200.0, 30000
    print(f"[1] energy bound and energy identity (P=1, N={n}, T={t_final}); gamma*P = gamma.")
    print(
        f"{'gP':>5} {'inf max|z|^2':>13} {'P-4/gP':>8} {'|dH|':>8} {'<Lz,z>_min':>11} {'2P-gP^2/2':>10}"
    )
    for g in (0.5, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0):
        peak_min, hdrift, luu = evolve(n, g, t_final, steps)
        bound = max(0.0, 1 - 4 / g)
        energy_pred = 2 - g / 2  # 2P - (gamma/2)P^2 at full dispersal, P=1
        note = "disperses" if peak_min < 0.05 else ("trapped" if peak_min > 0.3 else "partial")
        good = peak_min >= bound - 0.02 and hdrift < 1e-3
        if note == "disperses":
            good &= abs(luu - energy_pred) < 0.01
        ok &= good
        print(
            f"{g:>5.1f} {peak_min:>13.4f} {bound:>8.3f} {hdrift:>8.1e} {luu:>11.4f} "
            f"{energy_pred:>10.4f}  {note}{'' if good else '  FAIL'}"
        )

    times = (100.0, 200.0, 400.0, 800.0)
    print("[2] onset below 4: power within 10 sites of the seed (N=4096, dt=0.005)")
    print(f"{'gP':>5} " + " ".join(f"{'t=' + str(int(t)):>8}" for t in times))
    p10 = {}
    for g in (3.4, 3.5, 3.6, 3.7):
        p10[g] = local_power(4096, g, 0.005, times)
        print(f"{g:>5.2f} " + " ".join(f"{v:>8.4f}" for v in p10[g]))
    rate = [v * t for v, t in zip(p10[3.4], times, strict=True)]
    checks = {
        "3.4 decays like 1/t": max(rate) < 1.05 * min(rate),
        "3.5 decays more slowly, P10(800) >= 0.02": p10[3.5][-1]
        >= max(0.02, 1.5 * p10[3.5][0] / 8),
        "3.6 holds about a third": all(0.30 <= v <= 0.36 for v in p10[3.6]),
        "3.7 holds nearly half": all(0.43 <= v <= 0.49 for v in p10[3.7]),
    }
    for name, good in checks.items():
        ok &= good
        print(f"    {name}: {'ok' if good else 'FAIL'}")

    if long_run:
        times = tuple(400.0 * k for k in range(1, 9))
        print("[3] long runs (N=16384, dt=0.005, T=3200)")
        print(f"{'gP':>5} " + " ".join(f"{'t=' + str(int(t)):>8}" for t in times))
        p34 = local_power(16384, 3.4, 0.005, times)
        p35 = local_power(16384, 3.5, 0.005, times)
        for g, vals in ((3.4, p34), (3.5, p35)):
            print(f"{g:>5.2f} " + " ".join(f"{v:>8.4f}" for v in vals))
        rate = [v * t for v, t in zip(p34, times, strict=True)]
        late = p35[1:]  # t >= 800
        checks = {
            "3.4 keeps the 1/t decay to t=3200": max(rate) < 1.1 * min(rate),
            "3.5 keeps 2-7% within ten sites for t in [800, 3200]": all(
                0.02 <= v <= 0.07 for v in late
            ),
            "3.5 rebounds after t=1600": max(p35[4:]) > 2 * p35[1],
        }
        for name, good in checks.items():
            ok &= good
            print(f"    {name}: {'ok' if good else 'FAIL'}")

    print("=" * 78)
    print(
        "RESULT:",
        "PASS -- full dispersal energy-forbidden above 4 (energy-sharp); numerically dispersal up to "
        "gamma*P ~ 3.4 and a localized remnant from ~3.5"
        if ok
        else "FAIL",
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
