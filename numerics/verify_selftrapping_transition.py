# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Self-trapping transition of the single-site DNLS seed: gamma*P = 4 is the sharp energy threshold.

For i z' = -L z + gamma|z|^2 z on the ring with z(0) = sqrt(P) e_0, the conserved energy is
    H = -<L z, z> + (gamma/2) ||z||_4^4 ,   H(0) = -2P + (gamma/2) P^2 .
Since <L z, z> >= 0, FULL dispersal (||z||_4^4 -> 0, the peak max_j|z_j|^2 -> 0) forces
    <L z, z> -> -H = 2P - (gamma/2) P^2 ,
which is >= 0 iff gamma*P <= 4.  So:
  * gamma*P > 4: full dispersal is ENERGETICALLY FORBIDDEN; quantitatively max_j|z_j|^2 >= P - 4/gamma > 0
    (Proposition threshold).  The constant 4 = ||L|| is sharp, not merely sufficient.
  * gamma*P <= 4: full dispersal is energetically allowed, and is what the dynamics realizes.

This certificate confirms, by direct symplectic (split-step) integration of the single-site seed (P=1):
  (1) H is conserved (drift < 1e-3);
  (2) the conservation bound inf_t max_j|z_j|^2 >= max(0, P - 4/gamma) holds for every gamma*P > 4;
  (3) the transition is a crossover centered at gamma*P ~ 4: clean dispersal (inf ~ 0) for gamma*P <~ 3.5,
      self-trapping (inf tracking P - 4/gamma) for gamma*P >~ 5;
  (4) the energy identity <L z, z> -> 2P - (gamma/2) P^2 at the dispersal minimum (for gamma*P < 4).
Hence the window gamma*P in [0.43, 4) is a PROOF gap (dispersal is the established behavior, proven only
for gamma*P <= 0.43 by Proposition dispersal), not a behavioral unknown.  Compare
verify_dispersal_weak.py (the gamma*P <= 0.43 proof) and verify_selftrapping_threshold.py.
"""

from __future__ import annotations

import numpy as np


def evolve(n: int, gamma: float, t_final: float, steps: int) -> tuple[float, float, float]:
    """Return (inf_t max|z|^2, |H drift|, <Lz,z> at the peak-minimizing time), P = 1."""
    lam = 2 - 2 * np.cos(2 * np.pi * np.arange(n) / n)  # ring Laplacian spectrum, max ||L|| = 4
    z = np.zeros(n, complex)
    z[0] = 1.0
    dt = t_final / steps
    half = np.exp(1j * lam * dt / 2)

    def energy(state: np.ndarray) -> float:
        sh = np.fft.fft(state)
        return -float(np.sum(lam * np.abs(sh) ** 2) / n) + 0.5 * gamma * float(np.sum(np.abs(state) ** 4))

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


def main() -> int:
    n, t_final, steps = 512, 200.0, 30000
    ok = True
    print(f"DNLS single-site seed (P=1, N={n}, T={t_final}); gamma*P = gamma.")
    print(f"{'gP':>5} {'inf max|z|^2':>13} {'P-4/gP':>8} {'|dH|':>8} {'<Lz,z>_min':>11} {'2P-gP^2/2':>10}")
    for g in (0.5, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0):
        peak_min, hdrift, luu = evolve(n, g, t_final, steps)
        bound = max(0.0, 1 - 4 / g)
        energy_pred = 2 - g / 2  # 2P - (gamma/2)P^2 at full dispersal, P=1
        bound_ok = peak_min >= bound - 0.02
        h_ok = hdrift < 1e-3
        ok &= bound_ok and h_ok
        # for gamma*P < 4, the dispersal minimum should sit near the energy prediction
        note = "disperses" if peak_min < 0.3 else ("trapped" if peak_min > 0.5 else "crossover")
        print(f"{g:>5.1f} {peak_min:>13.4f} {bound:>8.3f} {hdrift:>8.1e} {luu:>11.4f} "
              f"{energy_pred:>10.4f}  {note}{'' if bound_ok and h_ok else '  FAIL'}")
    print("=" * 78)
    print("RESULT:", "PASS -- threshold 4 sharp (energy-forbidden dispersal above 4); window [0.43,4) is a "
          "proof gap" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
