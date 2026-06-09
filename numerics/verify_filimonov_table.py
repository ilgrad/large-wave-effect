# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Consistency check against Filimonov's published splash table (prior art).

The 1992 C. R. Acad. Sci. note [F92] established the "splash" amplitudes P_N for the string with N beads;
they are reproduced in the monograph [AAD21] (Andrianov-Awrejcewicz-Danishevskyy, CRC Press 2021), Table
10.1.  We check the prior-art boundary used in the paper's "Relation to prior work":

  * that published table grows at the SEGMENT constant 4/pi^2 (the Myshkis-Filimonov [MF03] log constant
    for the Dirichlet string), NOT at our RING constant 1/pi;
  * our ring ceiling U_N = (1/2N) sum_{r=1}^{N-1} csc(pi r/N) grows at 1/pi.

So the documented F92 result is the string/segment splash magnitude (P_N -> infinity at rate (4/pi^2)ln N),
while the ring law A_N ~ (1/pi) ln N and the arithmetic of WHICH N saturate are the present contribution.
This script does NOT reproduce F92's exact P_N (a stress, with F92's own initial data); it certifies the
growth CONSTANT that the prior-art claim turns on.
"""

from __future__ import annotations

import numpy as np
from numpy import log, pi, sin

# [AAD21] Table 10.1 -- Filimonov's splash amplitudes P_N for the string with N beads.
BOOK_TABLE = {8: 1.7561, 16: 2.0645, 32: 2.3468, 64: 2.6271, 128: 2.9078, 256: 3.1887}


def ring_ceiling(n: int) -> float:
    return float(sum(1.0 / sin(pi * r / n) for r in range(1, n)) / (2 * n))


def main() -> int:
    ns = np.array(sorted(BOOK_TABLE))
    ps = np.array([BOOK_TABLE[n] for n in ns])
    seg_slope = float(np.polyfit(log(ns), ps, 1)[0])

    rc = np.array([ring_ceiling(int(n)) for n in ns])
    ring_slope = float(np.polyfit(log(ns), rc, 1)[0])

    c_seg, c_ring = 4 / pi**2, 1 / pi
    print(f"[F92/AAD21 Table 10.1]  splash P_N slope d/d(lnN) = {seg_slope:.4f}   (4/pi^2 = {c_seg:.4f})")
    print(f"[our ring]   U_N slope d/d(lnN)              = {ring_slope:.4f}   (1/pi   = {c_ring:.4f})")
    print(f"  segment/ring constant ratio = {seg_slope / ring_slope:.3f}  (4/pi^2 : 1/pi = {c_seg / c_ring:.3f})")

    seg_ok = abs(seg_slope - c_seg) < 0.02            # table grows at the segment constant
    ring_ok = abs(ring_slope - c_ring) < 0.01         # ring ceiling grows at 1/pi
    distinct = seg_slope - ring_slope > 0.05          # the two regimes are genuinely different
    ok = seg_ok and ring_ok and distinct

    print("=" * 70)
    print("RESULT:", "PASS -- F92 table is the segment 4/pi^2 regime; ring is 1/pi" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
