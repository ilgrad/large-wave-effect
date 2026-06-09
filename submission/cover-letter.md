Ilia Gradina
iliagradina@gmail.com

To the Editors,

I am pleased to submit the manuscript

  **Large-wave amplification on periodic discrete lattices: linear theory and a nonlinear extension**

for consideration. The paper concerns the *large-wave effect*: a localized impulse on a homogeneous
one-dimensional lattice (the ring of N beads) amplifies into a peak that grows with the system size — a
genuinely discrete resonance with no continuum counterpart, going back to the work of Filimonov, Myshkis
and Kurchanov.

**Main contribution.** Writing A_N = sup_t max_j |u_j(t)| for the unit-impulse velocity Green's function,
the paper:

1. **Settles the order for every N:** A_N ~ (1/π) ln N (precisely (1/π)ln N + O(ln ln ln N)). The lower
   bound rests on a short *palindromic* argument proving that the first φ(2N)/2 frequencies 2sin(πr/N) are
   Q-independent, so a (1−o(1)) fraction of the Bohr ceiling can always be Kronecker-aligned. This was the
   central open question for the effect.
2. **Classifies saturation:** A_N = U_N (the ceiling) if and only if N is prime or a power of two —
   equivalently, every integer frequency relation has coefficient-sum divisible by four — via an explicit
   root-of-unity obstruction valid for all composite N.
3. **Resolves the discrete Schrödinger amplification:** B_N = Θ(√N) in the ℓ∞ operator norm, with
   liminf B_N/√N = c₀/√2 and a *parity-dependent* limsup (the constant has no limit), the odd-subsequence
   value β_odd = 0.928… being an elliptic-integral constant with no elementary closed form.
4. **Extends and continues:** the order law to the Dirichlet segment and to every ring product C_N □ H; and
   a rigorous chain for the focusing discrete NLS (weak-nonlinear persistence, the modulational-instability
   band from the L_N spectrum, existence and stability of the on-site breather).

**Methodology and reproducibility.** The arguments combine almost-periodicity (Bohr/Kronecker–Weyl),
cyclotomic field theory, Toeplitz/circulant spectral analysis, and Bessel/Debye asymptotics. Every
quantitative claim is reproduced by a standalone, version-controlled script, and the number-theoretic core
is independently re-derived floating-point-free in GAP and PARI/GP, symbolically in FriCAS, and
machine-checked in Rocq. The complete repository is public.

**Scholarship.** I am careful to separate what is new from the prior Filimonov-group work: the segment
constant C = 4/π² is the published Myshkis–Filimonov (2003) value, and the qualitative ring effect for
prime/2^m is credited to that line. In the interest of full transparency I note one item I have not been
able to consult directly — Filimonov's 1992 C. R. Acad. Sci. note (its splash table is reproduced in the
2021 monograph of Andrianov, Awrejcewicz and Danishevskyy); the manuscript flags this and the novel content
(the arithmetic of *which* N, the free-ring setting, the order-theorem proof, the Schrödinger parity, and
the formal layer) accordingly.

**Open problems** are stated honestly and reduced to clean statements (a sharp B_N upper bound on t>N/2; a
rigorous growth lower bound for the deficit; the strongly nonlinear regime; condition (B) for general
quasi-1D lattices).

The work is original, not under consideration elsewhere, and has no competing interests. I would be glad to
suggest referees from the discrete spectral theory / almost-periodic functions / mathematical physics
communities upon request. Thank you for your consideration.

Sincerely,
Ilia Gradina
