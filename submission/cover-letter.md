Ilia Gradina
iliagradina@gmail.com
ORCID: https://orcid.org/0009-0008-0911-8765

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
   root-of-unity obstruction valid for every composite N that is not a power of two. The obstruction is
   quantitative: U_N − A_N ≥ 1/(16p²) for p the least odd prime factor, above 0.043 whenever 3 | N, and
   asymptotic to π³/(16N²) for N = 2p.
3. **Locates the infinite-time saturation of the discrete Schrödinger amplification:** B_N = Θ(√N)
   in the ℓ∞ operator norm, with liminf B_N/√N ≥ c₀/√2 proved; B_N = √N exactly when the orbit closure (of
   dimension φ(N)/2) contains a flat phase vector, which holds for N = p and N = p² (the ε-uniform mixing
   of Godsil–Mullin–Roy and Cao–Zhu) and fails for every even N ≥ 6 (new for the even squares; for the
   other even N it also follows from Monterde's bipartite conditions) and for N = 15, 21 and 27 (exact
   certificates, computer-assisted: Nullstellensatz identities for 15 and 21, a Positivstellensatz on the
   9-torus for 27; for C_21 and C_27 this is new even for exact uniform mixing: Gray–Patel–Young name C_21 as
   the next open case, and Cao–Zhu leave C_27 open), so limsup B_N/√N = 1. On the
   dispersive time scale t = O(N) the constant splits by parity; the odd-N value β_odd = 0.928… is an elliptic-integral constant with no elementary closed form, and its
   windowed sharp upper bound is reduced to one residual inequality (proved by computer-assisted Arb
   certificates for every configuration with at most three copies, validated beyond) and an s-uniform transfer.
4. **Extends and continues:** the order law to the Dirichlet segment and to every ring product C_N □ H; and
   a rigorous chain for the focusing discrete NLS (weak-nonlinear persistence, the modulational-instability
   band from the L_N spectrum, no dispersal of a single-site excitation above γP = 4 and on-site trapping
   above γP = 8, existence of breathers and orbital stability of the site-centred one in the anti-continuum
   regime, with computer-assisted verdicts at one frequency).

**Methodology and reproducibility.** The arguments combine almost-periodicity (Bohr/Kronecker–Weyl),
cyclotomic field theory, Toeplitz/circulant spectral analysis, and Bessel/Debye asymptotics. Every
quantitative claim is reproduced by a standalone, version-controlled script, and the number-theoretic core
is independently re-derived floating-point-free in GAP and PARI/GP, symbolically in FriCAS, and
machine-checked in Rocq. The complete repository is public.

**Scholarship.** I am careful to separate what is new from the prior Filimonov-group work: the segment
constant C = 4/π² is the published Myshkis–Filimonov (2003) value, and the qualitative ring effect for
prime/2^m is credited to that line. Filimonov's 1992 C. R. Acad. Sci. note, compared in detail in the
manuscript, already contains the Kronecker identification of the ring ceiling under Q-independence, the
prime/2^m independence by the same cyclotomic argument (our proof is a streamlined reproof), the
even-composite dependence, and the unbounded growth. The new content is the rate (1/π) ln N for every N,
the saturation classification with its mod-4 criterion, the closed-form rank, the Schrödinger results, and
the formal layer.

**Open problems** are stated honestly and reduced to clean statements (the remaining odd composite N in the B_N
saturation pattern, and the sharp windowed B_N bound on N/2 < t ≤ SN/2; a
rigorous growth lower bound for the deficit; the strongly nonlinear regime; condition (B) for general
quasi-1D lattices).

The work is original, not under consideration elsewhere, and has no competing interests. I would be glad to
suggest referees from the discrete spectral theory / almost-periodic functions / mathematical physics
communities upon request. Thank you for your consideration.

Sincerely,
Ilia Gradina
