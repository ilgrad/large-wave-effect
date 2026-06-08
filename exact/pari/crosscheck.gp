\\ Independent PARI/GP cross-check of the cyclotomic facts used in the linear theory.
\\ A second implementation (distinct from GAP and from the Python core) of the Q-rank
\\ rank_Q {2 sin(pi r/N)} = 1/2 phi(2N), plus the cyclotomic-degree identities behind
\\ the prime, 2^m and 2p independence proofs.
\\
\\   gp -q exact/pari/crosscheck.gp

\\ Q-rank of {zeta^r - zeta^{-r} : r=1..N-1} in Q(zeta_{2N}), built from PARI polynomials
\\ (z^{-r} = z^{2N-r} since zeta^{2N} = 1).
freqrank(N) =
{
  my(c = polcyclo(2*N), d = poldegree(c), z = Mod(x, c), rows = List());
  for(r = 1, N-1, listput(rows, Vecrev(lift(z^r - z^(2*N - r)), d)));
  matrank(matrix(#rows, d, i, j, rows[i][j]));
}

ok = 1;
for(N = 3, 120, if(eulerphi(2*N)/2 != freqrank(N), ok = 0; print("  MISMATCH N=", N)));
print("rank_Q {2 sin(pi r/N)} = 1/2 phi(2N), 3<=N<=120 : ", if(ok, "PASS", "FAIL"));

t1 = 1; for(i = 1, 5, p = prime(i+1); if(polcyclo(2*p) != subst(polcyclo(p), x, -x), t1 = 0));
print("Phi_{2p}(x) = Phi_p(-x), odd p (Theorem prime)   : ", if(t1, "PASS", "FAIL"));

t2 = 1; for(m = 1, 6, if(polcyclo(2^(m+1)) != x^(2^m) + 1, t2 = 0));
print("Phi_{2N}(x) = x^N + 1, N = 2^m (Theorem two)     : ", if(t2, "PASS", "FAIL"));

t3 = 1; for(i = 1, 6, p = prime(i+1); if(subst(polcyclo(4*p), x, 1) != 1, t3 = 0));
print("Phi_{4p}(1) = Phi_p(-1) = 1, odd p (Prop 2p)     : ", if(t3, "PASS", "FAIL"));

t4 = 1; for(N = 3, 40, if(abs(prod(r = 1, N-1, 2*sin(Pi*r/N)) - N) > 1e-9, t4 = 0));
print("prod_{r} 2 sin(pi r/N) = N (det' L_N = N^2)      : ", if(t4, "PASS", "FAIL"));
quit
