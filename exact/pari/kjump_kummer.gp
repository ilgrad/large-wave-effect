\\ Exact (Kummer/Besicovitch) proof of saturation A_N = U_N for the next-nearest-neighbour circulant
\\ C_p({1,2}) at collision-free primes p.  Run:  gp -q exact/pari/kjump_kummer.gp
\\
\\ The wave frequencies are omega_r = sqrt(alpha_r), r = 1..(p-1)/2, with
\\     alpha_r = 4 - 2 cos(2 pi r/p) - 2 cos(4 pi r/p)
\\             = 4 - (z^r + z^-r) - (z^2r + z^-2r),   z = zeta_p,
\\ a totally positive element of Q(zeta_p).  Saturation A_N = U_N holds iff all (distinct) frequencies are
\\ Q-linearly independent (then Kronecker aligns every sin(omega_r t) -> 1 simultaneously).  By the
\\ Besicovitch criterion, sqrt(alpha_1),...,sqrt(alpha_k) are linearly independent over K = Q(zeta_p) -- hence
\\ over Q -- iff NO nonempty subset product prod_{r in T} alpha_r is a square in K.  This is an EXACT test
\\ (nfeltispower over the cyclotomic field), strictly stronger than the floating-point PSLQ evidence.

checkp(p) = {
  my(nf = nfinit(polcyclo(p)), z = Mod(x, polcyclo(p)), k = (p-1)/2, anysq = 0);
  my(alpha = vector(k, r, Mod(4 - (z^r + z^(-r)) - (z^(2*r) + z^(-2*r)), polcyclo(p))));
  forvec(e = vector(k, i, [0,1]),
    if(vecsum(e) == 0, next);
    my(P = Mod(1, polcyclo(p)));
    for(i = 1, k, if(e[i], P = P * alpha[i]));
    if(nfeltispower(nf, lift(P), 2), anysq = 1; break));
  printf("  C_%d({1,2}): %d distinct freqs, all 2^%d-1 subset products non-square -> %s\n",
    p, k, k, if(anysq, "FAIL (relation exists)", "PASS (Q-independent, saturation proven)"));
  !anysq;
}

allok = 1;
\\ collision-free primes (S and -S do not cover all nonzero residues mod p, i.e. p > 5)
for(j = 1, 6, allok = allok && checkp([7, 11, 13, 17, 19, 23][j]));
print("RESULT: " if(allok, "PASS", "FAIL") " -- saturation proven exactly for NNN circulant, primes 7..23");
quit;
