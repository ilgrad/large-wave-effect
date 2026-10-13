\\ Exact (Kummer/Besicovitch) proof of saturation A_N = U_N for the next-nearest-neighbour circulant
\\ C_p({1,2}) at collision-free primes p, and for the NNN ring (g = 1/2) at N = 15.
\\ Run:  gp -q -f exact/pari/kjump_kummer.gp </dev/null
\\
\\ C_p({1,2}): the wave frequencies are omega_r = sqrt(alpha_r), r = 1..(p-1)/2, with
\\     alpha_r = 4 - 2 cos(2 pi r/p) - 2 cos(4 pi r/p)
\\             = 4 - (z^r + z^-r) - (z^2r + z^-2r),   z = zeta_p,
\\ a totally positive algebraic integer of Q(zeta_p).  Saturation A_N = U_N holds if all (distinct)
\\ frequencies are Q-linearly independent (then Kronecker aligns every sin(omega_r t) -> 1 simultaneously).
\\ By the Besicovitch criterion, sqrt(alpha_1),...,sqrt(alpha_k) are linearly independent over K = Q(zeta_p)
\\ -- hence over Q -- if NO nonempty subset product prod_{r in T} alpha_r is a square in K.
\\
\\ Certificate (quadratic characters).  Let q = 1 (mod m) be prime and w an element of order m in F_q; each
\\ z = w^j, gcd(j, m) = 1, is the image of zeta_m at a prime ideal of residue degree one above q.  If
\\ prod_{r in T} alpha_r = beta^2 with beta in K, then beta is integral, so prod_{r in T} alpha_r(z) is a square
\\ mod q.  Hence the indicator vector of T lies in the kernel of the F_2-matrix whose rows are
\\ [alpha_r(z) is a non-residue mod q]_r over such primes (with all alpha_r(z) != 0).  Rank k excludes every
\\ nonempty T.  This is exact; for small p it is cross-checked by testing all 2^k - 1 subset products with
\\ nfeltispower over the cyclotomic field.
\\
\\ NNN ring at N = 15: omega_r = 2 sin(pi r/15) sqrt(2 + cos(2 pi r/15)), r = 1..7.  All 2 sin(pi r/15) lie in
\\ K = Q(zeta_60), and 4 (2 + cos(2 pi r/15)) = 8 + (y^r + y^-r)*2 with y = zeta_60^4 is integral in the same
\\ square class.  If these seven are independent in K*/K*^2, the seven radicals are K-linearly independent,
\\ so a rational relation sum c_r omega_r = 0 forces every c_r 2 sin(pi r/15) = 0, i.e. c = 0.

charcert(m, k, alpha, qmax = 10^6) = {
  my(rows = List(), q = 1, M);
  while(q < qmax,
    q += m;
    if(!isprime(q), next);
    my(w = znprimroot(q)^((q - 1)/m));
    for(j = 1, m - 1,
      if(gcd(j, m) != 1, next);
      my(v = vector(k, r, alpha(w^j, r)));
      if(vecmin(apply(t -> abs(lift(t)), v)) == 0, next);
      listput(rows, vector(k, r, (1 - kronecker(lift(v[r]), q))/2)));
    if(#rows < k, next);
    M = matrix(#rows, k, i, r, rows[i][r]);
    if(matrank(Mod(M, 2)) == k, return(q)));
  0;
}

exhaustive(p) = {
  my(P = polcyclo(p), nf = nfinit(P), z = Mod(x, P), k = (p - 1)/2, anysq = 0);
  my(alpha = vector(k, r, 4 - (z^r + z^(-r)) - (z^(2*r) + z^(-2*r))));
  forvec(e = vector(k, i, [0, 1]),
    if(vecsum(e) == 0, next);
    my(Q = Mod(1, P));
    for(i = 1, k, if(e[i], Q = Q * alpha[i]));
    if(nfeltispower(nf, lift(Q), 2), anysq = 1; break));
  !anysq;
}

kjump(z, r) = 4 - (z^r + z^(-r)) - (z^(2*r) + z^(-2*r));
nnn15(z, r) = 8 + 2*(z^(4*r) + z^(-4*r));
nnn12(z, r) = 8 + 2*(z^(4*r) + z^(-4*r));  \\ z of order 48: zeta_12 = zeta_48^4

allok = 1;
\\ collision-free primes (S and -S do not cover all nonzero residues mod p, i.e. p > 5)
plist = select(isprime, [7 .. 97]);
{
  for(i = 1, #plist,
    my(p = plist[i], k = (p - 1)/2, q = charcert(p, k, kjump));
    allok = allok && q > 0;
    printf("  C_%d({1,2}): %d distinct freqs, character matrix of rank %d (primes q <= %d) -> %s\n",
      p, k, k, q, if(q, "PASS", "FAIL")));
}
\\ cross-check of the certificate on the small primes by the exhaustive square test
cross = [7, 11, 13, 17];
{
  for(i = 1, #cross,
    my(ok = exhaustive(cross[i]));
    allok = allok && ok;
    printf("  C_%d({1,2}): exhaustive 2^%d - 1 subset products non-square: %s\n",
      cross[i], (cross[i] - 1)/2, if(ok, "PASS", "FAIL")));
}
q15 = charcert(60, 7, nnn15);
allok = allok && q15 > 0;
{
  printf("  NNN ring N=15: 7 radicals independent in Q(zeta_60)*/squares (primes q <= %d) -> %s\n",
    q15, if(q15, "PASS", "FAIL"));
}
\\ negative control: at N = 12 the NNN frequencies satisfy omega_3 = omega_6 = 2, so rank 6 must never occur
q12 = charcert(48, 6, nnn12, 10^5);
allok = allok && q12 == 0;
print("  NNN ring N=12 (relation omega_3 = omega_6): rank 6 not reached for q < 10^5 -> ", if(q12, "FAIL", "PASS"));
print("RESULT: " if(allok, "PASS", "FAIL") " -- saturation proven exactly for C_p({1,2}), primes 7..97, and the NNN ring at N=15");
quit;
