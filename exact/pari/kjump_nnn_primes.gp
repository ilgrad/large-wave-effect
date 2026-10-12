\\ Exact inputs of Proposition prop:kjump-prime: C_p({1,2}) saturates, A_p = U_p, for every prime p >= 7.
\\
\\   gp -q -f exact/pari/kjump_nnn_primes.gp </dev/null
\\
\\ zeta = zeta_p, c_r = zeta^r + zeta^-r, alpha_r = 4 - c_r - c_2r, k = (p-1)/2.  The proof uses:
\\   (1) alpha_r = -zeta^-r (1 - zeta^r)^2 (3 + c_r), with zeta^-r = (zeta^(-r(p+1)/2))^2 a square,
\\       so alpha_r = -(3 + c_r) modulo squares of Q(zeta_p);
\\   (2) prod_{r=1..k} (3 + c_r) = F_p (Fibonacci), every factor positive;
\\   (3) F_p is not a square and p does not divide F_p (F_p = (p/5) mod p);
\\   (4) c_r - c_s = (zeta^r - zeta^s)(1 - zeta^(-r-s)).
\\ (1) and (2) are checked as exact identities in Z[x]/Phi_p and (3) for every prime 7 <= p <= 211;
\\ (4), a Laurent-polynomial identity, for p <= 43.  As an independent check of the conclusion, (5)
\\ certifies that no nonempty subset product of the alpha_r is a square: the Legendre symbols of the
\\ alpha_r at degree-1 primes q = 1 (mod p) have full F_2-rank k.  A mutated alpha (4 replaced by 5) must break (1).

default(debugmem, 0);
default(parisizemax, 10^9);

ident(p, four) =
{
  my(T = polcyclo(p), z = Mod(x, T), k = (p-1)/2, ok = 1, P = Mod(1, T));
  for(r = 1, k,
    my(c = z^r + z^-r, al = four - c - (z^(2*r) + z^(-2*r)));
    ok = ok && al == -z^-r * (1 - z^r)^2 * (3 + c) && z^-r == (z^(-r*(p+1)/2))^2;
    P *= 3 + c;
    if(p <= 43, for(s = r + 1, k, ok = ok && c - (z^s + z^-s) == (z^r - z^s) * (1 - z^(-r-s)))));
  ok && P == fibonacci(p) && real(prod(r = 1, k, 3 + 2*cos(2*Pi*r/p))) > 0;
}

fibfacts(p) = my(F = fibonacci(p)); !issquare(F) && F % p && Mod(F, p) == kronecker(p, 5);

charrank(p) =
{
  my(k = (p-1)/2, ncol = k + 24, M = matrix(k, ncol), cnt = 0, q = 1);
  while(cnt < ncol,
    q += p;
    if(!isprime(q), next);
    my(g = Mod(znprimroot(q), q)^((q-1)/p), col = vector(k), good = 1);
    for(r = 1, k,
      my(a = 4 - (g^r + g^-r) - (g^(2*r) + g^(-2*r)));
      if(a == 0, good = 0; break);
      col[r] = kronecker(lift(a), q) == -1);
    if(good, cnt++; for(r = 1, k, M[r, cnt] = col[r])));
  matrank(Mod(M, 2)) == k;
}

{
my(np = 0, okI = 1, okF = 1, okC = 1);
forprime(p = 7, 211,
  np++;
  okI = okI && ident(p, 4);
  okF = okF && fibfacts(p);
  okC = okC && charrank(p));
my(mut = !ident(7, 5) && !ident(13, 5));
printf("[1,2,4] alpha_r = -(3+c_r) mod squares, prod = F_p (%d primes 7..211), c_r - c_s (p <= 43): %s\n",
  np, if(okI, "PASS", "FAIL"));
printf("[3] F_p not a square, F_p = (p/5) mod p: %s\n", if(okF, "PASS", "FAIL"));
printf("[5] Legendre-character F_2-rank k (no subset product a square): %s\n", if(okC, "PASS", "FAIL"));
printf("mutation alpha with 5 in place of 4 rejected: %s\n", if(mut, "PASS", "FAIL"));
print("RESULT: ", if(np == 44 && okI && okF && okC && mut, "PASS", "FAIL"),
  " -- C_p({1,2}) saturation inputs verified for primes 7..211");
}
quit;
