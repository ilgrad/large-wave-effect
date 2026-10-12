\\ Exact check of B_N = sqrt(N) for the prime squares N = 4, 9, 25, 49 (Remark rem:bn-psq).
\\
\\   gp -q -f exact/pari/bn_prime_square.gp </dev/null
\\
\\ B_N = sup_t sum_k |K_N(k,t)|, K_N(k,t) = (1/N) sum_r exp(-i t lambda_r) zeta_N^(r k),
\\ lambda_r = 2 - zeta_N^r - zeta_N^-r.  By Kronecker the closure of t -> (t lambda_r)_r mod 2 pi is the
\\ subtorus {theta : sum_r k_r theta_r = 0 mod 2 pi for every integer relation sum_r k_r lambda_r = 0}.
\\ For N = p^2 the eigenvalues are Q-dependent, so the prime-ring argument does not apply.  Here a
\\ point theta_r = 2 pi j_r / M of that subtorus is exhibited for each N and checked exactly:
\\   (1) j_0 = 0 and j_r = j_(N-r) (lambda_0 = 0, lambda_(N-r) = lambda_r);
\\   (2) sum_r k_r j_r = 0 (mod M) for every vector k of a Z-basis (matkerint) of the relation lattice;
\\   (3) |K(k)|^2 = 1/N for every k, in exact cyclotomic arithmetic.
\\ Then sum_k |K(k)| = sqrt(N) at that point, so B_N = sqrt(N) (with B_N <= sqrt(N) from unitarity).
\\ For odd p the point is searched in the ansatz j(p a + b) = pi(b) a + h(b) (mod p), M = p, pi a
\\ permutation with pi(0) = 0, pi(-b) = -pi(b); for N = 4 all j in (Z/4)^4 are searched, M = 4.
\\ Controls: no flat point for N = 6 on the (2 pi/6) Z grid, and test (2) rejects 12 of the 18 flat
\\ points of N = 9 on the (2 pi/9) Z grid (they are off the orbit closure).

relbasis(N) =
{
  my(P = polcyclo(N), m = N\2);
  matkerint(matconcat(vector(m, r, Vecrev(lift(Mod(2 - x^r - x^(N-r), P)), poldegree(P))~)));
}

\\ (1)-(3) for a phase vector j (j[r+1] = j_r, r = 0..N-1) with theta_r = 2 pi j_r / M
flatpoint(N, M, j, B) =
{
  my(L = lcm(N, M), P = polcyclo(L), z = Mod(x, P), m = N\2);
  if(j[1] % M, return(0));
  for(r = 1, N-1, if((j[r+1] - j[N-r+1]) % M, return(0)));
  for(c = 1, #B, if(sum(r = 1, m, B[r,c] * j[r+1]) % M, return(0)));
  for(k = 0, N-1,
    my(K = sum(r = 0, N-1, z^((-(L/M) * j[r+1] + (L/N) * r * k) % L)));
    my(Kc = sum(r = 0, N-1, z^(((L/M) * j[r+1] - (L/N) * r * k) % L)));
    if(K * Kc != N, return(0)));
  1;
}

\\ first flat point of the ansatz for N = p^2, p odd; [] if none
ansatz(p, B) =
{
  my(N = p^2, res = []);
  forperm(p, sg,
    my(pv = vector(p, b, sg[b] - 1), bad = (pv[1] != 0));
    for(b = 1, p-1, if((pv[p-b+1] + pv[b+1]) % p, bad = 1));
    if(bad, next);
    forvec(hv = vector((p-1)/2, i, [0, p-1]),
      my(h = vector(p));
      for(b = 1, (p-1)/2, h[b+1] = hv[b]; h[p-b+1] = (hv[b] - pv[b+1]) % p);
      my(j = vector(N, r1, my(r = r1 - 1); (pv[r % p + 1] * (r \ p) + h[r % p + 1]) % p));
      if(flatpoint(N, p, j, B), res = j; break(2))));
  res;
}

ok = 1;
{
my(B = relbasis(4), found = []);
forvec(v = vector(3, i, [0, 3]), my(j = [0, v[1], v[2], v[3]]); if(flatpoint(4, 4, j, B), found = j; break));
printf("N=4:  relation lattice rank %d, flat closure point theta = 2 pi j/4, j = %s -> %s\n",
  #B, found, if(#found, "PASS", "FAIL"));
ok = ok && #found && #B == 1;
}
{
foreach([3, 5, 7], p,
  my(B = relbasis(p^2), j = ansatz(p, B));
  printf("N=%d: relation lattice rank %d = (p-1)/2: %d; flat closure point theta = 2 pi j/%d, j = %s -> %s\n",
    p^2, #B, #B == (p-1)/2, p, j, if(#j, "PASS", "FAIL"));
  ok = ok && #j && #B == (p-1)/2);
}
\\ control: N = 6 is not a prime square; the same ansatz shape (M = 6, every j) finds no flat point
{
my(B = relbasis(6), any = 0);
forvec(v = vector(3, i, [0, 5]), my(j = [0, v[1], v[2], v[3], v[2], v[1]]); if(flatpoint(6, 6, j, B), any = 1; break));
printf("control N=6: no flat point with phases in (2 pi/6) Z -> %s\n", if(!any, "PASS", "FAIL"));
ok = ok && !any;
}
\\ control: the relation test bites -- with phases in (2 pi/9) Z, N = 9 has 18 flat points, 6 on the closure
{
my(B = relbasis(9), nf = 0, nc = 0);
forvec(v = vector(4, i, [0, 8]),
  my(j = vector(9, r1, my(r = r1 - 1); if(r == 0, 0, v[min(r, 9 - r)])));
  if(flatpoint(9, 9, j, []), nf++; if(flatpoint(9, 9, j, B), nc++)));
printf("control N=9: flat points with phases in (2 pi/9) Z: %d, of which on the closure: %d -> %s\n",
  nf, nc, if(nf == 18 && nc == 6, "PASS", "FAIL"));
ok = ok && nf == 18 && nc == 6;
}
print("RESULT: ", if(ok, "PASS", "FAIL"), " -- B_N = sqrt(N) certified for N = 4, 9, 25, 49");
quit;
