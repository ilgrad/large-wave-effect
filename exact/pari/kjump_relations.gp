\\ Exact Q-relation lattice of the k-jump circulant frequencies, and the mod-4 ceiling test on it.
\\
\\   gp -q -f exact/pari/kjump_relations.gp </dev/null
\\
\\ C_N(S) has frequencies omega_r = sqrt(alpha_r), r = 1..floor(N/2), with
\\     alpha_r = sum_{s in S} (2 - zeta^(s r) - zeta^(-s r)),   zeta = zeta_N,  alpha_r in K = Q(zeta_N).
\\ Decision procedure (no floating point in any verdict):
\\   (1) drop exact duplicates alpha_a = alpha_b (the collisions);
\\   (2) split the distinct alpha_r into square classes of K^*/K^*2: different Legendre-symbol rows at
\\       degree-1 primes q = 1 (mod N) prove different classes; equal rows are settled by nfroots;
\\   (3) square roots from different classes are K-linearly independent (Kummer), so every Q-relation
\\       lives inside one class; with class representative b, y_r = omega_r/omega_b lies in K and the
\\       relations are the integer kernel of the coordinate matrix of the y_r (matkerint: a Z-basis).
\\ The ceiling theorem then gives A_N = U_N when every basis relation has coefficient sum = 0 (mod 4).
\\ Floating point is used only to orient the sign of y_r and to re-check each relation at 150 digits.

default(debugmem, 0);
default(parisizemax, 2*10^9);
default(realprecision, 150);

freq(N, S, r) = sqrt(sum(i = 1, #S, 2 - 2*cos(2*Pi*S[i]*r/N)));
alphas(N, S) = my(T = polcyclo(N, 'y), z = Mod('y, T)); vector(N\2, r, sum(i = 1, #S, 2 - z^(S[i]*r) - z^(-S[i]*r)));

\\ Legendre row of every alpha_r at ncol degree-1 primes q (q = 1 mod N, image of zeta of order N mod q)
legrows(N, S, ncol) =
{
  my(k = N\2, M = matrix(k, ncol), cnt = 0, q = 1);
  while(cnt < ncol,
    q += N;
    if(!isprime(q), next);
    my(g = Mod(znprimroot(q), q)^((q-1)/N), col = vector(k), ok = 1);
    for(r = 1, k,
      my(a = sum(i = 1, #S, 2 - g^(S[i]*r) - g^(-S[i]*r)));
      if(a == 0, ok = 0; break);
      col[r] = kronecker(lift(a), q));
    if(ok, cnt++; for(r = 1, k, M[r, cnt] = col[r])));
  M;
}

\\ Z-basis of the relation lattice of the distinct frequencies, as a list of [r-vector, coefficient-vector]
relations(N, S) =
{
  my(T = polcyclo(N, 'y), d = poldegree(T), al = alphas(N, S), k = N\2, M = legrows(N, S, k + 40));
  my(keep = List());
  for(r = 1, k,
    my(dup = 0);
    for(j = 1, #keep, if(al[keep[j]] == al[r], dup = 1; break));
    if(!dup, listput(keep, r)));
  my(reps = List(), mem = List(), ys = List(), out = List());
  for(t = 1, #keep,
    my(r = keep[t], placed = 0);
    for(c = 1, #reps,
      my(b = reps[c]);
      if(M[r,] != M[b,], next);
      my(rt = nfroots(T, 'x^2 - lift(al[r] / al[b])));
      if(#rt == 0, next);
      my(yv = Vecrev(lift(rt[1]), d)~, yn = subst(lift(rt[1]), 'y, exp(2*Pi*I/N)));
      if(real(yn) < 0, yv = -yv);
      listput(mem[c], r); listput(ys[c], yv); placed = 1; break);
    if(!placed, listput(reps, r); listput(mem, List([r])); listput(ys, List([Vecrev(1, d)~]))));
  for(c = 1, #reps,
    if(#mem[c] < 2, next);
    my(Y = Mat(Vec(ys[c])), L = matkerint(Y * denominator(Y)), m = Vec(mem[c]));
    for(j = 1, #L,
      my(v = L[,j]~, chk = sum(i = 1, #m, v[i] * freq(N, S, m[i])));
      if(abs(chk) > 1e-140, error("relation fails numerically at N = ", N));
      listput(out, [m, v])));
  Vec(out);
}

mod4ok(R) = {for(i = 1, #R, if(vecsum(R[i][2]) % 4, return(0))); 1};

\\ (A) S = {1,3}, 5 <= N <= 126: relations exist exactly at N = 28, 56, 84, 112 (two each), all mod-4 clean
ok = 1;
{
my(relN = List(), allm4 = 1);
for(N = 5, 126,
  my(R = relations(N, [1,3]));
  if(#R, listput(relN, [N, #R]); allm4 = allm4 && mod4ok(R)));
my(okA = Vec(relN) == [[28, 2], [56, 2], [84, 2], [112, 2]] && allm4);
printf("[A] S={1,3}, N<=126: relations at [N, rank] %s; every basis sum = 0 mod 4: %d -> %s\n",
  Vec(relN), allm4, if(okA, "PASS", "FAIL"));
ok = ok && okA;
}

\\ (B) N = 28: the lattice is spanned by w1 - w3 - w9 + w14 and -w5 + w11 - w13 + w14, w14 = 2 sqrt 2
{
my(R = relations(28, [1,3]), m = R[1][1], e1 = vector(#m), e2 = vector(#m));
my(at = r -> select(x -> x == r, m, 1)[1]);
e1[at(1)] = 1; e1[at(3)] = -1; e1[at(9)] = -1; e1[at(14)] = 1;
e2[at(5)] = -1; e2[at(11)] = 1; e2[at(13)] = -1; e2[at(14)] = 1;
my(w = r -> freq(28, [1,3], r));
my(v1 = w(1) - w(3) - w(9) + w(14), v2 = -w(5) + w(11) - w(13) + w(14), v3 = w(14) - 2*sqrt(2));
my(okB = R[2][1] == m && mathnf(Mat([R[1][2]~, R[2][2]~])) == mathnf(Mat([e1~, e2~]))
  && abs(v1) < 1e-140 && abs(v2) < 1e-140 && abs(v3) < 1e-140);
printf("[B] N=28: lattice = <w1-w3-w9+w14, -w5+w11-w13+w14>, w14 = 2 sqrt 2 (residuals %.1e %.1e %.1e) -> %s\n",
  abs(v1), abs(v2), abs(v3), if(okB, "PASS", "FAIL"));
ok = ok && okB;
}

\\ (C) control N = 27: the lattice is zero, and the same two forms (w14 = w13 there) are far from zero
{
my(w = r -> freq(27, [1,3], r));
my(v1 = w(1) - w(3) - w(9) + w(14), v2 = -w(5) + w(11) - w(13) + w(14));
my(okC = #relations(27, [1,3]) == 0 && abs(v1) > 1e-3 && abs(v2) > 1e-3);
printf("[C] N=27: no relation; the N=28 forms give %.4f and %.4f -> %s\n", v1, v2, if(okC, "PASS", "FAIL"));
ok = ok && okC;
}

\\ (D) S = {1,2} and S = {1,2,3}, 5 <= N <= 126: no relation beyond the collisions
{
foreach([[1,2], [1,2,3]], S,
  my(bad = List());
  for(N = 5, 126, if(#relations(N, S), listput(bad, N)));
  printf("[D] S=%s, N<=126: N with a relation beyond collisions %s -> %s\n", S, Vec(bad), if(#bad, "FAIL", "PASS"));
  ok = ok && !#bad);
}

\\ (E) control, the ring S = {1}, 5 <= N <= 48: independent exactly at primes and powers of two, and the
\\     N = 12 lattice (paper: w5 = w1 + w3, w6 = 2 w2) has a basis sum != 0 mod 4
{
my(indep = List(), R12 = relations(12, [1]));
for(N = 5, 48, if(#relations(N, [1]) == 0, listput(indep, N)));
my(want = select(N -> isprime(N) || N == 1 << valuation(N, 2), [5..48]));
my(okE = Vec(indep) == want && #R12 == 2 && !mod4ok(R12));
printf("[E] ring S={1}, N<=48: independent at %s; N=12 rank %d, mod-4 obstructed -> %s\n",
  Vec(indep), #R12, if(okE, "PASS", "FAIL"));
ok = ok && okE;
}

print("RESULT: ", if(ok, "PASS", "FAIL"), " -- k-jump relation lattices decided exactly, N <= 126");
quit;
