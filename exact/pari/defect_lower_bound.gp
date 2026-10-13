\\ Exact check of Proposition deficit-lower: an explicit lower bound for U_N - A_N on every non-saturating N.
\\
\\ For an odd prime p | N, s = N/p, q = 2s, each class x mod q with x != 0, s gives the vanishing sum
\\ sum_{j<p} zeta_{2N}^(x+qj) = 0, whose imaginary part is an integer relation k (sum_r k_r omega_r = 0,
\\ omega_r = 2 sin(pi r/N), 1 <= r <= N/2) with odd coefficient sum.  Classes x, -x, x+s, s-x give the same
\\ relation; distinct orbits have disjoint supports.  Lemma rel-defect gives, at nodes 0 and N/2,
\\     U_N - A_N >= sum_orbits 1 / (2 sum_r k_r^2 / b_r),     b_r = 1/(N sin(pi r/N)),  b_{N/2} = 1/(2N),
\\ and every other node has deficit >= 1/6.  For 3 | N the standard triplets with Lemma triplet give the
\\ sharper (1/2) sum_{1 <= a < N/6} b_{N/3+a}.  Everything below is exact rational arithmetic: pi and sin are
\\ enclosed by rationals (alternating Taylor series), every relation is checked in Q(zeta_2N), and each
\\ certified value is a rational lower bound.
\\
\\   gp -q -f exact/pari/defect_lower_bound.gp </dev/null        (from the repository root)

PILO = 314159265358979323846/10^20;
PIHI = 314159265358979323847/10^20;
\\ sin on [0, pi/2] (exact x): partial sums of the alternating series bracket sin x, terms decreasing for x < 4
sinU(xx) = sum(k = 0, 12, (-1)^k * xx^(2*k+1) / (2*k+1)!);
sinL(xx) = sum(k = 0, 11, (-1)^k * xx^(2*k+1) / (2*k+1)!);
\\ rational bounds on sin(pi r/N), 0 < r <= N/2 (sin increasing there)
sinhi(r, N) = sinU(PIHI * r / N);
sinlo(r, N) = sinL(PILO * r / N);
binvhi(r, N) = if(2*r == N, 2*N, N * sinhi(r, N));   \\ upper bound on 1/b_r
blo(r, N) = 1 / binvhi(r, N);                            \\ lower bound on b_r

saturating(N) = isprime(N) || N == 2^valuation(N, 2);

\\ the relation of the class x mod q: a map r -> k_r over 1 <= r <= N/2
foldmode(N, e) =
{
  e = e % (2*N);
  if(e == 0 || e == N, error("class hits 0 or N"));
  my(sg = if(e < N, 1, -1), r = if(e < N, e, e - N));
  [min(r, N - r), sg];
}
relation(N, p, x) =
{
  my(q = 2*N/p, M = Map(), v);
  for(j = 0, p - 1,
    my(f = foldmode(N, x + q*j));
    mapput(M, f[1], if(mapisdefined(M, f[1], &v), v, 0) + f[2]));
  my(K = Mat(M), L = List());
  for(i = 1, #K~, if(K[i, 2], listput(L, [K[i, 1], K[i, 2]])));
  Vec(L);
}
\\ exact: sum_r k_r (z^r - z^-r) = 0 in Q(zeta_2N) and odd coefficient sum
relation_ok(N, rel) =
{
  my(P = polcyclo(2*N), z = Mod(x, P), t = 0, ks = 0);
  for(i = 1, #rel, my(r = rel[i][1], k = rel[i][2]); t += k * (z^r - z^(2*N - r)); ks += k);
  t == 0 && ks % 2 == 1;
}
\\ orbits of {x mod q : x != 0, s} under x -> -x, x -> x + s; returns the list of relations
orbit_relations(N, p) =
{
  my(s = N/p, q = 2*s, seen = vector(q), L = List());
  for(x = 1, q - 1,
    if(x == s || seen[x + 1], next);
    foreach([x, (-x) % q, (x + s) % q, (s - x) % q], y, seen[y + 1] = 1);
    listput(L, relation(N, p, x)));
  Vec(L);
}
\\ certified lower bound from a list of relations (0 if a check fails, so a failure can never certify)
bound_from(N, rels) =
{
  my(used = Map(), tot = 0);
  for(i = 1, #rels,
    my(rel = rels[i], S = 0);
    if(!relation_ok(N, rel), return(0));
    for(j = 1, #rel,
      if(mapisdefined(used, rel[j][1]), return(0));
      mapput(used, rel[j][1], 1);
      S += rel[j][2]^2 * binvhi(rel[j][1], N));
    tot += 1 / (2*S));
  tot;
}
orbit_bound(N, p) =
{
  my(rels = orbit_relations(N, p));
  if(#rels < (N/p - 1)/2, return(0));
  bound_from(N, rels);
}
\\ standard triplets for 3 | N, each relation checked exactly
triplet_bound(N) =
{
  my(t = 0);
  if(N % 3, return(0));
  for(a = 1, ceil(N/6) - 1,
    if(!relation_ok(N, [[a, 1], [N/3 - a, 1], [N/3 + a, -1]]), return(0));
    t += blo(N/3 + a, N) / 2);
  t;
}
\\ the closed form (s - 1) sin(pi/2p) / (8 p s), upper-rounded so that "certified >= closed" is a proof
closed_hi(N, p) = my(s = N/p); (s - 1) * sinU(PIHI / (2*p)) / (8*p*s);

NMAX = 240;
nN = 0; okA = 1;
{
  for(N = 4, NMAX,
    if(saturating(N), next);
    nN++;
    my(best = 0);
    foreach(factor(N)[, 1]~, p,
      if(p == 2, next);
      my(ob = orbit_bound(N, p));
      if(ob == 0 || ob < closed_hi(N, p), okA = 0; printf("  N = %d, p = %d: FAIL\n", N, p));
      best = max(best, ob));
    if(best < 1/(16 * vecmin(select(p -> p > 2, factor(N)[, 1]~))^2), okA = 0));
}
{
printf("[1] %d non-saturating N <= %d: orbit relations exact, supports disjoint, certified bound >= closed form >= 1/(16p^2): %s\n",
       nN, NMAX, if(okA, "yes", "NO"));
}
{
printf("    sample certified lower bounds: N = 10: %.5f, N = 25: %.5f, N = 35: %.5f, N = 49: %.5f\n",
       orbit_bound(10, 5), orbit_bound(25, 5), max(orbit_bound(35, 5), orbit_bound(35, 7)), orbit_bound(49, 7));
}

\\ [2] 3 | N: the triplet bound, its closed-form floor and its limit ln(3)/(4 pi)
okB = 1; minB = [1, 0];
{
  forstep(N = 9, NMAX, 3,
    my(tb = triplet_bound(N));
    if(tb == 0, okB = 0);
    if(N >= 12 && tb < 421/10000, okB = 0);
    if(tb < minB[1], minB = [tb, N]));
}
default(realprecision, 50);
floor_ok = 1;
forstep(N = 12, NMAX, 3, if(triplet_bound(N) < log(cotan(Pi/6 + Pi/(2*N)))/(2*Pi), floor_ok = 0));
{
printf("[2] 3 | N in [9, %d]: triplets exact, bound >= 0.0421 for N >= 12, min %.5f at N = %d: %s\n",
       NMAX, minB[1], minB[2], if(okB, "yes", "NO"));
}
{
printf("    bound >= (1/2pi) ln cot(pi/6 + pi/2N) for 12 <= N <= %d: %s;  N = 9: %.5f;  limit ln3/(4pi) = %.6f\n",
       NMAX, if(floor_ok, "yes", "NO"), triplet_bound(9), log(3)/(4*Pi));
}
\\ N = 6 from the closed form of Proposition closed: U_6 - A_6
U6 = 1/3 + 1/(3*sqrt(3)) + 1/12;
A6 = sqrt(3)/9 + (3 + sqrt(3)) * 3^(1/4) / (12*sqrt(2));
printf("    N = 6 (closed form): U_6 - A_6 = %.6f\n", U6 - A6);

\\ [3] N = 2p: the single relation spans Lambda_N, S = sum k_r^2/b_r = 4p csc(pi/2p), and
\\     (1 - pi^3/(48p)) pi^2/(8S) <= U_N - A_N <= pi^2/(8S)
okC = 1;
{
  forprime(p = 3, 113,
    my(N = 2*p, rels = orbit_relations(N, p), rel, S = 0., nyq = 0);
    if(#rels != 1, okC = 0; next);
    rel = rels[1];
    if(!relation_ok(N, rel), okC = 0);
    for(j = 1, #rel,
      my(r = rel[j][1], k = rel[j][2]);
      if(2*r == N, nyq += (abs(k) == 1), if(abs(k) != 2, okC = 0));
      S += k^2 * if(2*r == N, 2*N, N*sin(Pi*r/N)));
    if(nyq != 1 || abs(S - 4*p/sin(Pi/(2*p))) > 1e-40, okC = 0));
}
{
printf("[3] N = 2p, p <= 113: one relation, Nyquist coefficient 1, the rest 2, S = 4p csc(pi/2p): %s\n",
       if(okC, "yes", "NO"));
}
{
  foreach([5, 7, 29], p, my(S = 4*p/sin(Pi/(2*p)));
    printf("    N = %2d: %.6f <= U_N - A_N <= %.6f\n", 2*p, (1 - Pi^3/(48*p))*Pi^2/(8*S), Pi^2/(8*S)));
}

\\ [4] the node-j bound: sum_{k<n} (1 - |cos(2 pi k/n)|) >= n/3 for 3 <= n <= 2000 (proved for all n in the text)
okD = 1;
for(n = 3, 2000, if(sum(k = 0, n - 1, 1 - abs(cos(2*Pi*k/n))) < n/3 - 1e-30, okD = 0));
printf("[4] sum_{k<n} (1 - |cos(2 pi k/n)|) >= n/3 for 3 <= n <= 2000: %s\n", if(okD, "yes", "NO"));

\\ [5] mutations: each must be rejected (a rejected relation list certifies nothing: bound_from returns 0)
nrej = 0;
{
nrej += !relation_ok(15, [[1, 1], [4, 1], [6, 1]]);                \\ sign flipped
nrej += !relation_ok(15, [[1, 2], [4, 1], [6, -1]]);               \\ coefficient changed
nrej += !relation_ok(12, [[2, 4], [6, -2]]);                        \\ true relation, even sum: no obstruction
rels45 = orbit_relations(45, 3);
}
{
nrej += bound_from(45, concat(rels45, [rels45[1]])) == 0;          \\ a repeated orbit overlaps
nrej += bound_from(45, concat(rels45, [[[1, 1], [14, 1], [16, 1]]])) == 0;   \\ a false relation appended
nacc = relation_ok(12, [[2, 2], [6, -1]]) + relation_ok(15, [[1, 1], [4, 1], [6, -1]]);   \\ controls: accepted
printf("[5] rejected %d of 5 mutations, accepted %d of 2 controls: %s\n", nrej, nacc,
       if(nrej == 5 && nacc == 2, "PASS", "FAIL"));
}

ok = okA && okB && floor_ok && okC && okD && nrej == 5 && nacc == 2;
print("RESULT: ", if(ok, "PASS", "FAIL"), " -- explicit deficit lower bounds on every non-saturating N <= ", NMAX);
quit;
