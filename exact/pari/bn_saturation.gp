\\ Exact checks for the infinite-time Schrodinger saturation B_N = sqrt(N) (Theorem bn-saturation).
\\ lambda_r = 2 - zeta_N^r - zeta_N^{-r}; H_N = orbit closure of t -> (exp(-i t lambda_r))_r.
\\ Everything below is integer / cyclotomic arithmetic; no floating point.
\\ Run: timeout 1200 gp -q -f exact/pari/bn_saturation.gp </dev/null
\\ PARI prints *** on an error but still exits 0, so check the output for *** as well as RESULT.

lamrows(N) = {
  my(P = polcyclo(N), d = poldegree(P));
  matrix(d, N\2, j, r, polcoef(lift(Mod(2 - 'x^r - 'x^((N - r) % N), P)), j - 1));
}
\\ integer relation lattice {n : sum_{r=1}^{N/2} n_r lambda_r = 0}, as columns
rellat(N) = matkerint(lamrows(N));
\\ |sum_r zeta_M^{e_r + (M/N) r k}|^2 == N for every k, exactly in Q(zeta_M)
flat(N, M, e) = {
  my(P = polcyclo(M), s, c, step = M / N);
  for(k = 0, N - 1,
    s = Mod(sum(r = 0, N - 1, 'x^((e[r + 1] + step*r*k) % M)), P);
    c = Mod(sum(r = 0, N - 1, 'x^((M - (e[r + 1] + step*r*k) % M) % M)), P);
    if(lift(s*c) != N, return(0)));
  1;
}
\\ fold a relation on Z_N onto the coordinates r = 1..N/2 (x_{-r} = x_r, x_0 = 1)
foldv(N, v) = {
  my(w = vector(N\2));
  for(r = 1, N - 1, if(v[r + 1], my(s = min(r, N - r)); w[s] += v[r + 1]));
  w~;
}
\\ folded coset relations prod_{coset j} x = prod_{coset 0} x for N = p^a, j = 1..N/p-1
cosetrels(N, p) = {
  my(q = N/p, cols = List(), v0 = vector(N, s, ((s - 1) % q) == 0));
  for(j = 1, q - 1, my(vj = vector(N, s, ((s - 1) % q) == j)); listput(cols, foldv(N, vj - v0)));
  matconcat(Vec(cols));
}
\\ periodic autocorrelation at shifts 1..M-1 of an integer vector
acf(v) = my(M = #v); vector(M - 1, t, sum(i = 0, M - 1, v[i + 1]*v[(i + t) % M + 1]));
\\ number of +-1 periodic complementary pairs (u, v) of length M with u_{-i} = u_i, v_{-1-i} = v_i
pgp(M) = {
  my(m = Map(), hits = 0, U, V);
  forvec(s = vector(M\2 + M%2, i, [0, 1]),
    V = vector(M, i, my(j = min(i, M + 1 - i)); 1 - 2*s[j]);
    mapput(m, acf(V), 1));
  forvec(s = vector(M\2 + 1, i, [0, 1]),
    U = vector(M, i, my(j = min(i, M + 2 - i)); 1 - 2*s[j]);
    if(mapisdefined(m, -acf(U)), hits++));
  hits;
}
ok1 = ok2 = ok3 = ok4 = ok5 = 0;

{
  my(bad = select(N -> matrank(lamrows(N)) != eulerphi(N)/2, [3..200]));
  print("[1] rank_Q{lambda_r} == phi(N)/2 for 3 <= N <= 200: ", if(#bad, Str("FAIL ", bad), "yes"));
  ok1 = !#bad;
}
{
  my(pr = 1);
  forprime(p = 3, 101,
    if(#rellat(p) != 0, pr = 0);
    if(!flat(p, p, vector(p, r, (r - 1)^2 % p)), pr = 0));
  print("[2] primes 3..101: H_p = full torus and Gauss chirp exactly flat: ", if(pr, "yes", "FAIL"));
  ok2 = pr;
}
{
  my(good = 1);
  foreach([[9, 3], [25, 5], [27, 3], [49, 7], [81, 3], [121, 11], [125, 5], [243, 3]], c,
    my(N = c[1], p = c[2]);
    if(mathnf(rellat(N)) != mathnf(cosetrels(N, p)), good = 0));
  print("[3] N = p^a in {9,25,27,49,81,121,125,243}: relation lattice = span of coset relations: ",
        if(good, "yes", "FAIL"));
  ok3 = good;
}
{
  my(fr = 1);
  foreach([3, 5, 7, 11, 13], p,
    my(N = p^2, h = (p + 1)/2, e, K);
    e = vector(N, s, my(r = s - 1, i = r\p, j = r%p); (-i*j - h*j) % p);
    if(e[1] != 0 || sum(r = 1, N - 1, e[r + 1] != e[N - r + 1]), fr = 0);
    if(!flat(N, N, p*e), fr = 0);
    K = rellat(N);
    for(c = 1, #K, if(sum(r = 1, N\2, K[r, c]*e[r + 1]) % p, fr = 0)));
  print("[4] p in {3,5,7,11,13}: even Frank point x_{pi+j} = zeta_p^(-ij-hj) exactly flat, in H_{p^2}: ",
        if(fr, "yes", "FAIL"));
  ok4 = fr;
}
{
  my(pg = 1, res = List(), h);
  for(M = 2, 18, h = pgp(M); listput(res, [M, h]); if((h > 0) != (M == 2), pg = 0));
  print("[5] symmetric periodic Golay pairs [M, #found], 2 <= M <= 18: ", Vec(res),
        "  only M = 2: ", if(pg, "yes", "FAIL"));
  ok5 = pg;
}
{
  my(ch = List());
  for(N = 3, 80,
    my(K = rellat(N), md = if(N % 2, N, 2*N), good = 1);
    for(c = 1, #K, if(sum(r = 1, N\2, K[r, c]*r^2) % md, good = 0; break));
    if(good, listput(ch, N)));
  print("[6] N <= 80 whose orbit closure contains the quadratic chirp: ", Vec(ch),
        "  (= primes: ", Vec(ch) == primes([3, 80]), ")");
}
print("RESULT: ", if(ok1 && ok2 && ok3 && ok4 && ok5, "PASS", "FAIL"));
quit
