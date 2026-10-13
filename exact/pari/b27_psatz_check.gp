\\ Independent check of the N = 27 Positivstellensatz certificate written by numerics/b27_psatz_certificate.py
\\ (paper, Theorem bn-saturation (4): H_27 has no flat point, so B_27 < sqrt 27).
\\ Nothing is taken from the generator but the certificate:
\\  (a) the coordinates c(r) of mu_r = 2 cos(2 pi r/27) in the basis 1, mu_1, ..., mu_8 are checked exactly in
\\      Q(zeta_27), together with the Q-linear independence of that basis;
\\  (b) R_l(X) = sum_k X^{c(k+l) - c(k)} is built from c and checked to be real on T^9 (R_l(1/X) = R_l(X));
\\  (c) the Gram matrix Q is expanded from one value per orbit of index pairs of B_3 (word length <= 3, in
\\      lexicographic order) under G = <multiplier 2> x <inversion>, and the multipliers from their orbits;
\\  (d) r = m^* Q m + sum_l g_l R_l + c is computed exactly, coefficient by coefficient, and
\\      lambda |B_3| + c > ||r||_1 is checked for the claimed lambda.  That lambda_min(Q) >= lambda is certified
\\      by numerics/verify_b27_psatz.py (an exact Cholesky residual) and is not re-checked here;
\\  (e) a perturbed Gram value and a perturbed multiplier must be rejected.
\\ Run from the repository root: gp -q -f exact/pari/b27_psatz_check.gp </dev/null   (a few minutes;
\\ B27_CERT=<file> checks another certificate).  PARI prints *** on an error but still exits 0 and goes on
\\ with the next line, so every flag below starts at 0 and the verdict tests each one for == 1.

default(debugmem, 0);
default(parisizemax, "4G");
coordok = 0; indep = 0; realok = 0; countok = 0; r1 = -1; bound = 0; rej1 = 0; rej2 = 0;
path = getenv("B27_CERT");
if(path == 0, path = "exact/pari/b27_psatz_cert.gp.gz");
cert = readvec(path);
HD = cert[1];
if(HD[1] != 27 || HD[2] != 3, error("not an N = 27, degree-3 certificate"));
SC = HD[3]; CC = HD[4]; LAM = HD[5];
{
  VALS = List(); ROWS = List();
  for(i = 2, #cert,
    my(v = cert[i]);
    if(v[1] == 1, for(t = 2, #v, listput(VALS, v[t])),
       v[1] == 2, listput(ROWS, v[2 .. #v]),
       error("unknown record")));
  VALS = Vec(VALS); ROWS = Vec(ROWS);
}
countok = #VALS == HD[6] && #ROWS == HD[7];
N = 27;

\\ (a) mu-coordinates: c(0) = 2 e_0, c(r) = e_r (1 <= r <= 8), c(9) = -e_0, c(9+i) = -(e_i + e_{9-i}), c(-r) = c(r)
cr(r) =
{
  my(v = vector(9), s = r % N);
  if(s > 13, s = N - s);
  if(s == 0, v[1] = 2, s <= 8, v[s + 1] = 1, s == 9, v[1] = -1, v[s - 8] = -1; v[19 - s] = -1);
  v;
}
Phi = polcyclo(N, x);
mu(r) = Mod(x^(r % N) + x^((N - r) % N), Phi);
beta = vector(9, i, if(i == 1, Mod(1, Phi), mu(i - 1)));
coordok = prod(r = 0, N - 1, mu(r) == sum(i = 1, 9, cr(r)[i] * beta[i]));
indep = matrank(matrix(9, poldegree(Phi), i, k, polcoef(lift(beta[i]), k - 1))) == 9;
print("mu-coordinates exact: ", coordok, ", 1, mu_1, ..., mu_8 independent: ", indep);

\\ (b) the autocorrelations, as matrices [exponent, multiplicity]
key(e) = sum(i = 1, 9, (e[i] + 32) * 64^(9 - i));
K0 = key(vector(9));
{
  RLS = vector(13);
  for(l = 1, 13,
    my(M = Map(), o, e);
    for(k = 0, N - 1,
      e = cr(k + l) - cr(k);
      if(mapisdefined(M, e, &o), mapput(M, e, o + 1), mapput(M, e, 1)));
    RLS[l] = Mat(M));
  realok = 1;
  for(l = 1, 13,
    my(T = RLS[l], M = Map(), o);
    for(i = 1, #T[, 1], mapput(M, T[i, 1], T[i, 2]));
    for(i = 1, #T[, 1], if(!mapisdefined(M, -T[i, 1], &o) || o != T[i, 2], realok = 0)));
}
print("R_l real on T^9: ", realok);
RKEY = vector(13, l, vector(#RLS[l][, 1], i, key(RLS[l][i, 1]) - K0));
RCF = vector(13, l, RLS[l][, 2]~);
lfold(l) = my(s = l % N); min(s, N - s);

\\ (c) B_3, the group, the orbit expansion
gens = vector(26, t, if(t <= 13, cr(t) - cr(0), cr(0) - cr(t - 13)));
{
  my(seen = Map(), layer = [vector(9)], allv = List([vector(9)]), nxt, f);
  mapput(seen, vector(9), 1);
  for(d = 1, 3,
    nxt = List();
    for(a = 1, #layer, for(t = 1, 26,
      f = layer[a] + gens[t];
      if(!mapisdefined(seen, f), mapput(seen, f, 1); listput(nxt, f); listput(allv, f))));
    layer = Vec(nxt));
  B3 = vecsort(Vec(allv));
}
n = #B3;
if(vecmax(apply(e -> vecmax(abs(e)), B3)) >= 16, error("B_3 exponents out of range"));
idx = Map(); for(i = 1, n, mapput(idx, B3[i], i));
M2 = matrix(9, 9); M2[1, 1] = 1;
for(b = 1, 8, M2[, b + 1] = cr(2 * b)~);
if(!prod(r = 0, N - 1, M2 * cr(r)~ == cr(2 * r)~), error("M2 is not the multiplier 2"));
Gm = vector(18, t, (-1)^((t - 1) % 2) * M2^((t - 1) \ 2));
P = vector(18, t, vectorsmall(n, i, mapget(idx, (Gm[t] * B3[i]~)~)));
cB = vector(n, i, key(B3[i]));

\\ (d) the residual: integer coefficients (scaled by 2^SC) accumulated in a map from monomial keys
resid1(vals, rows) =
{
  my(S = Map(), vis = vectorsmall(n^2), nv = 0, k = 0, o, v, p, q, id, kk, row, j, e, y, l, kf, F);
  for(i = 1, n, for(jj = i, n,
    if(!vis[(i - 1) * n + jj],
      k++; v = vals[k];
      for(t = 1, 18,
        p = P[t][i]; q = P[t][jj];
        id = (p - 1) * n + q;
        if(!vis[id], vis[id] = 1; nv++; kk = cB[q] - cB[p] + K0;
          if(mapisdefined(S, kk, &o), mapput(S, kk, o + v), mapput(S, kk, v)));
        id = (q - 1) * n + p;
        if(!vis[id], vis[id] = 1; nv++; kk = cB[p] - cB[q] + K0;
          if(mapisdefined(S, kk, &o), mapput(S, kk, o + v), mapput(S, kk, v)))))));
  if(k != #vals || nv != n^2, error("the pair orbits do not match the certificate"));
  for(rr = 1, #rows,
    row = rows[rr]; j = row[1]; e = row[2 .. 10]; y = row[11];
    for(t = 1, 18,
      l = lfold(2^((t - 1) \ 2) * j); kf = key((Gm[t] * e~)~);
      for(u = 1, #RKEY[l],
        kk = kf + RKEY[l][u]; v = y * RCF[l][u];
        if(mapisdefined(S, kk, &o), mapput(S, kk, o + v), mapput(S, kk, v)))));
  if(mapisdefined(S, K0, &o), mapput(S, K0, o + CC), mapput(S, K0, CC));
  F = Mat(S);
  sum(i = 1, #F[, 1], abs(F[i, 2]));
}
bound = LAM * n * 2^SC + CC;
r1 = resid1(VALS, ROWS);
print("|B_3| = ", n, ", pair orbits ", #VALS, ", multiplier orbits ", #ROWS);
print("||r||_1 = ", r1 * 1. / 2^SC, ", lambda |B_3| + c = ", bound * 1. / 2^SC);

\\ (e) mutations
bad = VALS; bad[1] += 10^5 * 2^SC;
rej1 = resid1(bad, ROWS) >= bound;
badr = ROWS; badr[1][11] += 1000 * 2^SC;
rej2 = resid1(VALS, badr) >= bound;
print("mutations rejected: ", rej1 + rej2, "/2");
{
  if(coordok == 1 && indep == 1 && realok == 1 && countok == 1 && CC >= 0 && type(r1) == "t_INT"
     && r1 >= 0 && r1 < bound && rej1 == 1 && rej2 == 1,
    print("RESULT: PASS (given lambda_min(Q) >= ", LAM, ")"),
    print("RESULT: FAIL"));
}
quit;
