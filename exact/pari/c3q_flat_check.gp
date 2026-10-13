\\ Independent check of the N = 3q Nullstellensatz certificates written by numerics/bn_flat_certificate_3q.py
\\ (paper, Theorem bn-saturation (4): H_15 and H_21 have no flat point, so B_15 < sqrt 15, B_21 < sqrt 21).
\\ Nothing is taken from the generator but the certificate:
\\  (a) the y-parametrization of H_N is rebuilt: the exponent pattern of every y_j - y_{q-1} is even, vanishes
\\      at n = 0 up to a constant, and has DFT zero at every nonzero non-unit (exact, in Q(zeta_N)); the
\\      patterns span the whole space V of Theorem bn-saturation (1), whose dimension is computed here;
\\  (b) the autocorrelations A(b, s) of the 3 x q array are computed from their definition, and checked to
\\      equal the formulas I_s, II_s displayed in the paper;
\\  (c) the certificate is expanded as  sum (multiplier) * A(shift)  over nonzero shifts and must equal 1
\\      exactly in Z[L], L = Z^q / Z(1, ..., 1), after clearing the common denominator.
\\ Run: C3Q_CERT=/path/to/certificate.txt gp -q -f exact/pari/c3q_flat_check.gp </dev/null
\\ PARI prints *** on an error but still exits 0, so check the output for *** as well as RESULT.

\\ the N = 21 expansion needs about 60 MB of PARI stack; -f leaves parisizemax at 0 (no growth).  debugmem 0
\\ silences the stack-growth warnings, which print *** like an error
default(debugmem, 0);
default(parisizemax, "2G");
path = getenv("C3Q_CERT");
if(path == 0, error("set C3Q_CERT to the file written by numerics/bn_flat_certificate_3q.py --out"));
cert = readvec(path);
q = cert[1][1]; D = cert[1][2]; DEN = cert[1][3];
N = 3 * q;
if(!isprime(q) || q == 3, error("q must be a prime other than 3"));
ev(a) = vector(q, j, j - 1 == a % q);
\\ CRT: n <-> (b, a) = (n mod 3, n mod q).  Exponent vector of the array entry x(b, a): x(1,a) = y_a, x(2,a) = y_{-a}, x(0,a) = 1/(y_a y_{-a})
X(b, a) = {
  b = b % 3; a = a % q;
  if(b == 1, ev(a), b == 2, ev(-a), -ev(a) - ev(-a));
}

\\ (a) the parametrization
pattern(j) = vector(N, n, my(e = X((n - 1) % 3, (n - 1) % q)); e[j + 1] - e[q]);
Phi = polcyclo(N, z);
dft(f, k) = Mod(sum(n = 0, N - 1, f[n + 1] * z^(k * n % N)), Phi);
nonunits = select(k -> gcd(k, N) > 1, [1 .. N - 1]);
even(f) = prod(n = 1, N - 1, f[n + 1] == f[N - n + 1]);
\\ f = pattern - pattern(0): the constant shift only moves the DFT at k = 0
fixed(j) = my(p = pattern(j)); p - vector(N, n, p[1]);
{
  parok = prod(j = 0, q - 2, my(f = fixed(j)); even(f) && f[1] == 0
    && prod(i = 1, #nonunits, dft(f, nonunits[i]) == 0));
}
\\ V = {f real : f even, f(0) = 0, f^ = 0 at every nonzero non-unit}: rational conditions, coefficient-wise
condrows = List();
{
  for(n = 1, N - 1, listput(condrows, vector(N, m, (m == n + 1) - (m == N - n + 1))));
  listput(condrows, vector(N, m, m == 1));
  for(i = 1, #nonunits,
    my(cols = vector(N, n, Vecrev(lift(Mod(z^(nonunits[i] * (n - 1) % N), Phi)), poldegree(Phi))));
    for(r = 1, poldegree(Phi), listput(condrows, vector(N, n, cols[n][r]))));
}
dimV = N - matrank(matconcat(Vec(condrows)~));
span = matrank(Mat(vector(q - 1, j, fixed(j - 1)~)));
parok = parok && dimV == q - 1 && span == q - 1;

\\ Laurent polynomials in Z[L] as Maps: canonical exponent vector (minimum 0) -> integer coefficient
canon(e) = e - vector(q, j, vecmin(e));
addto(~M, e, c) = {
  my(k = canon(e), v);
  if(mapisdefined(M, k, &v), mapput(~M, k, v + c), mapput(~M, k, c));
}
terms(M) = my(T = Mat(M)); select(t -> t[2] != 0, vector(matsize(T)[1], i, [T[i, 1], T[i, 2]]));
same(M1, M2) = terms(M1) == terms(M2);

\\ (b) autocorrelations A(b, s) = sum_n x(n + (b, s)) * conj x(n), conj y = 1/y on the torus
A(b, s) = {
  my(M = Map());
  for(bb = 0, 2, for(a = 0, q - 1, addto(~M, X(bb + b, a + s) - X(bb, a), 1)));
  M;
}
Iform(s) = {
  my(M = Map());
  for(a = 0, q - 1,
    addto(~M, ev(a + s) + ev(a) + ev(-a), 1);
    addto(~M, -ev(a + s) - ev(a) - ev(-a), 1);
    addto(~M, ev(-a - s) - ev(a), 1));
  M;
}
IIform(s) = {
  my(M = Map());
  for(a = 0, q - 1,
    addto(~M, ev(a) + ev(-a) - ev(a + s) - ev(-a - s), 1);
    addto(~M, ev(a + s) - ev(a), 1);
    addto(~M, ev(a) - ev(a + s), 1));
  M;
}
formok = prod(s = 0, q - 1, same(A(1, s), Iform(s))) && prod(s = 1, (q - 1) \ 2, same(A(0, s), IIform(s)));

\\ (c) the certificate.  u acts by a -> u a; inv is y -> 1/y.  The multipliers are orbit sums:
\\   kind 0:  sum over the (Z/q)^* x {1, inv} orbit of m, times A(1, 0);
\\   kind 1:  for each unit u, the u-image of {m, inv m}, times A(1, u);
\\   kind 2:  for u = 1 .. (q-1)/2, the u-image of {m, inv m, (-1).m, inv (-1).m}, times A(0, u).
act(u, e) = my(o = vector(q)); for(a = 0, q - 1, o[(u * a) % q + 1] = e[a + 1]); o;
inv(e) = canon(-e);
Gorbit(m) = Set(concat(vector(q - 1, u, [canon(act(u, m)), inv(act(u, m))])));
TA10 = terms(A(1, 0));
TA1 = vector(q - 1, u, terms(A(1, u)));
TA0 = vector((q - 1) \ 2, u, terms(A(0, u)));
\\ add c * sum_{f in orbit} y^(act(u, f)) * T to M
addprod(~M, orbit, u, T, c) = {
  for(i = 1, #orbit, my(f = act(u, orbit[i]));
    for(k = 1, #T, addto(~M, f + T[k][1], c * T[k][2])));
}
addentry(~M, ent, c) = {
  my(kind = ent[1], m = ent[2]);
  if(kind == 0, addprod(~M, Gorbit(m), 1, TA10, c),
     kind == 1, my(o = Set([canon(m), inv(m)])); for(u = 1, q - 1, addprod(~M, o, u, TA1[u], c)),
     kind == 2, my(o = Set([canon(m), inv(m), canon(act(q - 1, m)), inv(act(q - 1, m))]));
                for(u = 1, (q - 1) \ 2, addprod(~M, o, u, TA0[u], c)),
     error("unknown kind"));
}
S = Map();
for(j = 2, #cert, addentry(~S, cert[j], cert[j][3]));
TS = terms(S);
idok = DEN != 0 && #TS == 1 && TS[1][1] == vector(q) && TS[1][2] == DEN;
\\ negative control: raising one numerator by 1 adds that entry's column, which must be nonzero
P = Map(); addentry(~P, cert[2], 1);
neg = #terms(P) > 0;

ok = parok && formok && idok && neg;
{
  printf("  N = %d: y-parametrization spans V (dim %d = q-1), DFT zero at the non-units: %s\n",
    N, dimV, if(parok, "ok", "FAIL"));
  printf("  autocorrelations equal the displayed I_s, II_s: %s\n", if(formok, "ok", "FAIL"));
  printf("  certificate: %d orbit coefficients, ball radius %d, denominator of %d digits\n",
    #cert - 1, D, #digits(abs(DEN)));
  printf("  sum (multiplier) * A(shift) == DEN * 1 exactly in Z[L]: %s; perturbed certificate rejected: %s\n",
    if(idok, "yes", "no"), if(neg, "yes", "no"));
}
{
  print("RESULT: ", if(ok, Str("PASS -- 1 lies in the ideal of the C_", N,
    " autocorrelations (checked in PARI/GP), so H_", N, " has no flat point"),
    "FAIL -- the certificate does not check"));
}
quit;
