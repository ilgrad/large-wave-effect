\\ Independent check of the N = 15 Nullstellensatz certificate written by exact/julia/c15_flat_certificate.jl:
\\     D * y^4 = sum_j N_j * m_j * f_{l_j}   (paper, Theorem bn-saturation (4): H_15 has no flat point).
\\ The flatness polynomials f_1..f_7 are rebuilt here from scratch, and the identity is checked by exact
\\ polynomial arithmetic over Z.
\\ Run: C15_CERT=/path/to/certificate.gp gp -q -f exact/pari/c15_flat_check.gp </dev/null
\\ PARI prints *** on an error but still exits 0, so check the output for *** as well as RESULT.

\\ coordinates of mu_n = 2 cos(2 pi n/15) = D_n(mu_1) in the basis (1, mu_1, mu_2 = t^2 - 2, mu_3 = t^3 - 3t)
Pm = t^4 - t^3 - 4*t^2 + 4*t + 1;
coord(n) = {
  my(d = lift(Mod(2*polchebyshev(n, 1, t/2), Pm)), a = vector(4, k, polcoef(d, k - 1, t)));
  [a[1] + 2*a[3], a[2] + 3*a[4], a[3], a[4]];
}
C = vector(15, n, coord(n - 1));
mono(e) = x^e[1] * y^e[2] * z^e[3] * u^e[4];
\\ f_l: the autocorrelation R_l = sum_k X^{c(k+l) - c(k)} times the least monomial making it a polynomial
cleared(l) = {
  my(E = vector(15, k, C[(k - 1 + l) % 15 + 1] - C[k]));
  my(lo = vector(4, j, vecmin(vector(15, k, E[k][j]))));
  sum(k = 1, 15, mono(E[k] - lo));
}
f = vector(7, l, cleared(l));
\\ every f_l has 15 distinct monomials, all with coefficient 1 (15 terms, coefficient sum 15)
nterms(p) = if(type(p) != "t_POL", p != 0, sum(i = 0, poldegree(p), nterms(polcoef(p, i))));
fok = prod(l = 1, 7, nterms(f[l]) == 15 && substvec(f[l], [x, y, z, u], [1, 1, 1, 1]) == 15);
path = getenv("C15_CERT");
if(path == 0, error("set C15_CERT to the certificate file written by exact/julia/c15_flat_certificate.jl"));
cert = readvec(path);
D = cert[1];
S = sum(j = 2, #cert, my(c = cert[j]); c[6] * mono(c[2..5]) * f[c[1]]);
\\ negative control: changing one coefficient by 1 must break the identity
neg = S + mono(cert[2][2..5]) * f[cert[2][1]] != D * y^4;
ok = fok && D != 0 && S == D * y^4 && neg;
{
  printf("  flatness polynomials rebuilt: %s; certificate: %d terms, D of %d digits\n",
    if(fok, "ok", "FAIL"), #cert - 1, #digits(abs(D)));
  printf("  sum_j N_j m_j f_{l_j} == D*y^4 exactly: %s; perturbed certificate rejected: %s\n",
    if(S == D * y^4, "yes", "no"), if(neg, "yes", "no"));
}
{
  print("RESULT: " if(ok, "PASS -- y^4 lies in the C_15 flatness ideal (checked in PARI/GP), so H_15 has no flat point",
    "FAIL -- the certificate does not check"));
}
quit;
