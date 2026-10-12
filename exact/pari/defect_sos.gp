\\ Exact certificates for the composite defect U_N - A_N, odd N = 9, 15, 21 (Remark rem:defect).
\\
\\   gp -q -f exact/pari/defect_sos.gp </dev/null        (from the repository root)
\\
\\ A_N = sup_t G_N(0,t), G_N(0,t) = sum_{r=1}^m b_r sin(omega_r t), omega_r = 2 sin(pi r/N), b_r = 2/(N omega_r),
\\ m = (N-1)/2, U_N = sum_r b_r.  The data (exact/pari/defect_sos_certs.gp, written by
\\ numerics/build_defect_sos_certificates.py) give an integer matrix C, a monomial basis B and a Gaussian-
\\ rational Gram matrix Q.  Everything below is exact rational arithmetic:
\\   (1) omega_1..omega_M (M = phi(2N)/2) are Q-independent and omega_r = sum_s C[r,s] omega_s for every r,
\\       checked on i*omega_r = zeta^r - zeta^-r in Q(zeta_2N).  By Kronecker the closure of the orbit
\\       (omega_r t)_r is {(<C_r, psi>)_r : psi in T^M}, so A_N = max_psi F(psi), F = sum_r b_r sin<C_r, psi>.
\\   (2) b_r is enclosed in [blo_r, bhi_r]: x_r = 2cos(2 pi r/N) is isolated by Sturm counts among the m+1
\\       distinct roots of 2 T_N(x/2) - 2, and b_r^2 = 4/(N^2 (2 - x_r)).
\\   (3) Q is positive definite: every pivot of its fraction-free LDL* (Bareiss) is real and positive.
\\   (4) With z_s = e^{i psi_s}, v = (z^{B_a})_a and Ft the polynomial F with weights bt_r = (blo_r+bhi_r)/2,
\\       r := c - Ft - v*Qv with c = tr Q is computed coefficient by coefficient.  On the torus
\\       F <= Ft + sum_r (bhi_r - blo_r) and Ft = c - v*Qv - r <= c + ||r||_1, so
\\           U_N - A_N >= sum_r blo_r - c - ||r||_1 - sum_r (bhi_r - blo_r) =: delta.
\\ The script proves delta >= claim for each N, then runs mutations that must all be rejected.

read("exact/pari/defect_sos_certs.gp");
default(realprecision, 80);

\\ (2) Sturm enclosures of b_1..b_m, width about 10^-40
benclose(N) =
{
  my(m = (N-1)/2, f = 2*subst(polchebyshev(N), x, x/2) - 2, D = 10^40, iv = vector(m+1), blo = vector(m), bhi = vector(m));
  if(polsturm(f) != m + 1, error("2T_N(x/2)-2 should have m+1 distinct real roots"));
  iv[1] = [2, 2];
  for(k = 1, m,
    my(xk = 2*cos(2*Pi*k/N));
    iv[k+1] = [floor(xk*D - 1)/D, ceil(xk*D + 1)/D]);
  for(k = 0, m,
    if(polsturm(f, iv[k+1]) != 1, error("interval without a unique root"));
    if(k > 0 && iv[k+1][2] >= iv[k][1], error("intervals overlap")));
  \\ m+1 disjoint, decreasing intervals each holding one of the m+1 roots: interval k holds 2cos(2 pi k/N)
  for(r = 1, m,
    my(L = 4/(N^2*(2 - iv[r+1][1])), H = 4/(N^2*(2 - iv[r+1][2])), K = 10^80);
    blo[r] = sqrtint(floor(L*K))/10^40;
    bhi[r] = (sqrtint(ceil(H*K)) + 1)/10^40;
    if(blo[r]^2 > L || bhi[r]^2 < H, error("square-root rounding")));
  [blo, bhi];
}

\\ (1) subtorus relations and prefix independence in Q(zeta_2N)
relations_ok(N, C) =
{
  my(m = (N-1)/2, M = eulerphi(2*N)/2, cyc = polcyclo(2*N), z = Mod(x, cyc), w);
  if(matsize(C) != [m, M], return(0));
  w = vector(m, r, z^r - z^(-r));
  if(matrank(matrix(M, poldegree(cyc), s, j, polcoeff(lift(w[s]), j-1))) != M, return(0));
  for(r = 1, m, if(sum(s = 1, M, C[r,s]*w[s]) != w[r], return(0)));
  1;
}

\\ (3) positive definiteness by fraction-free Gaussian elimination: pivot k is the leading minor Delta_k
pd_ok(A) =
{
  my(n = #A, W = A, prev = 1, p);
  for(k = 1, n,
    p = W[k,k];
    if(imag(p) != 0 || p <= 0, return(0));
    for(i = k+1, n, for(j = k+1, n, W[i,j] = (p*W[i,j] - W[i,k]*W[k,j])/prev));
    prev = p);
  1;
}

\\ certify(cert, sgn): sgn = 1 is the real objective; sgn = -1 flips sin -> -sin (a mutation).
\\ Returns [verdict, delta, reason].
certify(cert, sgn) =
{
  my(N = cert[1], C = cert[2], B = cert[3], den = cert[4], tri = cert[5], claim = cert[6]);
  my(n = matsize(B)[1], m = (N-1)/2, Qi = matrix(n, n), t = 0, S = Map(), T = Map(), keys, v, rd, l1 = 0, bb, bt, c, delta);
  if(N % 2 == 0, error("odd N only"));
  if(!relations_ok(N, C), return([0, 0, "subtorus relation fails"]));
  for(i = 1, n, for(j = 1, i, t++; Qi[i,j] = tri[t]; Qi[j,i] = conj(tri[t])));
  for(i = 1, n, if(imag(Qi[i,i]) != 0, return([0, 0, "Q not Hermitian"])));
  if(!pd_ok(Qi), return([0, 0, "Q not positive definite"]));
  bb = benclose(N);
  bt = vector(m, r, (bb[1][r] + bb[2][r])/2);
  for(a = 1, n, for(k = 1, n,
    my(d = B[k,] - B[a,]);
    mapput(S, d, if(mapisdefined(S, d, &v), v, 0) + Qi[a,k]/den)));
  for(r = 1, m,
    my(p = C[r,], q = -C[r,]);
    mapput(T, p, if(mapisdefined(T, p, &v), v, 0) + sgn*I*bt[r]/2);
    mapput(T, q, if(mapisdefined(T, q, &v), v, 0) - sgn*I*bt[r]/2));
  c = mapget(S, vector(matsize(B)[2]));
  keys = setunion(Set(Mat(S)[,1]~), Set(Mat(T)[,1]~));
  for(e = 1, #keys,
    my(d = keys[e]);
    if(d == vector(#d), next);
    rd = if(mapisdefined(T, d, &v), v, 0) - if(mapisdefined(S, d, &v), v, 0);
    l1 += abs(real(rd)) + abs(imag(rd)));
  delta = vecsum(bb[1]) - c - l1 - vecsum(bb[2] - bb[1]);
  [delta >= claim, delta, if(delta >= claim, "proved", "bound below the claim"), l1];
}

npass = 0;
print("Exact SOS certificates for U_N - A_N (claims proved in exact arithmetic):");
{
  for(i = 1, #defect_certs,
    my(cert = defect_certs[i], res = certify(cert, 1));
    printf("  N = %2d  |B| = %2d  ||r||_1 <= %.1e  delta >= %.10f  claim %.8f : %s\n", cert[1],
           matsize(cert[3])[1], if(#res > 3, res[4], 0), floor(res[2]*10^10)/10^10, cert[6],
           if(res[1], "PASS", Str("FAIL (", res[3], ")")));
    npass += res[1]);
}
\\ a data file that failed to load leaves defect_certs a bare variable: count, do not trust silence
ok = type(defect_certs) == "t_VEC" && #defect_certs == 3 && npass == 3;

\\ Mutations: each must be rejected.  Raise every claim by 10^-7; perturb Q, C and the objective at N = 9.
print("Mutations (each must be rejected):");
nmut = 0; nrej = 0;
mutate(label, cert, sgn) =
{
  my(res = certify(cert, sgn));
  nmut++; if(!res[1], nrej++);
  printf("  %-44s: %s\n", label, if(res[1], "ACCEPTED (bad)", Str("rejected (", res[3], ")")));
}
{
  for(i = 1, #defect_certs,
    my(cert = defect_certs[i]);
    cert[6] += 1/10^7;
    mutate(Str("N = ", cert[1], ": claim + 1e-7"), cert, 1));
  my(c9 = defect_certs[1], cert, tri, j);
  cert = c9; tri = cert[5]; tri[#tri] = -tri[#tri]; cert[5] = tri;
  mutate("N = 9: last diagonal entry of Q negated", cert, 1);
  cert = c9; tri = cert[5]; j = 2; tri[j] = tri[j]*1001/1000; cert[5] = tri;
  mutate("N = 9: Q[2,1] scaled by 1.001", cert, 1);
  cert = c9; cert[2][matsize(cert[2])[1], 1] += 1;
  mutate("N = 9: C[m,1] + 1 (false relation)", cert, 1);
  mutate("N = 9: objective sin -> -sin", c9, -1);
}
okm = nmut == 7 && nrej == nmut;
print("  rejected ", nrej, " of ", nmut, " mutations: ", if(okm, "PASS", "FAIL"));
print("RESULT: ", if(ok && okm, "PASS", "FAIL"), " -- defect lower bounds certified for N = 9, 15, 21");
quit;
