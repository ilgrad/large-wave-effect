\\ Exact check of Lemma cleancomb: for 3 | N the support-3 relations omega_a + omega_b = omega_c
\\ (omega_r = 2 sin(pi r/N), 1 <= a < b < c <= N/2) are exactly a + b = N/3, c = N/3 + a, plus, when 30 | N,
\\ the single extra triple (N/10, N/6, 3N/10) coming from the minimal relation (R5:R3).
\\
\\ A 60-digit scan with tolerance 1e-40 finds every candidate (a true relation evaluates to ~1e-58);
\\ each candidate is then decided exactly: the relation is the vanishing of the cyclotomic integer
\\ z^a + z^(N-a) + z^b + z^(N-b) + z^(N+c) + z^(2N-c), z = zeta_{2N}, tested modulo Phi_{2N}.
\\
\\   gp -q exact/pari/support3_relations.gp

default(realprecision, 60);

isrel(N, a, b, c) = {
  my(P = polcyclo(2*N));
  (x^a + x^(N-a) + x^b + x^(N-b) + x^(N+c) + x^(2*N-c)) % P == 0;
}

NMAX = 330;
nstd = 0; extra = List(); spurious = 0;
{
forstep(N = 3, NMAX, 3,
  my(s = vector(N\2, r, sin(Pi*r/N)));
  for(a = 1, N\2, for(b = a+1, N\2,
    my(v = s[a] + s[b]);
    if(v > 1 + 1e-30, break);
    for(c = b+1, N\2,
      if(abs(s[c] - v) < 1e-40,
        if(!isrel(N, a, b, c), spurious++; next);
        if(a + b == N/3 && c == N/3 + a, nstd++, listput(extra, [N, a, b, c])))))));
}
want = vector(NMAX\30, k, [30*k, 3*k, 5*k, 9*k]);
ok1 = (Vec(extra) == want) && (spurious == 0);
printf("[1] 3 | N <= %d: %d standard solutions, extra = %s, spurious candidates = %d\n", NMAX, nstd, Vec(extra), spurious);
printf("    extras are exactly (N/10, N/6, 3N/10) for 30 | N: %s\n", if(ok1, "yes", "NO"));

\\ every standard solution is exact (spot-check the identity sin(pi/3 + x) = sin x + sin(pi/3 - x) in Q(zeta_2N))
ok2 = 1;
{
forstep(N = 3, 120, 3, for(a = 1, ceil(N/6) - 1,
  if(!isrel(N, a, N/3 - a, N/3 + a), ok2 = 0)));
}
printf("[2] standard triplets (a, N/3-a, N/3+a) vanish exactly, 3 | N <= 120: %s\n", if(ok2, "yes", "NO"));

print("RESULT: ", if(ok1 && ok2, "PASS", "FAIL"));
quit
