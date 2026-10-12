\\ Exact layer for B_p = sqrt(p) on prime rings (Proposition prop:bn-prime of the paper).
\\
\\   gp -q -f exact/pari/bn_prime_flat.gp </dev/null
\\
\\ B_N = sup_t sum_k |K_N(k,t)|, K_N(k,t) = (1/N) sum_r exp(-i t lambda_r) zeta_N^(r k),
\\ lambda_r = 2 - zeta_N^r - zeta_N^-r.  For an odd prime p the proof has two exact inputs:
\\   (1) lambda_1..lambda_m, m = (p-1)/2, are Q-linearly independent (so Kronecker aligns the phases);
\\   (2) at the aligned phases exp(-i t lambda_r) = zeta_p^(r^2) the kernel is a Gauss sum with
\\       |sum_r zeta_p^(r^2 + r k)|^2 = p for every k, i.e. |K_p(k)| = 1/sqrt(p) and sum_k |K_p(k)| = sqrt(p).
\\ Both are checked here in exact cyclotomic arithmetic.  Part (3) is numerical: explicit times t with
\\ ||K_p(.,t)||_1 / sqrt(p) within 1e-5 of 1 for p = 5, 7, and the attained value 1 at p = 3, t = 2 pi/9.

\\ (1) Q-rank of {lambda_r} in exact coordinates of Q(zeta_p) (power basis modulo Phi_p)
lamrank(p) =
{
  my(c = polcyclo(p), z = Mod(x, c));
  matrank(matrix((p-1)/2, p-1, r, j, polcoeff(lift(2 - z^r - z^(p-r)), j-1)));
}

\\ (2) |G_k|^2 == p exactly for all k, G_k = sum_r zeta^(r^2 + r k); conj(zeta^e) = zeta^(-e)
gaussflat(p) =
{
  my(c = polcyclo(p), z = Mod(x, c), ok = 1);
  for(k = 0, p-1,
    my(g = sum(r = 0, p-1, z^((r^2 + r*k) % p)), gb = sum(r = 0, p-1, z^((-(r^2 + r*k)) % p)));
    if(lift(g * gb) != p, ok = 0));
  ok;
}

\\ np counts the primes actually checked, so a loop that dies early cannot report PASS
ok1 = 1; ok2 = 1; np = 0;
{
  forprime(p = 3, 53,
    if(lamrank(p) != (p-1)/2, ok1 = 0; print("  rank deficit at p = ", p));
    if(!gaussflat(p), ok2 = 0; print("  Gauss sum not flat at p = ", p));
    np++);
}
ok1 = ok1 && np == 15; ok2 = ok2 && np == 15;
print("(1) rank_Q {lambda_r} = (p-1)/2, ", np, " odd primes p <= 53: ", if(ok1, "PASS", "FAIL"));
print("(2) |sum_r zeta^(r^2+rk)|^2 = p exactly, all k, p <= 53   : ", if(ok2, "PASS", "FAIL"));

\\ (3) numerical witnesses.  ratio(p,t) = ||K_p(.,t)||_1 / sqrt(p).
default(realprecision, 50);
ratio(p, t) =
{
  my(lam = vector(p, r, 2 - 2*cos(2*Pi*(r-1)/p)));
  sum(k = 0, p-1, abs(sum(r = 0, p-1, exp(-I*t*lam[r+1] + 2*Pi*I*r*k/p)))) / p / sqrt(p);
}
\\ Kronecker search: fix the r = 1 phase exactly on the target, t = 2 pi (n - 1/p)/lambda_1, and score the
\\ remaining targets -t lambda_r = 2 pi r^2/p (mod 2 pi) by sum (1 - cos(error)); then refine t locally.
witness(p, nmax) =
{
  my(m = (p-1)/2, lam = vector(m, r, 2 - 2*cos(2*Pi*r/p)), be = oo, bt = 0.);
  for(n = 1, nmax,
    my(t = 2*Pi*(n - 1/p)/lam[1], e = 0.);
    for(r = 2, m, e += 1 - cos(t*lam[r] + 2*Pi*r^2/p));
    if(e < be, be = e; bt = t));
  my(h = 1e-3, best = ratio(p, bt));
  for(it = 1, 40,
    my(up = ratio(p, bt + h), dn = ratio(p, bt - h));
    if(up > best, bt += h; best = up, if(dn > best, bt -= h; best = dn, h /= 2)));
  [bt, best];
}
w3 = ratio(3, 2*Pi/9);
w5 = witness(5, 20000);
w7 = witness(7, 400000);
printf("(3) p = 3: ratio at t = 2 pi/9          = %.30f\n", w3);
printf("    p = 5: ratio at t = %.10f = %.12f\n", w5[1], w5[2]);
printf("    p = 7: ratio at t = %.10f = %.12f\n", w7[1], w7[2]);
ok3 = abs(w3 - 1) < 1e-40 && w5[2] > 1 - 1e-5 && w7[2] > 1 - 1e-5 && w5[2] <= 1 && w7[2] <= 1;
print("    B_3 attained; B_5/sqrt5, B_7/sqrt7 within 1e-5 of 1      : ", if(ok3, "PASS", "FAIL"));
print("RESULT: ", if(ok1 && ok2 && ok3, "PASS", "FAIL"), " -- B_p = sqrt(p) inputs verified");
quit;
