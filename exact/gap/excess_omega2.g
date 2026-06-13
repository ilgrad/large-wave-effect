# Exact relation-lattice structure for the OPEN case omega(m) = 2 (odd part m = p*q, two distinct
# primes) of the excess lemma  A_N <= L_pre + O(1)  (Remark rem:defect of paper/sections/linear.tex).
#
# The excess lemma is PROVED when omega(m) <= 1 (odd part a prime power): N=2p has a single dependent
# mode omega_p=2 of weight 1/(2N), and N=4p two dependent modes; the dependent band is O(1) so the
# unconditional sandwich |A_N - L_pre| <= U_N - L_pre = (1/pi)ln(m/phi(m)) + O(1) already closes it.
# This script characterizes, in EXACT cyclotomic arithmetic, why omega(m)=2 is structurally different:
# the dependent band now carries Theta(N) modes, and the relation lattice has a rich chained structure.
#
# Findings (reproduced below):
#   (A) rank Lambda_N = floor(N/2) - phi(2N)/2, and the lattice is UNIMODULAR onto its span: the Smith
#       normal form is all 1's (no torsion).  So the obstruction is NOT a torsion/index phenomenon.
#   (B) A minimal generating set consists of "chained" relations with alternating coefficients
#       (+1,+1,-1,-1,+1,...).  The shortest support is governed by the SMALLEST prime factor:
#         * 3 | N  ->  support-3 generators  omega_a + omega_b = omega_c  with  a + b = N/3,
#           because 2 sin(pi(a+b)/2N) = 2 sin(pi/6) = 1 collapses the product-to-sum
#           sin(pi a/N)+sin(pi b/N) = 2 sin(pi(a+b)/2N) cos(pi(a-b)/2N) into a single sine.
#         * 5 | N (and 3 not | N)  ->  support-5 generators (the golden-ratio angle 2 sin(pi/10)).
#         * the second prime supplies a long "comb" generator of support ~ 2*(smaller prime)-1.
#   (C) the SHORT chained generators all have  sum_r k_r  in {1, 3} mod 4 (never 0); some long "comb"
#       generators do hit sum = 0 mod 4, but since saturation needs ALL generators = 0 mod 4, the short
#       ones alone violate the mod-4 ceiling criterion (Theorem thm:ceiling), so A_N < U_N -- these N
#       are non-saturating, as required.
#
# The "prime power" step of the omega(m)<=1 proof is exactly (A)+(B) degenerating: for a prime power
# odd part the dependent band is O(1) generators of O(1) total weight, the sandwich is O(1/N), and the
# excess lemma is trivial.  For omega(m)=2 the SAME sandwich is (1/pi)ln(m/phi(m)) = Theta(1) and grows
# with further prime factors, so the lemma stops being free: it becomes a genuine SIGNED-cancellation
# statement on the relation-locked dependent band (the open core; see numerics/verify_excess_omega2.py).
#
#   gap -q --nointeract exact/gap/excess_omega2.g

Read("exact/gap/relation_lattice.g");

LW_OddPart := function(N)
  local m;
  m := N;
  while m mod 2 = 0 do m := m / 2; od;
  return m;
end;;

LW_Omega := N -> Length(Filtered(PrimeDivisors(LW_OddPart(N)), p -> true));;

# Exact cyclotomic residual of a putative relation sum_i co[i] * omega_{idx[i]} = 0; returns 0 iff it
# holds (omega_r and zeta^r - zeta^-r are proportional, so this is the genuine real relation).
LW_RelResidual := function(N, idx, co)
  local z, s, i;
  z := E(2 * N); s := 0;
  for i in [1 .. Length(idx)] do s := s + co[i] * (z ^ idx[i] - z ^ (-idx[i])); od;
  return s;
end;;

# Smith normal form diagonal (nonzero invariant factors) of the relation-lattice basis.
LW_SNFDiag := function(N)
  local B, snf, i, d;
  B := LW_RelationLattice(N);
  if Length(B) = 0 then return []; fi;
  snf := SmithNormalFormIntegerMat(B);
  d := [];
  for i in [1 .. Minimum(Length(snf), Length(snf[1]))] do
    if snf[i][i] <> 0 then Add(d, snf[i][i]); fi;
  od;
  return d;
end;;

# Short generators (LLL-reduced), sorted by support size, with their support / coeffs / sum mod 4.
LW_ShortGens := function(N)
  local B;
  B := LLLReducedBasis(LW_RelationLattice(N)).basis;
  Sort(B, function(a, b) return Number(a, x -> x <> 0) < Number(b, x -> x <> 0); end);
  return B;
end;;

# ---- (A) ranks and Smith normal form -----------------------------------------------------------
Print("==== omega(m)=2 relation lattice: rank + Smith normal form (unimodular onto span) ====\n");
Print("  N : m=floor(N/2)  M=phi(2N)/2  rank(Lambda)  m-M  SNF-diag\n");
LW_OK := true;;
for N in [15, 21, 33, 35, 45, 55, 30, 60, 42, 66, 70, 90] do
  if LW_Omega(N) <> 2 then continue; fi;
  B := LW_RelationLattice(N);;
  m := Int(N / 2);; M := Phi(2 * N) / 2;;
  snf := LW_SNFDiag(N);;
  if Length(B) <> m - M then LW_OK := false; fi;
  if Length(snf) > 0 and Set(snf) <> [1] then LW_OK := false; fi;
  Print("  ", String(N, 3), " : ", String(m, 11), String(M, 12),
        String(Length(B), 13), String(m - M, 6), "   ", snf, "\n");
od;
Print("  rank = m - phi(2N)/2 and SNF all-ones for every case: ", LW_OK, "\n\n");

# ---- (B) shortest-support generators and their arithmetic origin -------------------------------
Print("==== shortest generators: support set by smallest prime; chained alternating coeffs ====\n");
for N in [15, 21, 33, 45, 35, 55] do
  if LW_Omega(N) <> 2 then continue; fi;
  gens := LW_ShortGens(N);;
  pmin := Minimum(PrimeDivisors(LW_OddPart(N)));;
  shortsupp := Minimum(List(gens, k -> Number(k, x -> x <> 0)));;
  Print("  N=", String(N, 3), "  primes=", PrimeDivisors(LW_OddPart(N)),
        "  smallest p=", pmin, "  shortest support=", shortsupp, "\n");
  for k in gens do
    supp := Filtered([1 .. Length(k)], r -> k[r] <> 0);;
    if Length(supp) = shortsupp then
      Print("    idx=", supp, " co=", List(supp, r -> k[r]),
            " sum_mod4=", Sum(k) mod 4);
      if shortsupp = 3 then
        Print("   [a+b=", supp[1] + supp[2], "=N/3? ", supp[1] + supp[2] = N / 3, "]");
      fi;
      Print("\n");
    fi;
  od;
od;

# ---- verify the support-3 mechanism exactly (the factor-3 collapse 2 sin(pi/6) = 1) -------------
Print("\n==== exact verification: support-3 generators are genuine sine identities (residual = 0) ==\n");
LW_ID := true;;
for trip in [[15, [1, 4, 6]], [15, [2, 3, 7]], [21, [1, 6, 8]], [21, [3, 4, 10]],
             [33, [1, 10, 12]], [45, [1, 14, 16]]] do
  res := LW_RelResidual(trip[1], trip[2], [1, 1, -1]);;
  if res <> 0 then LW_ID := false; fi;
  Print("  N=", String(trip[1], 3), " omega_", trip[2][1], "+omega_", trip[2][2],
        "=omega_", trip[2][3], "  (a+b=", trip[2][1] + trip[2][2], "=N/3)  residual=", res, "\n");
od;
Print("  all support-3 identities exact: ", LW_ID, "\n");

# verify support-5 identities for 5*q with 3 not dividing N (golden-ratio angle)
Print("\n==== exact verification: support-5 generators for 5|N, 3 not | N (residual = 0) ====\n");
LW_ID5 := true;;
for quint in [[35, [1, 6, 8, 13, 15]], [35, [2, 5, 9, 12, 16]], [55, [1, 10, 12, 21, 23]]] do
  res := LW_RelResidual(quint[1], quint[2], [1, 1, -1, -1, 1]);;
  if res <> 0 then LW_ID5 := false; fi;
  Print("  N=", String(quint[1], 3), " idx=", quint[2], " co=[1,1,-1,-1,1]  residual=", res, "\n");
od;
Print("  all support-5 identities exact: ", LW_ID5, "\n");

# ---- (C) mod-4 verdict: a short generator has nonzero coeff sum mod 4 => non-saturating ---------
# (Saturation needs ALL generators = 0 mod 4; the short chained ones violate it.  The long comb
# generators MAY have sum = 0 mod 4 -- that does not rescue saturation.)
Print("\n==== mod-4 ceiling verdict: some generator has sum != 0 mod 4 => A_N < U_N (non-saturating) =\n");
LW_NS := true;;
for N in [15, 21, 33, 35, 45, 55, 30, 60, 42, 66, 70, 90] do
  if LW_Omega(N) <> 2 then continue; fi;
  crit := LW_Criterion(N);;
  if crit.saturates then LW_NS := false; fi;
  Print("  N=", String(N, 3), "  saturates=", crit.saturates,
        "  (node-0 generator sums mod4 = ", crit.summ0, ")\n");
od;
Print("  every omega(m)=2 case is non-saturating (A_N < U_N): ", LW_NS, "\n");

Print("\n==== RESULT ====\n");
if LW_OK and LW_ID and LW_ID5 and LW_NS then
  Print("PASS: rank=m-phi(2N)/2, SNF unimodular, chained generators exact, all non-saturating.\n");
  Print("The excess lemma's 'prime power' step = the dependent band being O(1) generators of O(1)\n");
  Print("weight; at omega(m)=2 it is Theta(N) modes of weight (1/pi)ln(m/phi(m)) = Theta(1), so the\n");
  Print("free sandwich no longer closes the lemma -- it becomes a signed-cancellation statement.\n");
  QUIT_GAP(0);
else
  Print("FAIL\n");
  QUIT_GAP(1);
fi;
