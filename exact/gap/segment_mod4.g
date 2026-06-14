# Segment ceiling criterion (the Dirichlet analogue of the ring mod-4 criterion, Theorem ceiling).
#
# The Dirichlet segment tridiag(-1,2,-1) has the SIMPLE spectrum lambda_k^Dir = 4 sin^2(pi k/(2(N+1))),
# k = 1..N, frequencies mu_k = 2 sin(pi k/(2(N+1))).  In cyclotomic coordinates 2i mu_k = z^k - z^{-k},
# z = zeta_M, M = 4(N+1).  The relation lattice is Lambda^Dir = { k in Z^N : sum_r k_r mu_r = 0 }, computed
# exactly by NullspaceIntMat of the cyclotomic coordinate matrix.
#
# SEGMENT CEILING CRITERION:  A_N^Dir = U_N^Dir  iff  sum_r k_r = 0 (mod 4) for every k in Lambda^Dir.
# (Same Kronecker reachability argument as the ring: U_N^Dir is realised iff the unique maximiser
# phi_* = (pi/2,...,pi/2) lies in the orbit closure, i.e. <k,phi_*> = (pi/2) sum_r k_r = 0 (mod 2pi),
# i.e. sum_r k_r = 0 (mod 4), for every relation k.  The simple spectrum has NO antipodal node, unlike the
# ring's even-N branch.)
#
# This script CERTIFIES that the mod-4 criterion holds exactly for N+1 prime or a power of two, i.e. the
# criterion bridges the rank deficiency to A_N^Dir < U_N^Dir for composite N+1 (the step thm:segrank states).
# Run:  gap -q --nointeract exact/gap/segment_mod4.g

LW_SegFreqCoords := function(N)
  local M, z;
  M := 4 * (N + 1); z := E(M);
  return List([1 .. N], k -> CoeffsCyc(z^k - z^(-k), M));
end;;

# Lambda^Dir and the criterion verdict (true = "every relation has Sum(k) = 0 mod 4").
LW_SegCriterion := function(N)
  local B, k, sat, viol;
  B := NullspaceIntMat(LW_SegFreqCoords(N));
  sat := true; viol := fail;
  for k in B do
    if Sum(k) mod 4 <> 0 then sat := false; if viol = fail then viol := k; fi; fi;
  od;
  return rec(relrank := Length(B), saturates := sat, witness := viol);
end;;

LW_IsPrimeOrTwoPower := M -> IsPrimeInt(M) or (M = 2 ^ LogInt(M, 2));;

allok := true;
Print("N   N+1  rk(Lambda^Dir)  mod4-criterion   N+1 prime/2^m   match\n");
for N in [1 .. 40] do
  c := LW_SegCriterion(N);
  exp := LW_IsPrimeOrTwoPower(N + 1);
  m := (c.saturates = exp);
  allok := allok and m;
  Print(String(N, 2), " ", String(N + 1, 4), " ", String(c.relrank, 12), "  ",
        String(c.saturates, 10), "  ", String(exp, 14), "   ", m, "\n");
od;
Print("\nRESULT: ", allok,
      "  -- segment mod-4 criterion <=> (N+1 prime or 2^m) for N<=40; so for composite N+1 the\n");
Print("        relation lattice carries a relation with Sum(k) != 0 (mod 4), forcing A_N^Dir < U_N^Dir.\n");
QUIT_GAP(0);
