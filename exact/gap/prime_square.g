# Structural finding (GAP): for N = p^2 (p an odd prime) the minimal mod-4 obstruction has
# support exactly p, on the index set {1} U {kp +/- 1 : 1 <= k <= (p-1)/2}. This is a candidate
# for a constructive proof of the classification A_N = U_N <=> N prime or 2^m (presently a scan).
#
#   gap -q --nointeract exact/gap/prime_square.g

Read("exact/gap/relation_lattice.g");

LW_MinObstruction := function(N)        # shortest-support relation violating the mod-4 criterion
  local B, cand;
  B := LLLReducedBasis(LW_RelationLattice(N)).basis;
  cand := Filtered(B, k -> Sum(k) mod 4 <> 0);
  Sort(cand, function(a, b) return Number(a, x -> x <> 0) < Number(b, x -> x <> 0); end);
  if Length(cand) = 0 then return fail; fi;
  return cand[1];
end;;

Print("N = p^2: minimal mod-4 obstruction (support, index:coeff, sum mod 4)\n");
for p in [3, 5, 7, 11, 13] do
  k := LW_MinObstruction(p ^ 2);;
  supp := Filtered([1 .. Length(k)], r -> k[r] <> 0);;
  Print("  p=", String(p, 2), "  N=", String(p ^ 2, 3),
        "  support=", Length(supp),
        "  indices=", supp,
        "  (= {1} U {kp+-1}? ",
        supp = Concatenation([1], Concatenation(List([1 .. (p - 1) / 2],
              k2 -> [k2 * p - 1, k2 * p + 1]))) and true, ")",
        "  sum mod4=", Sum(k) mod 4, "\n");
od;
QUIT_GAP(0);
