# Independent exact verification (GAP) of the Schrodinger orbit-closure arithmetic
# (Theorem bn-saturation), with lambda_r = 2 - z^r - z^-r, z = E(N):
#   (1) rank_Q {lambda_r : 1 <= r <= N/2} = 1/2 phi(N)  for 3 <= N <= 150   (part 1)
#   (2) for N = p^a the integer relation lattice of the lambda_r, folded onto r = 1..N/2, is spanned by
#       the coset relations prod_{coset j} x = prod_{coset 0} x of the subgroup of order p   (part 2)
#
# A second exact implementation next to exact/pari/bn_saturation.gp.
#
#   gap -q -A --quitonbreak exact/gap/schrodinger_saturation.g

LambdaCoords := function(N)
  local B, rows, r;
  B := CanonicalBasis(CF(N));
  rows := [];
  for r in [1 .. Int(N/2)] do
    Add(rows, Coefficients(B, 2 - E(N)^r - E(N)^(-r)));
  od;
  return rows;
end;;

FoldedCosetRels := function(N, p)
  local q, rels, j, v, r, s;
  q := N/p;
  rels := [];
  for j in [1 .. q-1] do
    v := List([1 .. Int(N/2)], i -> 0);
    for r in [1 .. N-1] do
      s := Minimum(r, N-r);
      if r mod q = j then v[s] := v[s] + 1; fi;
      if r mod q = 0 then v[s] := v[s] - 1; fi;
    od;
    Add(rels, v);
  od;
  return rels;
end;;

# HermiteNormalFormIntegerMat keeps zero rows, so compare the nonzero parts
NonzeroHNF := M -> Filtered(HermiteNormalFormIntegerMat(M), r -> not IsZero(r));;

rankok := ForAll([3 .. 150], N -> RankMat(LambdaCoords(N)) = Phi(N)/2);;
Print("[gap] rank_Q{lambda_r} = phi(N)/2 for 3 <= N <= 150: ", rankok, "\n");

latok := ForAll([[9, 3], [25, 5], [27, 3], [49, 7], [81, 3], [121, 11], [125, 5]],
  c -> NonzeroHNF(NullspaceIntMat(LambdaCoords(c[1]))) = NonzeroHNF(FoldedCosetRels(c[1], c[2])));;
Print("[gap] N = p^a in {9,25,27,49,81,121,125}: relation lattice = coset relations: ", latok, "\n");

if rankok and latok then
  Print("RESULT: PASS\n"); QUIT_GAP(0);
else
  Print("RESULT: FAIL\n"); QUIT_GAP(1);
fi;
