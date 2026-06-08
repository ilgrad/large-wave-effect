# Exact cyclotomic core for the ring frequencies omega_r = 2 sin(pi r / N).
#
# We never divide by i: the relations among {omega_r} coincide with the relations among
# {zeta^r - zeta^{-r}} (zeta = zeta_{2N}), since multiplying by the nonzero scalar 1/(2i)
# preserves Q-linear dependence. This mirrors src/large_wave_effect/cyclotomic.py exactly,
# but in GAP's exact cyclotomic arithmetic -- an independent reimplementation.
#
# Read this file from the repository root:  Read("exact/gap/relation_lattice.g");

# Integer coordinates of (zeta^r - zeta^{-r}) in Q(zeta_{2N}), for the distinct positive
# frequencies r = 1 .. floor(N/2).
LW_FreqCoords := function(N)
  local z;
  z := E(2 * N);
  return List([1 .. Int(N / 2)], r -> CoeffsCyc(z^r - z^(-r), 2 * N));
end;;

# Q-rank of the full frequency set {2 sin(pi r/N) : r = 1 .. N-1} (Theorem: = 1/2 phi(2N)).
LW_QRankFull := function(N)
  local z;
  z := E(2 * N);
  return RankMat(List([1 .. N - 1], r -> CoeffsCyc(z^r - z^(-r), 2 * N)));
end;;

# Integer basis of the relation lattice Lambda_N = { k in Z^m : sum_r k_r omega_r = 0 }.
LW_RelationLattice := function(N)
  return NullspaceIntMat(LW_FreqCoords(N));
end;;

# The mod-4 ceiling criterion (Theorem ceiling). Returns a record with both branches and the
# saturation verdict A_N = U_N. node-0 branch: sum_r k_r = 0 (mod 4); antipodal (even N, node
# N/2): sum_r k_r + 2 sum_r r k_r = 0 (mod 4). Checking a generating set suffices because each
# branch is a homomorphism Lambda_N -> Z/4.
LW_Criterion := function(N)
  local B, k, r, s0, sh, sat0, sath;
  B := LW_RelationLattice(N);
  sat0 := true; sath := true; s0 := []; sh := [];
  for k in B do
    Add(s0, Sum(k) mod 4);
    Add(sh, (Sum(k) + 2 * Sum([1 .. Length(k)], r -> r * k[r])) mod 4);
    if Sum(k) mod 4 <> 0 then sat0 := false; fi;
    if (Sum(k) + 2 * Sum([1 .. Length(k)], r -> r * k[r])) mod 4 <> 0 then sath := false; fi;
  od;
  return rec(relrank := Length(B), summ0 := s0, summA := sh,
             saturates := sat0 or (N mod 2 = 0 and sath));
end;;

LW_IsPrimeOrTwoPower := N -> IsPrimeInt(N) or (N = 2 ^ LogInt(N, 2));;
