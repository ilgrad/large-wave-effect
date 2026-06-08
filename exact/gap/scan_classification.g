# Independent exact verification (GAP) of the linear-theory arithmetic:
#   (1) rank_Q {2 sin(pi r/N) : 1<=r<=N-1} = 1/2 phi(2N)            (Theorem qrank)
#   (2) A_N = U_N  <=>  N prime or N = 2^m                          (Theorem ceiling + scan)
#   (3) every composite N that is not a power of two carries a mod-4 obstruction.
#
# Reproduces, in exact cyclotomic arithmetic and independently of the Python core,
# the saturation scan -- extended here well past the N<=160 of the paper.
#
#   gap -q --nointeract exact/gap/scan_classification.g
# Set LW_MAXN below (default 400; the paper reports the verified range N<=800, ~35 s).

Read("exact/gap/relation_lattice.g");
if not IsBound(LW_MAXN) then LW_MAXN := 400; fi;

LW_rankok := true;; LW_sat := [];; LW_t0 := Runtime();;
for N in [3 .. LW_MAXN] do
  if LW_QRankFull(N) <> Phi(2 * N) / 2 then
    LW_rankok := false; Print("RANK MISMATCH at N = ", N, "\n");
  fi;
  if LW_Criterion(N).saturates then Add(LW_sat, N); fi;
od;

Print("==== exact GAP verification, N = 3 .. ", LW_MAXN, " ====\n");
Print("rank_Q = 1/2 phi(2N) for all N        : ", LW_rankok, "\n");
Print("saturating set = {prime} U {2^m}      : ",
      LW_sat = Filtered([3 .. LW_MAXN], LW_IsPrimeOrTwoPower),
      "   (", Length(LW_sat), " values)\n");
Print("elapsed                               : ", (Runtime() - LW_t0) / 1000.0, " s\n");
if LW_rankok and LW_sat = Filtered([3 .. LW_MAXN], LW_IsPrimeOrTwoPower) then
  Print("RESULT: PASS\n"); QUIT_GAP(0);
else
  Print("RESULT: FAIL\n"); QUIT_GAP(1);
fi;
