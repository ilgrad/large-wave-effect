(* Machine-checked arithmetic core of the saturation classification (Theorem classification),
   uniform in N -- not a finite scan.  Self-contained: Rocq 9.x standard library only.

   For composite N = p*s that is not a power of two, take an odd prime factor p (so s = N/p >= 2).
   The imaginary part of  sum_{j<p} zeta_{2N}^{1 + 2 s j} = 0  is a relation in Lambda_N whose p
   reduced indices e_j = (1 + 2 s j) mod 2N avoid {0, N}; every reduction sign is therefore +-1, the
   coefficient sum S is a sum of p (odd many) signs, hence odd, and an odd S violates BOTH mod-4
   branches of the ceiling criterion, giving A_N < U_N.

   This file proves, for ALL N, the number theory behind "indices avoid {0,N}" (so signs are +-1)
   and the combinatorial fact "odd-length +-1 sum blocks both branches".  What is supplied from
   outside (proved on paper, checked exactly by GAP/PARI and numerics/verify_classification_proof.py)
   is only that the constructed vector is a genuine frequency relation -- the cyclotomic input. *)

From Stdlib Require Import ZArith Lia List.
Import ListNotations.
Open Scope Z_scope.

(* (1) 1 + 2 s j is odd, while 2N is even: no index ever reduces to 0. *)
Lemma index_odd : forall s j, Z.odd (1 + 2*s*j) = true.
Proof. intros. rewrite Z.odd_add, Z.odd_mul, Z.odd_mul. now destruct (Z.odd s), (Z.odd j). Qed.

(* (2) The crux: (p s - 1) mod (2 s) = s - 1 for p odd, s > 0.  Nonzero for s >= 2, so 2s does not
   divide N - 1 = p s - 1, i.e. no index reduces to N = p s. *)
Lemma ps_minus1_mod : forall p s, Z.odd p = true -> 0 < s -> (p*s - 1) mod (2*s) = s - 1.
Proof.
  intros p s Hp Hs. apply Z.odd_spec in Hp. destruct Hp as [m Hm].
  replace (p*s - 1) with ((s - 1) + m*(2*s)) by lia.
  rewrite Z.mod_add by lia. apply Z.mod_small. lia.
Qed.

Lemma index_ne_N : forall p s, Z.odd p = true -> 2 <= s -> (p*s - 1) mod (2*s) <> 0.
Proof. intros p s Hp Hs. rewrite ps_minus1_mod by (assumption || lia). lia. Qed.

(* (3) An odd value violates the mod-4 criterion on both branches. *)
Lemma mod4_ne0_of_odd : forall x, Z.odd x = true -> x mod 4 <> 0.
Proof.
  intros x Hodd Hmod. apply Z.mod_divide in Hmod; [| lia]. destruct Hmod as [c Hc].
  assert (Heven : Z.even x = true) by (apply Z.even_spec; exists (2*c); lia).
  rewrite <- Z.negb_odd, Hodd in Heven. discriminate.
Qed.

Lemma odd_blocks_both : forall S T, Z.odd S = true -> S mod 4 <> 0 /\ (S + 2*T) mod 4 <> 0.
Proof.
  intros S T HS. split.
  - now apply mod4_ne0_of_odd.
  - apply mod4_ne0_of_odd. now rewrite Z.odd_add, Z.odd_mul, HS.
Qed.

(* (4) The coefficient sum: a sum of +-1 signs has the parity of its length. *)
Fixpoint zsum (l : list Z) : Z := match l with [] => 0 | x :: t => x + zsum t end.

Lemma zsum_parity : forall l, (forall x, In x l -> Z.odd x = true) ->
  Z.odd (zsum l) = Nat.odd (length l).
Proof.
  induction l as [| a l IH]; intros H.
  - reflexivity.
  - simpl. rewrite Z.odd_add, (H a (or_introl eq_refl)), IH
      by (intros x Hx; apply H; now right).
    now rewrite Nat.odd_succ, <- Nat.negb_odd; destruct (Nat.odd (length l)).
Qed.

(* CAPSTONE (uniform in N): a single relation whose signs are all +-1 and whose support count is odd
   blocks BOTH mod-4 branches of the ceiling criterion -- hence A_N < U_N.  With index_ne_N (the
   construction's reduced indices avoid N, so all signs are +-1) and p odd, this is the arithmetic
   core of the classification for every N. *)
Theorem single_relation_blocks_both : forall (signs : list Z) (T : Z),
  (forall x, In x signs -> x = 1 \/ x = -1) -> Nat.Odd (length signs) ->
  (zsum signs) mod 4 <> 0 /\ (zsum signs + 2*T) mod 4 <> 0.
Proof.
  intros signs T Hpm Hodd. apply odd_blocks_both.
  rewrite zsum_parity by (intros x Hx; destruct (Hpm x Hx) as [E|E]; rewrite E; reflexivity).
  now apply Nat.odd_spec.
Qed.

(* ----- Order theorem (Theorem order): the prefix-independence arithmetic core. -----
   A relation among the first K frequencies yields an anti-palindromic Q of degree 2K with
   Phi_2N | Q, so phi(2N) <= 2K; and 2K = phi(2N) is impossible because an anti-palindromic
   polynomial equal to a (palindromic) multiple of Phi_2N must vanish.  We formalize both the
   "anti-palindromic + palindromic => 0" step and the resulting index bound K >= phi/2 + 1
   (phi := phi(2N) is even).  The cyclotomic input (Phi_2N | Q, hence phi <= 2K) is supplied
   externally, exactly as for the classification. *)

(* A coefficient that equals both +b (palindromic) and -b (anti-palindromic) is zero. *)
Lemma pal_and_antipal_zero : forall a b : Z, a = b -> a = - b -> a = 0.
Proof. intros. lia. Qed.

(* Hence no relation lives in the first phi/2 indices: from phi <= 2K and 2K <> phi (phi even),
   the largest index K satisfies phi/2 + 1 <= K. *)
Lemma prefix_index_bound : forall K phi : Z,
  Z.even phi = true -> phi <= 2*K -> 2*K <> phi -> phi/2 + 1 <= K.
Proof.
  intros K phi Hev Hle Hne. apply Z.even_spec in Hev. destruct Hev as [m Hm].
  assert (phi/2 = m) by (rewrite Hm; rewrite Z.mul_comm; apply Z.div_mul; lia).
  lia.
Qed.
