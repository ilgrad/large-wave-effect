(* Staggered-current lemma (paper: Lemma lem:stagger-leading) -- algebraic kernel, machine-checked.

   Context (linear DNLS on the line, single-site seed):
     a_n(t) = exp(2 i t) * (-i)^n * J_n(2t),   J_n real.
   The phase ratio (-i)^(n+2)/(-i)^n = (-i)^2 = -1 is REAL (verified in Maxima,
   numerics/verify_staggered_current.mac), and J_n is real, so
     a_{n+2} = ( -J_{n+2}(2t)/J_n(2t) ) * a_n
   is a REAL scalar multiple of a_n.

   This file machine-checks the remaining algebraic step: if w is a real scalar multiple of z, then the
   bilinear "current"  Im( conj(z) * w )  vanishes.  With z = a_n, w = a_{n+2}, the distance-2 current
   Im(conj(a_n) a_{n+2}) = 0, hence the RHS of the staggered-momentum identity (paper, eq:stagger)
     d/dt Re sum_n conj(a_n) a_{n+2}  =  -gamma sum_n (|a_n|^2 - |a_{n+2}|^2) Im(conj(a_n) a_{n+2})
   vanishes on the linear flow -- the leading (order-gamma) staggered-momentum drift is identically zero. *)

From Stdlib Require Import Reals.
Open Scope R_scope.

(* complex numbers as pairs of reals (Coquelicot is not installed; stay Stdlib-only, as in the rest of the layer) *)
Record C : Set := mkC { re : R ; im : R }.

Definition cmul (z w : C) : C :=
  mkC (re z * re w - im z * im w) (re z * im w + im z * re w).
Definition cconj (z : C) : C := mkC (re z) (- im z).
Definition cscale (r : R) (z : C) : C := mkC (r * re z) (r * im z).

(* the bilinear current  Im( conj(z) * w ) *)
Definition current (z w : C) : R := im (cmul (cconj z) w).

(* a real scalar multiple carries no current *)
Lemma current_real_multiple (z : C) (r : R) :
  current z (cscale r z) = 0.
Proof. unfold current, cmul, cconj, cscale; simpl; ring. Qed.

(* hence: if w is any real multiple of z, the distance-2 current vanishes *)
Corollary staggered_current_vanishes (z w : C) (r : R) :
  w = cscale r z -> current z w = 0.
Proof. intros ->; apply current_real_multiple. Qed.
