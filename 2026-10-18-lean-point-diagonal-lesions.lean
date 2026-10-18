import TworowD4Kernel.PointDiagonalNoGo

/-!
# Lesion suite for `TworowD4Kernel.PointDiagonalNoGo`

Each lesion DELETES one hypothesis and then proves the NEGATION of what is left, with an
explicit counterexample.  Refuting the lesioned statement grades the STATEMENT; merely
watching the original proof script fail grades the SCRIPT (5 of 6 lesions did that on
10-17 with `Unknown identifier`).  A lesion that cannot be refuted is a candidate free
strengthening and is marked EXPECT-SURVIVES.
-/

open TworowD4Kernel.PointDiagonal

namespace Lesion

/-! ## L1 — `no_single_cell_move` without `h : ∀ m ∈ S, Psi χ {m} = 1` -/
example : ¬∀ (χ : ℕ → ℚ) (S : Finset ℕ), Psi χ S = 1 := by
  intro h
  have := h (fun _ => 2) {0}
  norm_num [Psi] at this

/-! ## L2 — `psi_prod_insert` without `hm : m ∉ S` -/
example : ¬∀ (χ : ℕ → ℚ) (m : ℕ) (S : Finset ℕ),
    Psi χ (insert m S) = Psi χ {m} * Psi χ S := by
  intro h
  have := h (fun _ => 2) 0 {0}
  norm_num [Psi] at this

/-! ## L3 — `moved_cell_card_le_one` without `hother` -/
theorem prod_two : Psi (fun _ => (2 : ℚ)) ({0, 1} : Finset ℕ) = 4 := by
  rw [Psi, Finset.prod_pair (by norm_num : (0 : ℕ) ≠ 1)]
  norm_num

example : ¬∀ (χ : ℕ → ℚ) (S : Finset ℕ), Psi χ S ≠ 1 → S.card ≤ 1 := by
  intro h
  have hle := h (fun _ => 2) {0, 1} (by rw [prod_two]; norm_num)
  rw [Finset.card_pair (by norm_num : (0 : ℕ) ≠ 1)] at hle
  omega

/-! ## L4 — `moved_cell_card_le_one` without `hS : Psi χ S ≠ 1` -/
example : ¬∀ (χ : ℕ → ℚ) (S : Finset ℕ), (∀ T, T ≠ S → Psi χ T = 1) → S.card ≤ 1 := by
  intro h
  have hle := h (fun _ => 1) {0, 1} (by intro T _; simp [Psi])
  rw [Finset.card_pair (by norm_num : (0 : ℕ) ≠ 1)] at hle
  omega

/-! ## L5 — `transposition_invariance_on_diagonal` without `hf : f i = f j` -/
example : ¬∀ (f : ℕ → ℚ) (i j : ℕ) (S : Finset ℕ),
    ∏ m ∈ S.map (Equiv.swap i j).toEmbedding, f m = ∏ m ∈ S, f m := by
  intro h
  have := h (fun n => (n : ℚ)) 0 1 {0}
  norm_num [Equiv.swap_apply_left] at this

/-! ## L6 — `psi_determined_by_singletons` without `h` -/
example : ¬∀ (χ χ' : ℕ → ℚ) (S : Finset ℕ), Psi χ S = Psi χ' S := by
  intro h
  have := h (fun _ => 1) (fun _ => 2) {0}
  norm_num [Psi] at this

/-! ## L7 — `no_point_diagonal_moves_exactly_a_pair` without `hab : a ≠ b` -/
example : ¬∀ (χ : Unit → ℚ) (a b : Unit),
    ¬(Psi χ {a, b} ≠ 1 ∧ ∀ T, T ≠ ({a, b} : Finset Unit) → Psi χ T = 1) := by
  intro h
  refine h (fun _ => 2) () () ⟨?_, ?_⟩
  · rw [Finset.pair_eq_singleton]
    norm_num [Psi]
  · intro T hT
    rw [Finset.pair_eq_singleton] at hT
    have hT' : T = ∅ := by
      by_contra h0
      obtain ⟨x, hx⟩ := Finset.nonempty_iff_ne_empty.mpr h0
      refine hT (Finset.eq_singleton_iff_unique_mem.mpr ⟨?_, fun y _ => Subsingleton.elim y ()⟩)
      rwa [Subsingleton.elim ((() : Unit)) x]
    simp [hT', Psi]

/-! ## L8 — `konckiFamily_ne_std` without `hab : a ≠ b`.  EXPECT-SURVIVES. -/
theorem konckiFamily_ne_std_no_hab {L : Type*} [DecidableEq L] {K : Type*} [Field K]
    {a b : L} {c : K} (hc : c ≠ 0) : konckiFamily a b c a ≠ (std a : L → K) := by
  intro h
  rw [konckiFamily_self] at h
  rcases eq_or_ne b a with rfl | hba
  · have hb := congrFun h b
    simp only [Pi.add_apply, Pi.smul_apply, std, Pi.single_eq_same, smul_eq_mul,
      mul_one] at hb
    exact hc (by simpa using hb)
  · have hb := congrFun h b
    simp only [Pi.add_apply, Pi.smul_apply, std, Pi.single_eq_of_ne hba, Pi.single_eq_same,
      smul_eq_mul, mul_one, zero_add] at hb
    exact hc hb

/-! ## L9 — `konckiFamily_not_diagonal` without `hab : a ≠ b`.  Refuted by an explicit
witness: at `a = b` the transvection IS diagonal, with scaling `1 + c`. -/
example : ¬∀ (a b : Finset ℕ) (c : ℚ), c ≠ 0 →
    ∀ d : Finset ℕ → ℚ, ¬∀ x : Finset ℕ, konckiFamily a b c x = d x • (std x : Finset ℕ → ℚ) := by
  intro h
  refine h ∅ ∅ 1 one_ne_zero (fun x => if x = ∅ then 2 else 1) ?_
  intro x
  rcases eq_or_ne x ∅ with rfl | hx
  · rw [konckiFamily_self]
    funext y
    by_cases hy : y = (∅ : Finset ℕ) <;>
      simp [std, Pi.single_apply, hy] <;> norm_num
  · rw [konckiFamily_of_ne _ _ _ hx]
    funext y
    simp [std, Pi.single_apply, hx]

/-! ## L10 — `konckiFamily_not_diagonal` without `hc : c ≠ 0` -/
example : ¬∀ (a b : Finset ℕ) (c : ℚ), a ≠ b →
    ∀ d : Finset ℕ → ℚ, ¬∀ x : Finset ℕ, konckiFamily a b c x = d x • (std x : Finset ℕ → ℚ) := by
  intro h
  refine h ∅ {0} 0 (by decide) (fun _ => 1) ?_
  intro x
  rcases eq_or_ne x ∅ with rfl | hx
  · rw [konckiFamily_self]
    funext y
    simp [std, Pi.single_apply]
  · rw [konckiFamily_of_ne _ _ _ hx]
    funext y
    simp [std, Pi.single_apply]

/-! ## L11 — `cKoncki_ne_zero`, each of the three hypotheses in turn -/
example : cKoncki (0 : ℚ) 2 2 = 0 := by norm_num [cKoncki]
example : cKoncki (1 : ℚ) 2 2 = 0 := by norm_num [cKoncki]
example : cKoncki (2 : ℚ) 1 1 = 0 := by norm_num [cKoncki]

end Lesion
