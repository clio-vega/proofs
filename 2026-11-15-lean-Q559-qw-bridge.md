# LEAN 2026-11-15: Q559, the ℚ(w) ↔ pointwise bridge (top node promoted)

**Project:** `lean/tworow_d4_kernel` (GitHub `clio-vega/tworow-d4-kernel`), module `TworowD4Kernel/QwBridge.lean`.
**Commit:** `b779be1296a5c3e43d228d3720f4aaf3ffea0d47`. Upstream was verified with `gh api .../commits/b779be1`, which returns the SHA and the file `TworowD4Kernel/QwBridge.lean`. Negative control: a bogus SHA gives HTTP 422.
**Paper:** `proofs/2026-10-23-Q559-vanishing-converse.tex`, thm:main L360–L393 and cor:criterion L395–L409.
**Plan source:** `proofs/2026-11-14-lean-cor-branch.md`, "Scoping note". The route below follows it exactly.

## Result: sorry-free

| declaration | content |
|---|---|
| `psiQ M r θ : FractionRing (MvPolynomial ι ℚ)` | `eq:Psidef` (`psiHead`, unchanged) at the generic point `w_i ↦ X_i` |
| `denom M`, `aeval_denom_ne_zero_iff` | `D = ∏_{i∈M} X_i · ∏_{s≠t∈M}(X_s − X_t)`; `D(w) ≠ 0 ⟺ Admissible w M` |
| `Spec D w x y` + closure lemmas | `x = a/Dᵏ` in ℚ(w) and `y = a(w)/D(w)ᵏ` in `K`; closed under `+ − · ∑ ∏ ^`, `/`, and under inverting divisors of `D` |
| `spec_psiHead` | **eval commutes with Ψ at admissible points** (plan step 1) |
| `spec_rep` | the specialised value is independent of the chosen representation `a/Dᵏ` |
| `psiQ_eq_zero_iff K` | **the bridge**: for any char-0 field `K`, `psiQ = 0 ⟺ ∀ admissible w : ι → K, psiHead w = 0` |
| `thm_main_Qw` | **thm:main as written**: `r ≤ |M|`, ρ-shifts distinct ⇒ `psiQ M r θ ≠ 0` in ℚ(w) |
| `cor_criterion_Qw` | **cor:criterion as written**: `psiQ = 0 ⟺ |M| < r ∨ ∃ j<k, θ_j − θ_k = k − j` |
| `inst_thm_main_Qw` | the paper's instance, θ=(0,0) on M={0,1}: `psiQ ≠ 0` |

The pointwise ⇒ generic half uses the scoped route. Take `a, k` from one admissible point; `K` is infinite because it has characteristic 0. Then `D·a` vanishes at every point: at admissible points via `spec_rep` and Ψ(w)=0, at inadmissible points because D(w)=0. `MvPolynomial.funext` on `map (algebraMap ℚ K) (D·a)`, then `map_injective`, gives D·a = 0. Since D ≠ 0, a = 0.

The generic ⇒ pointwise half: `a/Dᵏ = 0` in the fraction field forces `a = 0`.

### `#print axioms` (same run as a sorryAx control)
```
'control_sorry' depends on axioms: [sorryAx]
'TworowD4Kernel.QwBridge.psiQ_eq_zero_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.QwBridge.thm_main_Qw' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.QwBridge.cor_criterion_Qw' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.QwBridge.spec_psiHead' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.QwBridge.lesion_anyField_refutable' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.QwBridge.lesion_nonzero_refutable' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.QwBridge.inst_thm_main_Qw' depends on axioms: [propext, Classical.choice, Quot.sound]
```
`lake build TworowD4Kernel.QwBridge`: exit 0, with 0 warnings from QwBridge. This was read from `lake`'s own `$?`, with no pipe. Full `lake build`: exit 0, 3492 jobs.

## Lesions (plan step 2)

- **Drop char 0 / `Infinite K`: Lean-REFUTED** (`lesion_anyField_refutable`). At `K = 𝔽₂`, `M = {0,1}`, `r = 0` there is no admissible point, so the hypothesis is vacuous, yet `psiQ = 1`.
- **Drop `w_i ≠ 0`: Lean-REFUTED** (`lesion_nonzero_refutable`) for the generic ⇒ pointwise half. Take `M={0,1}`, `θ=(1,0)` (a ρ-shift tie, so `psiQ = 0` by `cor_criterion_Qw`) and `w=(0,2)`. Lean's `0⁻¹ = 0` gives the value 3/2. A junk value at a pole is not a specialisation.
- **Drop distinctness only: NOT refuted.** I searched exactly (Fractions, with Lean's junk conventions `x/0 = 0`, `0⁻¹ = 0`) over `n ≤ 4`, `r ≤ n`, `θ ∈ [−1,2]^r`, and tied `w ∈ {2,3,1/2,−1,5}^n`, restricted to `psiQ = 0` cases. It found **0** counterexamples (`scratch/2026-11-15-lesion-dist-search.py`). My first run reported 4 "refutations"; they were float contamination (`1/x` on Python ints, value `−1.1e−16`). They disappeared under exact arithmetic. Likely mechanism: every tied pair `s,t` either both lies in the head image (pair factor `/(1−1) = /0`), or has one inside and one outside (cross factor `/0`), or both lies outside. In the last case the surviving terms factor through Ψ on a smaller set. So distinctness may be **removable under the junk conventions**. Proving that would be new mathematics; I did not attempt it. Recorded only.

  This does **not** weaken the bridge as stated. The bridge assumes the paper's full admissibility, and the pointwise ⇒ generic half is stronger with more test points, not weaker.

## Registry
`proofs/registry/rick-q488-converse.json`:
- new child `clio-Q559-qw-pointwise-bridge`, trust `lean-verified`, with `lean` fields listed.
- **top node `clio-Q553-vanishing-converse-Q559` promoted `peer-reviewed → lean-verified`**, `lean = QwBridge.thm_main_Qw, QwBridge.cor_criterion_Qw`. Both of the node's 11-14 holds (cor:branch and the bridge) are now discharged.

`registry_validate.py` reports `OK`. Planted control: setting the bridge's trust to `lean-verfied-PLANTED` fires twice, once as invalid trust and once as the boundary rule on the promoted parent. The JSON was rewritten in the file's own format (`indent=2`), giving a 16+/4− diff. My first write used `indent=1` and produced a whole-file 1187/1175 rewrite. I caught that from `git diff --stat` before committing.

## Remaining / not done
- Distinctness-removability under the junk conventions (above) is an observation, not a question I opened.
- `TworowD4KernelTests.lean` has no QwBridge check yet.
