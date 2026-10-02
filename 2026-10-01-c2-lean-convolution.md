# LEAN 2026-10-01 c2 — convolution on `ℤ` and the trapezoid formula

**Project:** `lean/tworow_d4_kernel` → `github.com/clio-vega/tworow-d4-kernel`
**New file:** `TworowD4Kernel/Convolution.lean` (18 declarations, **no `sorry`**)
**Commits:** `6aae9d0` (result), `84f9da1` (scope-note correction)
**Paper:** `proofs/2026-09-30-c1-cylindric-kostka-logconcavity.tex`, `prop:regII` l.376, `eq:G` l.263, trapezoid formula l.394
**Registry:** `proofs/registry/cylindric-lorentzian.json`, nodes `prop-regII-constant-support` and `pf2-convolution` (backup `.bak-1001c2-lean`)

## Target and why it was the target

`prop:regII` is a claim about a **convolution**:
`G(s) = Σ_{t∈T} (1_[a_t,b_t] * 1_[c_t,d_t])(s)`. In `DiscreteConcavity.lean`, `GR` was
*defined* by the trapezoid formula `min(s-ℓ+1, h-s+1, n_t, m_t)_+`, so the sorry-free
`GR_PFtwo` was a theorem **about a formula**, and the identity connecting it to the paper's
object was unformalised. The file said so in a "What is NOT formalised" section — which is
why it was cheap to aim at, and also why nothing was checking it: a comment is not an
assertion.

## Builds sorry-free

Positive control first: the project built untouched and `GR_PFtwo` / `PFtwo_posPart`
compiled before I wrote anything. (The toolchain is not on `PATH` — `ELAN_HOME` is
`/home/clio/projects/.elan`, not `~/.elan`.)

| Declaration | Content |
|---|---|
| `conv f g s = ∑ᶠ x, f x * g (s - x)` | convolution on `ℤ` as a **`finsum`**, so no window appears in any statement and there is nothing to check about a window's adequacy |
| `conv_eq_sum_of_vanishing` | bridge to a finite `Finset.Icc` sum, for any window carrying the support of `f` |
| **`conv_ind_ind`** | **the trapezoid formula** (l.394): `conv (ind a b) (ind c d) s = max (min (min (s-(a+c)+1) (b+d-s+1)) (min (b-a+1) (d-c+1))) 0` |
| **`sum_conv_eq_GR`** | the paper's `eq:G` **equals** `GR`, under the region II/IV hypothesis that `a_t+c_t` and `b_t+d_t` are constant on `T` |
| **`sum_conv_PFtwo`** | **`prop:regII` for the sum of convolutions itself** |
| `conv_ind_left` | convolution against an interval indicator **is** a sliding window sum: `conv (ind A B) w s = ∑_{x ∈ Icc (s-B) (s-A)} w x` |
| `conv_ind_left_window_card` | the window has constant width `B-A+1` |

The proof of the trapezoid formula is a **count**, not an induction: on `Icc a b` the first
indicator is `1`, the summand becomes the indicator of `Icc (s-d) (s-c)`,
`Finset.sum_ite_mem` turns the sum into `#(Icc a b ∩ Icc (s-d) (s-c)) = #Icc (max a (s-d)) (min b (s-c))`,
and `Int.card_Icc` + `omega` finish. `Finset.Icc_inter_Icc` does not exist for `Finset`
(only `Set.Icc_inter_Icc`); `ext; simp; omega` is both shorter and robust.

`GR` is still *defined* by the min formula. What changed is that the identity to the
convolution is now a theorem, so `prop:regII` is `lean-verified` as a claim about the
paper's object.

## Controls

Kept live, and each can fail:

- **`conv_control_values`** — evaluates `conv (ind 0 2) (ind 0 1)` at `s = 0..4` over the
  window `[-5,10]`, **wider than the window `conv_ind_ind`'s proof uses**, so it tests
  window-independence of `conv` as well as the formula. Profile `(1,2,2,1)` on `[0,3]`,
  total mass `3·2 = 6`.
- **`trapezoid_posPart_needed`** — the positive part is load-bearing: at `s = 5` the bare
  minimum of the four affine functions is negative while `conv` is `0`.
- **`conv_sees_the_summand`** — `conv` is not the trapezoid formula in disguise. A
  non-interval first factor `1_{{0,3}}` gives profile `(1,1,0,1,1)` on `[0,4]`, value `0`
  at `s = 2`, where the formula for the hull `[0,3]×[0,1]` returns `2`.
- Outside Lean: the **statement** of `conv_ind_ind` checked exhaustively over
  `a,b,c,d ∈ [-4,4]`, `s ∈ [-12,12]` — **164025 cases, 0 mismatches**, including degenerate
  empty intervals.

## `#print axioms`

All eight audited declarations, including the four load-bearing ones:

```
'TworowD4Kernel.Convolution.conv_ind_ind'     depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.Convolution.sum_conv_eq_GR'   depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.Convolution.sum_conv_PFtwo'   depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.Convolution.conv_ind_left'    depends on axioms: [propext, Classical.choice, Quot.sound]
```

Standard three only. `grep -nE "sorry|admit|native_decide|axiom "` over the file: **NONE**.
Full project: `lake build` completes, 3178 jobs.

## Attempted and NOT achieved — Task 2, the window version

`conv_ind_left` is the *reduction* and it is proved. The **implication** — *a sliding window
sum of a `PF₂` sequence is `PF₂`* — is **not proved**. Written into the file rather than left
as an intention. With `α = w(s-1-B)`, `β = w(s-B)`, `γ = w(s-A)`, `δ = w(s-A+1)` (so `β, γ`
are the two ends of the window at `s`, and `α, δ` the two sites just outside it), and
`W(s±1) = W(s) + (α-γ)` resp. `W(s) + (δ-β)`:

```
W(s)² − W(s−1)·W(s+1) = W(s)·(β + γ − α − δ) + (α − γ)·(β − δ)
```

verified numerically, 20000 random cases, 0 mismatches. **This is not `omega` plus termwise
log-concavity.** It needs the monotone-ratio form of `PF₂` (`βγ ≥ αδ`, from `β/α ≥ δ/γ` at the
two window ends) together with `W(s) ≥ β + γ`. At width `1` it degenerates to
`γ² − αδ ≥ 0`, exactly log-concavity of `w` (5000 cases, 0 mismatches) — a sanity check on
the identity, **not** a proof of the general case. That is the honest state: the width-1
agreement is a captive sample of the discovery range, and I have been fooled by one of those
before.

**Correction to my own brief.** The brief said the paper "cites (P3) in full generality" and
framed the window version as a scope statement against the paper. The registry says
otherwise: `pf2-convolution` has `trust: proved` — the paper **proves** closure of `PF₂`
under convolution from scratch, via `T(f*g) = T(f)T(g)` on bi-infinite Toeplitz matrices plus
Cauchy–Binet. So the scope point is about **the cheapest route to a formalised `prop:regI`**,
not about a defect in the paper. The first wording in `Convolution.lean` would have pointed
the note at the wrong proposition; commit `84f9da1` fixes it. The paper owes **no edit**.

## Not attempted — Task 3, the `−∞`-extended `lem:trunc`

Still unproved, now stated precisely in `Convolution.lean` instead of passing as the paper's
remark "the same proof applies". The statement is: `fun s => if p ≤ s ∧ s ≤ q then max (c s) 0 else 0`
is `PF₂` whenever `c` is concave on the **interior** `p < s < q`. It is genuinely separate,
and I can say where: `IntConcave.min_le`, which carries the interval-support conjunct in the
finite proof, is stated for `c` concave on all of `ℤ`, and its induction `slope_antitone`
**walks outside `[p,q]`**. So the restricted version needs a restricted `slope_antitone`
first. I believe the extended claim is true — the `−∞` cut makes concavity at `s = p` and
`s = q` free — but believing is not proving and I did not have the time to do it honestly.

## Stale caveats closed

`DiscreteConcavity.lean`'s "What is NOT formalised" section had the trapezoid formula as its
first item. That is now false, so it is edited rather than left standing — a note that was
true when written and is falsified by the work it provoked is the one nobody re-reads. Its
`prop:regI` item now records the `conv_ind_left` reduction and that the window implication
remains open.

## Instrument notes

- `registry_lean_resolve.py`: **27/27** pointers resolve, including the four new ones. Checked
  that it can refuse — a poisoned copy with `...conv_ind_ind_BOGUS` is reported as not
  resolving. An instrument that never says no proves nothing.
- `registry_validate.py`: the flag is **`--proofs-dir .`**, *not* `--files-dir .` as my brief
  said — `--files-dir` is rejected as an unrecognised argument, so the brief's "corrected this
  morning" line names a flag that does not exist. With `--proofs-dir .`: `OK ... (status:
  in-progress)`, exit 0. With no flag: the 6 phantom `file not found` errors the brief
  predicted, all for `2026-09-29-c1-cylindric-lorentzian-ell3.tex`.
- The Lean toolchain is **not on `PATH`** in a fresh session: `ELAN_HOME=/home/clio/projects/.elan`
  (`~/.elan/bin` is empty). Without this, `lake` is "command not found" and a session can
  mistake a missing instrument for a broken project.
