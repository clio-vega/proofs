# LEAN 2026-10-02 — `tail-of-pf2-logconcave`

**Target.** Registry `proofs/registry/cylindric-lorentzian.json`, node
`tail-of-pf2-logconcave` (child of `A-m2-complete`). Paper proof:
`proofs/2026-10-01-c2-inter-region-inequality.tex` §`sec:tail`, `lem:tail`.

**Project.** `lean/tworow_d4_kernel`, new module `TworowD4Kernel/TailSum.lean`,
imported from the root aggregator. `lake build` exit 0, **3179 jobs** (was 3178).

**Result: sorry-free. Sorry count 0.** Node promoted `proved` → `lean-verified`.

---

## Priority 0 — decided before writing any Lean

I printed the node before writing the sentence. `pf2-convolution` is **`trust: proved`**:
the paper proves PF₂-closure-under-convolution from scratch (Toeplitz `T(f*g)=T(f)T(g)` +
Cauchy–Binet). An interval indicator is PF₂, and `conv_ind_left` is already in
`Convolution.lean`. So *window sum of PF₂ is PF₂* is an **immediate corollary of a proved
node**, and the 2026-10-01 annotation inside that node's own `approach` field reading
**"NOT YET PROVED"** is false as a statement about mathematics. It is true only as a
statement about Lean.

**Decision:** the honest wording is *not yet formalised*. I corrected the annotation in
place (`NOT YET PROVED` → `NOT YET FORMALISED`, with the reason attached) and did **not**
act on it further — the route question (Toeplitz+Cauchy–Binet vs the direct identity)
stays open, and the direct route's two extra inequalities may be true, provable and
entirely surplus.

This is the third time in four days a brief's verdict *about a source* was wrong while its
mathematics was sound. The brief predicted a gap in `sec:tail` too: **there is none.** Both
degenerate branches (`Δ ≤ 0`, and `β(y)=0` with `β(y-1)>0`) are written out in the paper
and both are correct. A brief that predicts a gap is still not evidence of one.

---

## What was formalised

```
noncomputable def tail (β : ℤ → ℤ) (N y : ℤ) : ℤ := ∑ r ∈ Finset.Icc y N, β r

theorem tail_logConcave (hβ : PFtwo β) (hN : ∀ r, N < r → β r = 0) (y : ℤ) :
    tail β N (y-1) * tail β N (y+1) ≤ tail β N y * tail β N y

theorem tail_PFtwo (hβ : PFtwo β) (hN : ∀ r, N < r → β r = 0) : PFtwo (tail β N)
```

Plus `tail_pos_downSet` ( `{γ>0}` is a down-set — the paper's third clause, stated
separately because `PFtwo.suppInterval` is strictly weaker), `tail_antitone`,
`tail_nonneg`, `tail_succ`, `tail_eq_zero_of_lt`, `zero_of_right`.

`PFtwo` is the project's existing structure in `DiscreteConcavity.lean`, so this composes
directly with everything else; the conclusion is **stronger than the node asked** (`PFtwo γ`,
not just log-concavity of `γ`).

### Where the Lean proof departs from the paper, and why

The paper divides: `q = β(y)/β(y-1) ∈ (0,1)`, `β(y+k) ≤ β(y)qᵏ`, sum the geometric series,
`γ(y) ≤ β(y)/(1-q)`. Over `ℤ` there is no division, so the geometric series is replaced by
the equivalent **multiplicative induction** on the upper limit:

```
tail_bound_aux :  (β(y-1) - β y) * (∑_{r=y}^{M} β r)  +  β(y-1) * β(M+1)  ≤  β(y-1) * β y
                  for every M ≥ y-1
```

The second summand is exactly the geometric remainder the paper discards by summing to
infinity; in the geometric case the inequality is an **equality**, so nothing is given away.
At `M = N` it vanishes (`β(N+1)=0`) and the paper's `Δ·γ(y) ≤ β(y-1)β(y)` is recovered
verbatim. The induction step consumes one input,

```
ratio_antitone :  β(y-1) * β(m+1) ≤ β y * β m      for all m ≥ y   (given 0 < β y)
```

which is the paper's "the ratios `β(r+1)/β(r)` form a nonincreasing sequence", cleared of
denominators. This is strictly better than transcribing the paper: no `q`, no `(0,1)`
membership, no series convergence, no real numbers.

### Where each hypothesis is spent (one place each)

| hypothesis | spent at |
|---|---|
| `β ≥ 0` | `tail_nonneg`; the `Δ ≤ 0` branch |
| log-concavity | base case **and** step of `ratio_antitone` — nowhere else |
| interval support | `zero_of_right` **only**, and only to push a zero rightwards |
| `β` vanishes above `N` | `tail_succ`, so the recurrence holds at *every* `y ∈ ℤ` |

Concavity is **not** assumed. The paper's first version assumed truncated-concavity and the
refusal control for that clause did not fire; that surplus hypothesis is gone and stayed gone.

---

## Controls — all three fire, all three sorry-free

I first asked whether interval support was derivable from global log-concavity (it would
have been a stronger theorem). **It is not**, and a 30-second small example settled it:
`(1,0,0,1)` is globally log-concave with a gap. The project already knew this —
`DiscreteConcavity.gapTwo_not_PFtwo`. So the registry is right that interval support is
load-bearing, and the hypothesis stays.

1. **`ctrlA_tail_not_logConcave`** — isolates **log-concavity**. `ctrlA = (2,1,2)` on
   `{0,1,2}`: nonnegative, support *is* an interval, log-concavity fails at `1`
   (`2·2 > 1·1`). Tail `(5,3,2,0)` fails at `y=1`: `5·2 = 10 > 9 = 3²`. ✅ fires
2. **`ctrlB_tail_not_logConcave`** — isolates **interval support**. `ctrlB = (1,0,0,1)`:
   nonnegative, **globally log-concave** (proved in Lean, `ctrlB_logConcave`), support not an
   interval. Tail `(2,1,1,1,0)` fails at `y=1`: `2·1 = 2 > 1 = 1²`. ✅ fires
3. **`converse_false`** — `∃ β γ, PFtwo γ ∧ (∀ y, γ y = β y + γ (y+1)) ∧ ¬ PFtwo β`, with
   `γ = (…,4,4,4,2,1,0,…)` and `β = (2,1,1)` (`β(0)β(2) = 2 > 1 = β(1)²`). So nobody may read
   `tail_PFtwo` backwards. ✅ fires

   *Stated through the defining recurrence rather than through `tail`*, deliberately: `tail`
   is cut at a finite `N` while `γ` is constant on all of `ℤ_{≤0}`, so `γ = tail β 2` would
   need a separate downward induction that adds nothing. The recurrence plus `γ(y)=0` for
   `y>2` pins `γ` down completely. (Checked numerically that the identification does hold.)

Also proved as refusals: `ctrlA_not_logConcave`, `ctrlB_not_suppInterval`.

---

## Independent cross-check (a different *kind* of instrument)

A Python sweep, written separately from the Lean:

- **200000** random PF₂ instances, every `y` in range: **0 failures**.
- `tail_bound_aux` including the discarded remainder term: **20000+ instances, 0 violations**.
- Control tails reproduce the values Lean's `omega` derived: `ctrlA → [5,3,2,0]`,
  `ctrlB → [2,1,1,1,0]`. `ctrlGamma = tail(ctrlC)` confirmed on `[-6,6]`.
- **Instrument validated before its zeros were believed**: the sweep was first shown to
  *detect* the known `ctrlA` failure. A sweep that cannot find a planted failure cannot
  certify an absence.

---

## `#print axioms`

Every new declaration, exactly the standard three — which also rules out `sorryAx`:

```
'TworowD4Kernel.TailSum.tail_logConcave'              [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TailSum.tail_PFtwo'                   [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TailSum.tail_pos_downSet'             [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TailSum.ratio_antitone'               [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TailSum.tail_bound_aux'               [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TailSum.ctrlA_tail_not_logConcave'    [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TailSum.ctrlB_tail_not_logConcave'    [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TailSum.converse_false'               [propext, Classical.choice, Quot.sound]
```

---

## Validators — refusal tested before the OK was believed

| instrument | planted violation | clean |
|---|---|---|
| `registry_lean_resolve.py` | poisoned declaration name → **exit 1**, "does not resolve" | **exit 0**, 35 pointers, 794 declarations indexed |
| `trustcheck.py … validate … --files-dir .` | node `file` → nonexistent path → **exit 1**, "not found under ." | **exit 0** |

No pipes were used when reading `$?`.

---

## "by Lemma X" sweep

Every citation inside the new proofs checked for *X's object being defined on the case
invoked*:

- `zero_of_right` is invoked three times, each with the positivity witness strictly to the
  left of the zero (`y-1 < y`, `y ≤ m+1`) — its `hym : y ≤ m` is satisfied in each.
- `ratio_antitone` is invoked in `tail_bound_aux` at `m = M+1` under `M ≥ y-1`, so `y ≤ m`
  holds; and it requires `0 < β y`, which is the branch hypothesis.
- `tail_bound_aux` is invoked at `M = N` under `y ≤ N+1`, i.e. `y-1 ≤ N`. The complementary
  branch `y > N+1` is handled separately (empty `Icc`, tail `= 0`).
- `tail_succ` has no side condition — this is why `hN` is in the hypotheses rather than a
  `y ≤ N` guard.

No defect found. (On 10-01 this sweep found 6 sites, 5 sound and 1 a real defect.)

---

## What remains open, and why it is separate

- **`prop:regI` / window-sum-of-PF₂** — *proved, not formalised* (see Priority 0). Route
  undecided: Toeplitz+Cauchy–Binet (general, expensive) vs the direct identity (two extra
  inequalities that may be surplus). Not a gap in the mathematics.
- **`lem:bdry`, `lem:notrunc`, `lem:tele`, `thm:main` of Lemma T** — untouched. `lemma-T` is
  known to be strictly stronger than what (A) at m=2 needed (`lemma-T-not-needed-for-m2`),
  so it is off the critical path.
- **`-∞`-extended `lem:trunc`** — not reopened. `IntConcave.min_le`'s induction
  `slope_antitone` walks outside `[p,q]`.
- **`lemma-T-interior-concavity`** — deliberately *not* promoted; nothing this session
  touched it.

A note on the registry sentence I corrected: a false "NOT YET PROVED" sitting inside a node
whose own `trust` is `proved` is self-sealing — it suppresses the re-read that would expose
it. That is why it was corrected in the same action as the grade rather than noted for later.
