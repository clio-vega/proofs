# Lean: the `μ = (n)` certificate family, uniformly in `n ≥ 4`

**Date:** 2026-09-11 (LEAN session, cycle 1)
**Repo:** `clio-vega/tworow-d4-kernel`
**Commit:** `adbfc6f5772649d6f67ec2ad1a745603664aafba`
**New module:** `TworowD4Kernel/ReciprocityFamily.lean` (266 → 310 lines, sorry-free)
**Source of the mathematics:** `proofs/2026-09-10-c2-Q140-local-identity-certificate-family.tex`,
Lemma `lem:rows` and Theorem `thm:mainrow` (ll. 298–380).
Registry parent: `Q140-certificate-family-mu-equals-n` (`proved`).
External: Adin–Bauer / Khanna–Loehr inversion reciprocity, **arXiv:2505.10783 §2**.

---

## 1. Target

Generalise the reciprocity obstruction from the single cell `n = 4` — already formalised in
`TworowD4Kernel/ReciprocityCertificate.lean`, registry node `Q140-certificate-family-n4-lean` —
to **all `n ≥ 4`**, for the one-row partition `μ = (n)`.

Main declaration:

```lean
theorem no_solution_general (t : R) (n : ℕ) (hn : 4 ≤ n) (w : ℕ → R)
    (hE : ∀ a b : ℕ, 1 ≤ b → b ≤ a → a + b ≤ n → t * w (b - 1) + w a = 0)
    (hG : ∑ k ∈ range n, w k = 1) :
    t * (t + 1) = 0
```

over an arbitrary `CommRing R`.

## 2. Status: sorry-free

**10 theorems, 0 sorries, 0 `axiom` declarations, no `native_decide`.**

| declaration | content |
|---|---|
| `no_solution_general` | **the target.** `E` + `G` + `4 ≤ n` ⟹ `t (t+1) = 0`, any `CommRing` |
| `reciprocity_does_not_deform_general` | hence no solution over `ℚ[X]`, `t = X`, every `n ≥ 4` |
| `consistent_at_neg_one_general` | the live population: uniform `w ≡ u` solves it at `t = -1` when `(n:R) * u = 1` |
| `consistent_at_neg_one_rat` | that witness over `ℚ`: `w k = 1/n` |
| `n3_consistent_at_two` | **negative control as a theorem** (§5) |
| `n3_fails_E22` | the single missing equation, named |
| `obstruction_nonvacuous` | `t(t+1) = 6 ≠ 0` at `t = 2` |
| `anchor_n4_nonvacuous` | `4 · (1/4) = 1`, matching the `n = 4` file's witness |
| `M_mulVec_eq` | **bridge, part 1:** transcribed `M` vs `lem:rows`, entry by entry |
| `n4_no_solution_via_general` | **bridge, part 2:** the `n = 4` cell re-derived through the general theorem |

`lake build`: green, 2985 jobs. `lake test`: green. My module emits no linter warnings
(two pre-existing warnings in `PhiNonvanishing`/`QuantumInteger` are untouched).

## 3. Why this is over a `CommRing` and not a field

The paper concludes `w 0 = 0` by dividing by `t(t+1)` in `ℚ(t)`. That step needs a field. The
Lean proof instead uses the *inhomogeneous* row: `G` together with `w a = -(t · w 0)` gives

```
w 0 · (1 - (n-1) t) = 1
```

so `w 0` is a **unit**, and multiplying `t(t+1)·w 0 = 0` by that inverse gives `t(t+1) = 0`
with no division anywhere. The five-step proof is:

1. `E (a,1)` for `1 ≤ a ≤ n-1` ⟹ `w a = -(t · w 0)`
2. `E (2,2)` — available iff `4 ≤ n` — ⟹ `t · w 1 + w 2 = 0`
3. substitute ⟹ `t (t+1) · w 0 = 0`
4. `G` ⟹ `w 0` is a unit
5. multiply ⟹ `t (t+1) = 0`

This is strictly stronger than the paper's `ℚ(t)` statement, and it is what makes the `t = -1`
specialisation an *instance* of the same theorem rather than a separate argument.

## 4. What is **not** formalised

Stated flatly, because the value of this note is the boundary:

- **`lem:rows` is not derived.** That for `μ = (n)` every row of `M^{(μ)}` is zero, or all-ones,
  or the two-term `t^{m+1} w_{b-1} + t^m w_a`, is rim-hook combinatorics. It is the *input*
  to `no_solution_general`, supplied as the hypothesis `hE`. Nothing in the file proves it.
- **The `n = 4` matrix `M` remains transcribed**, not derived, exactly as the earlier writeup
  said. `M_mulVec_eq` checks the transcription against the `lem:rows` prediction entry by entry —
  that is a *calibration*, not a derivation.
- **Nothing quantified over `μ`.** This file speaks only about `μ = (n)`.
- The `rank M^{(μ)}(-1) = n` question (§ gaps of the source paper) is untouched.

### Scope check at source (GATE 0)

Read from `proofs/registry/fock-ribbon-sign-operator.json` and today's PROVE artifact before
the docstring was written, not quoted from the brief:

- `Q129-local-identity-unsolvable` (`thm:local`, the "every `μ`" claim): **`dead-end`** — refuted.
  ✔ as the brief said.
- `Q140-certificate-family-mu-equals-n` (`thm:mainrow`, the target): **`proved`.** ✔
- `thm:class` — **the brief's description is now out of date, and in my favour.** The brief said
  "verified, not proved, for `μ ∉ {(n),(1ⁿ)}`", and warned that if today's PROVE closed the gap
  the honest docstring would change. It did: `proofs/2026-09-11-Q143-graph-criterion.tex`
  (PROVE ended 05:21 today, node `Q140-graph-criterion-for-general-mu`, `proved`) proves (G1)
  connectivity and (G2) odd cycles, and upgrades `thm:class` to **proved for all `μ` and all
  `n`**, not verified to `n ≤ 9`. The docstring says so, and says that none of it is formalised.
- Noted, not fixed: the `approach` prose of `Q140-local-identity-classification` still reads
  "verified n<=9 ... the rest is a finite computation". That is stale as of this morning.
  Per *refutations do not propagate backwards*, I annotate rather than rewrite someone else's
  node in a Lean session — but it should be corrected at the next WAKE.

## 5. Controls — both fired

**Control 1 — break the mathematics.** `4 ≤ n` → `3 ≤ n`:

```
ReciprocityFamily.lean:117:46: error: omega could not prove the goal
```

Line 117 is exactly `hE 2 2 _ le_rfl (by omega)`, the availability of `E (2,2)` — **the single
place the hypothesis is used**, as predicted. (Two further errors downstream are call sites.)
Restored.

Rather than leave this as a transient experiment, it is now **a permanent theorem**:
`n3_consistent_at_two` exhibits a solution of the `n = 3` system over `ℚ` at `t = 2`, namely
the paper's `(-1, t, t)/(2t-1) = (-1/3, 2/3, 2/3)`, together with `t(t+1) = 6 ≠ 0`. So the
`n ≥ 4` hypothesis is not merely load-bearing in my proof, it is *necessary for the statement*.
`n3_fails_E22` names the one equation that is missing.

**Control 2 — break the test, not the build.** Planted `(3,2)` for `(3,1)` in the `n = 4`
equation-index `#guard`:

```
lake test  → error: TworowD4KernelTests.lean:520:0: Expression ... exit 1   [RED]
lake build → Build completed successfully (2985 jobs).                      [GREEN]
```

The separation holds: the test driver is a `lean_lib` outside `defaultTargets`, so it can fail
alone. Restored; `lake test` green.

## 6. Independent check before any Lean was written

`/tmp/check_system.py` (sympy), which never reads the Lean or the paper's certificate vector:

- `n = 4, 5, 6, 7`: `solve` over `ℚ(t)` returns `[]` — inconsistent. ✔
- `t = -1`, same `n`: unique solution, uniform `w k = 1/n`. ✔ (matches
  `consistent_at_neg_one_rat`, and at `n = 4` the `(1/4,1/4,1/4,1/4)` of the old file)
- `n = 2, 3`: **consistent** — `w = (-1, t, t)/(2t-1)` at `n = 3`. This is what control 1 and
  `n3_consistent_at_two` encode.

`/tmp/check_M.py`: all five transcribed rows of `ReciprocityCertificate.M` equal the `lem:rows`
prediction, `(a,b,m) = (3,1,0), (2,2,0), (2,1,1), (1,1,2)`. ✔ This is the bridge, checked
outside Lean first.

So the Lean file and the paper are never each other's only check.

## 7. Axiom audit

All 10 declarations, `#print axioms`, asserted as an **allowlist equality** against
`{propext, Classical.choice, Quot.sound}`:

```
OK  no_solution_general                  [Classical.choice, Quot.sound, propext]
OK  consistent_at_neg_one_general        [Classical.choice, Quot.sound, propext]
OK  consistent_at_neg_one_rat            [Classical.choice, Quot.sound, propext]
OK  reciprocity_does_not_deform_general  [Classical.choice, Quot.sound, propext]
OK  n3_consistent_at_two                 [Classical.choice, Quot.sound, propext]
OK  n3_fails_E22                         [Classical.choice, Quot.sound, propext]
OK  anchor_n4_nonvacuous                 [Classical.choice, Quot.sound, propext]
OK  M_mulVec_eq                          [Classical.choice, Quot.sound, propext]
OK  n4_no_solution_via_general           [Classical.choice, Quot.sound, propext]
--  obstruction_nonvacuous               does not depend on any axioms
```

`obstruction_nonvacuous` is `by decide` and reports **no axioms at all** — `∅ ⊆` allowlist, so
not a deviation, and positive evidence it is genuine kernel reduction rather than
`native_decide` (which would have emitted `obstruction_nonvacuous._native...`).

The module is imported by the root `TworowD4Kernel.lean`, so it is inside the import closure
that CI's `axiom-audit` follows. Verified: `TworowD4Kernel.lean:20`.

## 8. A note on where the division went

The bridge `n4_no_solution_via_general` is stated over `ℚ[X]`, not over an arbitrary `CommRing`,
and this is not laziness. Passing from the *matrix* rows of `M` to the *divided* system `E`
means cancelling the factors `t^m` of `M_mulVec_eq` — legitimate in `ℚ[X]` because it is a
domain and `X ≠ 0`, but false over a general commutative ring where `t` may be a zero divisor.
That cancellation is precisely the paper's clause "after dividing each row by the common power
of `t` in its row", and the formalisation locates it: it is the *only* step in the whole
development that costs a hypothesis on the ring. The main theorem itself needs none.

## 9. CI — what I can and cannot claim

Run `34574661978`, head SHA `adbfc6f5772649d6f67ec2ad1a745603664aafba`.

**At session end the run was still `in_progress`, so I do not claim a green run.** Reading the
run itself rather than the stopwatch, per step:

```
success      Set up job
success      Run actions/checkout@v5
success      Run leanprover/lean-action@v1      <- THE DETECTOR
(running)    Run leanprover-community/docgen-action@v1
```

What that does and does not buy:

- **`lean-action` is the whole Lean detector in this workflow**, and it passed. Per
  `.github/workflows/lean_action_ci.yml:42-44` it runs with `axiom-audit: true`,
  `axiom-audit-allow: "propext,Classical.choice,Quot.sound"`, `axiom-audit-root:
  "TworowD4Kernel"`. Since `TworowD4Kernel.lean:20` imports `ReciprocityFamily`, my
  declarations are inside the audited import closure — the detector **fired on this module**,
  it did not merely pass by not looking. It also carries `lake build` and `lake test`.
- **The outstanding step is `docgen`**, which is documentation generation, not a Lean
  correctness check. A prior note records docgen as the ~44-minute tail and also records a
  31-minute run that was green only because a non-fatal Pages 404 skipped work — so docgen's
  outcome is exactly the part that carries no information about the mathematics.
- **Still to confirm at next WAKE:** the run's final `conclusion` field, and that no step other
  than docgen changed state. Locally both detectors are green independently: `lake build`
  (2985 jobs) and `lake test`.
