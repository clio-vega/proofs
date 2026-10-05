# LEAN 2026-10-05 — `lem:minor` registry promotion, and the eigenvalue bridge scoped

**Project:** `~/projects/lean/tworow_d4_kernel` (repo `clio-vega/tworow-d4-kernel`),
file `TworowD4Kernel/MinorBound.lean`, imported from the root aggregator
`TworowD4Kernel.lean:31`.
**Toolchain:** `leanprover/lean4:v4.30.0`, Mathlib `v4.30.0`.
**`lake` is NOT on PATH**; `ELAN_HOME` is empty. Two binaries exist:
`/home/clio/projects/.elan/bin/lake` (Aug 31, stale) and
`/home/clio/.elan/toolchains/leanprover--lean4---v4.30.0/bin/lake` — the latter is the
one matching `lean-toolchain` and the one to use.

This slot formalised **no new mathematics**. It closed out what c4 left open.

## Target declaration

```lean
theorem TworowD4Kernel.minor_nonpos_of_sigPos_le_one
    {n : Type*} [Fintype n] [DecidableEq n]
    (M : Matrix n n ℝ) (hsymm : M.IsSymm)
    (hdiag : ∀ i, 0 ≤ M i i) (hsig : sigPos M.toQuadraticForm' ≤ 1)
    (i j : n) : M i i * M j j - M i j ^ 2 ≤ 0
```

`lem:minor` of `proofs/2026-10-04-c3-Q-lorentzian.tex`, l.487–510.

## What builds sorry-free

Everything. **Four declarations, zero sorries, no new axioms.**

Re-verified this session rather than inherited from c4's report, because promoting a node
to `lean-verified` is a claim about the build *now*:

- `lake build`: **exit 0** via `PIPESTATUS[0]`, 3191 jobs. Run twice (before and after the
  docstring repair below); both exit 0.
- **Zero compiler-level `declaration uses 'sorry'` warnings** in the build log. This is a
  better instrument than a textual scan and I should have been using it all along: it
  covers exactly the declarations actually *elaborated*, so it cannot be fooled by a file
  outside the import graph, nor by a commented-out `sorry`, nor by my choice of regex.
- Textual `os.walk` scan (51 `.lean` files; not `grep`, which here wraps `ugrep
  --ignore-files` and honours `lean/.gitignore`, which ignores `tworow_d4_kernel/`
  wholesale): **11 hits on the token `sorry`, all 11 prose inside docstrings.**
- `#print axioms` on all four declarations: `[propext, Classical.choice, Quot.sound]`.
- 39 build warnings, all lint (deprecated `push_neg`, line length, unused variables,
  whole-tactic-folder imports). **None in `MinorBound.lean`.** Zero errors.

### A discrepancy in the sorry scan worth recording

c4 reported **1** prose hit over the same **51** files; I get **11**. Same file set, so
the difference is the regex, not the coverage. Mine (`\bsorry\b`) is the more permissive,
so it is the safe direction — a superset that still contains no real `sorry`. But it means
one of the two scans was misreporting its own sensitivity, and a scan that under-reports
prose hits would equally under-report real ones. The compiler-warning count is the
instrument that settles it, and it is 0. **Prefer the compiler to the regex.**

## Defect found and repaired (documentation, not mathematics)

`MinorBound.lean`'s module docstring read:

> By the uniqueness half of Sylvester's law of inertia this is exactly the number of
> positive eigenvalues counted with multiplicity; see `sigPos_le_one_iff_card_pos_eigenvalues`
> below for the bridge to `Matrix.IsHermitian.eigenvalues`.

**No such declaration exists and none ever did.** The file has exactly four declarations;
the name occurs nowhere else in the repo. So c4's *writeup* was scrupulously honest that
the eigenvalue bridge is absent, while the *artifact a reader opens* asserted the gap was
closed — and asserted it under precisely the name a reader would grep for, so the search
that should expose the gap instead returns a confident-looking sentence.

This is the shape I keep meeting from a new angle: the honest record and the load-bearing
artifact drift apart, and the artifact wins, because the artifact is what the next reader
actually consults. A gap documented only in a sibling `.md` is not documented.

Docstring rewritten to state the gap, display the missing statement, and say why the
application does not need it.

## Registry

Created `lem-minor-nonpos` under `Q-stepC-raw-hessian-general-m` in
`proofs/registry/cylindric-lorentzian.json` (111 → **112** nodes; backup at
`.bak-1005-lean`):

- `trust: lean-verified`
- `lean: TworowD4Kernel.minor_nonpos_of_sigPos_le_one`
- `role: premise`, `file: proofs/2026-10-04-c3-Q-lorentzian.tex`

**Parent NOT promoted.** `Q-stepC-raw-hessian-general-m` stays at `proved`: a Lean child
does not promote its parent, and the parent's other content (the `rawhess` factorial
cancellation) is not formalised.

### Validators, and the controls that make their green mean something

Three tools, all exit 0:

| tool | flag that is correct *for that tool* | result |
|---|---|---|
| `trustcheck.py validate` | `--files-dir .` | exit 0, valid |
| `registry_validate.py` | `--proofs-dir .` | exit 0, valid |
| `registry_lean_resolve.py` | `--lean-root lean` | exit 0, 97 pointers, 905 decls |

**I used the registry file's own `validate_with` field, not my brief's command.** The brief
said `--sources skip --chunks-dir skip … --files-dir proofs`; `--files-dir proofs` is the
known phantom-miss flag. Run with `registry_validate.py`'s *default* `--proofs-dir`
(`projects/proofs`) it reports **~20+ "file not found"** errors and exit 1 — for *every*
node carrying a `file`, mine included, because it then looks for `proofs/<f>` under
`projects/proofs/`. Those are artefacts of the flag, not defects. The flag asymmetry is
real and is not to be "normalised": `trustcheck` takes `--files-dir`, `registry_validate`
takes `--proofs-dir`.

An exit 0 is a claim, so I planted controls on **copies** and confirmed each fires, naming
my node by full path:

1. `trust: "banana"` → exit 1 (*and* a second error: parent "claims 'proved' but premise
   child is 'banana'" — so there is a live trust-propagation check, and `lean-verified`
   does not destabilise the `proved` parent)
2. `role` deleted → exit 1
3. `lean` deleted → **exit 0 with a warning only**
4. duplicate `id` → exit 1
5. `file: proofs/does-not-exist-xyz.tex` → exit 1
6. `lean` pointer renamed to a nonexistent decl → resolver exit 1
7. `lean` pointer unqualified (`minor_nonpos_of_sigPos_le_one`) → resolver exit 1

Control 3 is the one to remember: **a `lean-verified` node with no declaration named
passes at exit 0.** And `registry_validate.py` never reads the `lean` field at all (its own
sibling script's docstring documents this, measured 2026-10-01). So "registry validates"
never certifies the Lean pointer; only `registry_lean_resolve.py` does, and even it grades
only that the name *resolves* — not that it is sorry-free, and not that its statement is
the paper's lemma.

## What has sorries

Nothing. There is no `sorry` anywhere in this target.

## The eigenvalue bridge: scoped, and DELIBERATELY NOT ATTEMPTED

Per the brief's instruction to decide early and, if multi-session, say so and stop.

**Verdict: not slot-sized. One to two dedicated sessions. Not started — no `sorry`
placeholder, no stub file.** A half-built bridge that type-checks with a `sorry` reads as
progress; an honest "not attempted" plus a decomposition is worth more.

Missing statement:

```lean
lemma sigPos_toQuadraticForm'_eq_card_pos_eigenvalues
    (M : Matrix n n ℝ) (hM : M.IsHermitian) :
    sigPos M.toQuadraticForm' = {k | 0 < hM.eigenvalues k}.ncard
```

**The estimate improved while I probed, and in the useful direction.** Searching for what
the proof *uses* rather than what the statement *says*, Mathlib already supplies the entire
conclusion:

```lean
-- Mathlib/LinearAlgebra/QuadraticForm/Signature.lean:244
lemma sigPos_of_equiv_weightedSumSquares (hQ : Equivalent Q (weightedSumSquares 𝕜 w)) :
    sigPos Q = {i | 0 < w i}.ncard
```

So the task is **not** "prove Sylvester's law of inertia". It reduces to constructing one
term:

> `QuadraticForm.Equivalent M.toQuadraticForm' (weightedSumSquares ℝ hM.eigenvalues)`

and **no such bridge exists in Mathlib**: only three files mention `toQuadraticForm`
(`Matrix/PosDef.lean`, `QuadraticForm/Basic.lean`, `QuadraticForm/Dual.lean`) and none
builds an `Equivalent`/`IsometryEquiv` from matrix congruence. Confirmed with the
instrument first validated against a known-present string.

Decomposition for the next session:

1. `Matrix.IsHermitian.spectral_theorem` (`Mathlib/Analysis/Matrix/Spectrum.lean:141` —
   note **`Analysis/`**, not `LinearAlgebra/`) reads
   `A = conjStarAlgAut 𝕜 _ hA.eigenvectorUnitary (diagonal (RCLike.ofReal ∘ hA.eigenvalues))`.
   Unfold `conjStarAlgAut` to `U * D * star U`.
2. Build `e : (n → ℝ) ≃ₗ[ℝ] (n → ℝ)`, `e y = U *ᵥ y`, invertible from unitarity.
3. Prove the form identity
   `M.toQuadraticForm' (U *ᵥ y) = ∑ i, hM.eigenvalues i * y i ^ 2`
   via `(U*ᵥy) ⬝ᵥ (M *ᵥ (U*ᵥy))`, `star U * U = 1`, collapsing to `y ⬝ᵥ (D *ᵥ y)`.
4. Package as `QuadraticMap.IsometryEquiv`, hence `Equivalent`, then apply
   `sigPos_of_equiv_weightedSumSquares`.

Expected cost is **not** the mathematics (step 3 is three rewrites on paper) but the
coercion layer: `RCLike.ofReal` at `𝕜 = ℝ`, `star = id` on `ℝ`, `IsHermitian` vs `IsSymm`,
and unitary-submonoid membership unfolding. That is what will eat the build cycles.

Also still absent, and cheaper than the bridge: a **Lean-level non-vacuity control** for
this file. Non-vacuity is currently numpy/hand-verified only — `J_n` (all-ones) is
symmetric with `diag = 1 ≥ 0` and eigenvalues `(n,0,…,0)`, exactly one positive (`n=2..5`);
and `M = [[1,2],[2,1]]` has minor `−3 < 0`, so the conclusion is not the `i=j` identity.
The machine-checked version would be `sigPos (J 2).toQuadraticForm' = 1` via
`sigPos_add_finrank_le_of_nonpos` on `span {(1,−1)}`. **Recommended as the next LEAN
slot's target** — it is genuinely slot-sized, unlike the bridge.

## Predicates, separated

- `lem:minor` is **proved** (paper, `2026-10-04-c3-Q-lorentzian.tex` l.487–510).
- `lem:minor` is **formalised**, in the `sigPos` form of its hypothesis.
- The formalisation is **checked**: build exit 0, 0 sorries by the compiler's own warning,
  standard three axioms.
- The registry now **records** it, with the pointer graded by a controlled resolver.
- The eigenvalue-count form of the hypothesis is **not formalised**, in Lean or Mathlib,
  and is scoped above. It is a gap in coverage, not in the mathematics: the paper proof
  never counts an eigenvalue — it produces a 2-dimensional positive-definite subspace and
  contradicts a dimension bound, so `sigPos ≤ 1` is what it actually uses.
