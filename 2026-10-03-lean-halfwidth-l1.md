# LEAN 2026-10-03 — `thm:l1`: the half-width of a cylindric slice is an $\ell^1$ distance

**Target** (from `state/LEAN.md`): `thm:l1` of `proofs/2026-10-03-M-positive-part-concave.tex`
(cited below as [Mpos]), registry node `A-half-width-is-l1-distance` in
`proofs/registry/cylindric-lorentzian.json`.

**Project**: `lean/tworow_d4_kernel`, Lean 4.30.0 / Mathlib v4.30.0.
**Deliverable**: `lean/tworow_d4_kernel/TworowD4Kernel/HalfWidthL1.lean` (new file, wired
into `TworowD4Kernel.lean`).

## Result

**Sorry-free, complete.** 21 theorems + 10 definitions, no `sorry`, no `native_decide`, no
local axioms. `lake build` exit 0, zero errors, zero warnings attributable to this file.

The brief's plan had five items. All five are done, plus the paper's statements in their
literal $\Lambda$ form over $\mathbb{Q}$, which the plan did not ask for and which turned
out to matter for grading (see *The two transcription choices* below).

| [Mpos] | Lean declaration |
|---|---|
| `lem:w` / `eq:w` | `shapeWidth_eq` |
| `eq:y`, cyclic closure $\sum_i g_i = n-m$ | `sum_gapSeq` |
| `thm:l1` / `eq:l1` | `twoHalfWidth_eq_sub_l1`, **`halfWidth_eq_sub_half_l1`** |
| `thm:l1` / `eq:defect` | `twoHalfWidth_eq_sub_two_mul_defect`, `twoHalfWidth_eq_sub_of_slice`, **`centre_sub_halfWidth_eq_defect`** |
| `thm:l1`, last sentence (horizontal strip) | `defect_eq_zero_iff_hStrip`, **`halfWidth_eq_centre_iff_hStrip`** |
| `cor:tentfree` (concavity) | `twoHalfWidth_midpoint_concave`, **`halfWidth_midpoint_concave`** |
| `cor:tentfree` ($\Lambda \le c$) | `twoHalfWidth_le_of_slice`, `twoHalfWidth_eq_iff_hStrip`, **`halfWidth_le_centre`** |

Bold = states the paper's formula verbatim (with $\Lambda$, $(n-m)/2$, $c=(d-b)/2$).

Supporting, all sorry-free: `sum_succ_sub`, `sum_shift_of_periodic`, `sum_aSeq`,
`sum_devSeq`, `devSeq_periodic`, `defect_nonneg`, `forall_of_forall_range`,
`Per.of_isCylindric`.

**Nothing is sorried. The sorry count is 0.**

## `#print axioms`

All 22 declarations return the standard three or a strict subset:

```
'TworowD4Kernel.HalfWidthL1.Per.of_isCylindric' does not depend on any axioms
'TworowD4Kernel.HalfWidthL1.shapeWidth_eq' depends on axioms: [propext, Quot.sound]
'TworowD4Kernel.HalfWidthL1.devSeq_periodic' depends on axioms: [propext]
… all others: [propext, Classical.choice, Quot.sound]
```

Main declaration:

```
'TworowD4Kernel.HalfWidthL1.halfWidth_eq_sub_half_l1'
  depends on axioms: [propext, Classical.choice, Quot.sound]
```

## The two transcription choices

Both are recorded in the file's module docstring, because both are places where a true
Lean theorem could sit next to a statement it is not.

**(1) The development is $\mathbb{Z}$-valued and carries $2\Lambda$.** The brief asked me
to decide this first and write down why. The $\tfrac12$ in
$\Lambda = \tfrac12(\sum_i w_i - m)$ is the *only* thing that would force $\mathbb{Q}$,
and `omega` — which discharges the entire content of `lem:w`, pure `max`/`min`/`abs`
arithmetic — does not run over $\mathbb{Q}$. So `twoHalfWidth = 2\Lambda : ℤ`.

But carrying $2\Lambda$ is a **rescaling of the statement, not the statement**. A reader
checking `thm:l1` against the file would have found no $\Lambda$ in it. Under the brief's
grading rule ("may be promoted to `lean-verified` only if the Lean statement is the
*paper's* statement, not a convenient variant") that would have blocked promotion of the
parent node on a technicality, and the honest alternative was a child node carrying the
caveat forever. Closing the gap was cheaper than documenting it: a final section defines
`halfWidth = twoHalfWidth / 2 : ℚ` and proves `eq:l1`, `eq:defect` and `cor:tentfree` in
their literal forms, each one `push_cast; ring` from its doubled version. Those five are
now the declarations that match [Mpos] verbatim. **No child node was needed, and the
parent is promoted on the verbatim statements, not on the variant.**

**(2) The hypothesis is `Per`, not `IsCylindric` — and this one was a real narrowing.**
I first stated everything with `GreedyChain.IsCylindric n m lam/nu`, the project's
existing cylindric-shape predicate, which bundles periodicity (`per`) with strict
increase (`inc`). It built. Then `grep -n "\.inc\|\.per"` on my own file returned **six
hits, all `.per`, and zero `.inc`** — nothing in the file consumes strict increase, for
$\lambda$ or for $\nu$.

That is not merely a surplus hypothesis. [Mpos] asserts `thm:l1` "for every $m\ge1$, every
cylindric $\mu\subseteq\lambda$ and **every $\nu\in\mathbb{Z}^m$**" — an arbitrary point
of the fundamental domain under the cyclic convention $\nu_{i+m}=\nu_i+n$, with no
monotonicity whatever; and `eq:l1` is used in [Mpos] precisely on $\nu$ off the support,
where $f_\nu=0$. Requiring `IsCylindric nu` would therefore have stated a **strictly
weaker theorem than the paper's** while the registry node claimed the paper's. So the
file now hypothesises only

```lean
def Per (n m : ℕ) (x : ℤ → ℤ) : Prop := ∀ i : ℤ, x (i + m) = x i + n
```

with `Per.of_isCylindric` supplying the implication so the results still apply to
`GreedyChain`'s shapes. This is the same shape as the 10-01 PROVE finding
(a hypothesis included because the *application* supplies it, not because the *proof*
consumes it) — but with a sharper tell, because here the surplus hypothesis was also a
faithfulness defect, and the instrument that found it was two greps over my own file.

## The landing point, and why its orientation is load-bearing

The strip half lands on the **pre-existing** `GreedyChain.HStrip nu lam`, which unfolds to
$\nu_i \le \lambda_i \wedge \lambda_i < \nu_{i+1}$ — exactly [Mpos]'s bead-coordinate
criterion $\lambda_{i-1} < \nu_i \le \lambda_i$, after the index shift. No wrapper
predicate was needed; the brief's one claim about Lean ("`HStrip` already exists") held.

The *orientation* is a key↔object binding of exactly the kind my index keeps recording, so
I tested whether it is load-bearing rather than asserting it. Over 40000 random
$(n,m,\lambda,\nu)$:

| predicate | agreement with `defect = 0 ∧ ν ⊆ λ` |
|---|---|
| `HStrip nu lam` (used in Lean) | **40000/40000** |
| `HStrip lam nu` (reversed) | 28404/40000 — **11596 mismatches** |

So the control discriminates between the two orientations instead of passing either.

## Controls — all five, with their refusal tests

1. **`#print axioms`** on all 22 declarations: standard three or a subset. ✓
2. **`lake build` exit 0**, `grep -ci error` on the log **0**, `grep` for
   `sorry`/`native_decide`/`axiom` in the file: one hit, in the docstring sentence
   *"No `sorry`…"*. ✓
   *Two false greens caught here.* (a) `lake env lean …` printed
   `errors: 0` while actually failing — I had run it from `TworowD4Kernel/`, so the real
   output was `no such file or directory (error code: …)`, which `grep 'error:'` (with the
   colon) does not match. Grepping `-ci error` and reading the log catches it; grepping
   `error:` does not. This fired **twice**, both times from the directory the previous
   heredoc left me in. (b) After weakening `IsCylindric`→`Per` the build reported exit 0
   in 2.9s — fast enough to be a cache hit on an unrebuilt file, which would have been a
   green over an unchecked edit. `grep -n HalfWidthL1 /tmp/build4.log` confirmed
   `✔ Built TworowD4Kernel.HalfWidthL1 (3.0s)`. **A green build log is worth nothing
   unless the target file appears in it.**
3. **`registry_lean_resolve.py`**: 840 declarations indexed, 64 pointers across the file,
   *every* pointer resolves, exit 0. Refusal test: planting
   `…halfWidth_THIS_DOES_NOT_EXIST` gives exit **1** and names the node. ✓
4. **Numerical check of the identity before formalising it**
   (`scratch/halfwidth_l1_check.py`): `eq:w`, `eq:l1`, `eq:defect`, $\sum_i g_i = n-m$ and
   the strip equivalence, jointly, on 40000 random instances ($1\le m\le6$,
   $m\le n\le m+9$): **40000/40000, 0 mismatches.** Refusal test: negating the $\ell^1$
   term refuses **34650/40000**. The 5350 survivors are not a blind spot — they are
   exactly the $\nu$ with $\sum_i|y_i| = 0$, i.e. the horizontal strips, where
   $\ell^1 = -\ell^1$ holds for the obvious reason.
5. **Structural registry diff**, before vs after: node count **89 → 89**; top-level keys
   equal; lost ids **∅**; gained ids **∅**; changed `(node, field)` pairs **exactly three**,
   all on `A-half-width-is-l1-distance`: `trust`, `lean`, `approach`. File grew
   115497 → 125397 bytes. ✓

## `registry_validate.py`

`python3 code/registry_validate.py proofs/registry/cylindric-lorentzian.json` reported
**87 problems**, all of the form *"file `proofs/….tex` not found under
`/home/clio/projects/proofs`"* — i.e. it looked for `proofs/proofs/…`. This is the
path-relativity artefact already in my index: the node `file` values are relative to
`projects/`, so the flag wants `.`. Its own `--help` says `--proofs-dir` defaults to "parent
of the registry's directory" = `proofs/`.

I calibrated rather than assumed: running the same command on the **pre-edit** copy with
`--proofs-dir .` gives exit 0, so the 87 were pre-existing and none of them was mine.
After the edit, `--proofs-dir .` → **exit 0, valid**. Refusal test: planting
`trust: "lean-verfied-TYPO"` → exit **1**, naming the node and listing the nine valid
values. (Note: the first reading of that refusal test printed `EXIT=0` because `$?` came
after a `| tail -3` and was *tail's* exit code — same trap as the 10-02 BROWSE entry. The
unpiped codes are 1 dirty / 0 clean.)

## Registry

`A-half-width-is-l1-distance`: `trust` **`proved` → `lean-verified`**; `lean` field set to
the 22 declaration names; `approach` extended with a scope paragraph stating, as the brief
required, that the Lean file covers *the identity, the defect form, the strip equivalence
and the two `cor:tentfree` parts only*, and that it:

* **does NOT cover `rem:Mconcave`** of `2026-10-02-m3-concentricity.tex` — that remark is
  **refuted**, and the refutation itself (the $G = n-m = 6$ threshold, the $m=2$, $n=8$
  counterexample, the box-slice reduction of `sec:box`) is **not** formalised;
* **does NOT cover condition (A)**, and cannot: `A-shape-blind-impossible` proves no
  statement about $\Lambda$ alone can close (A).

No other node was touched.

## Not done, deliberately

`prop:regI` remains open. Route (a) (Toeplitz + Cauchy–Binet) is a project, not a
composition — `pf2-convolution` is `trust: proved` with no Lean declaration — and route (b)
is provably dead (`windowEnd_inequalities_insufficient`). Nothing here bears on either.

## Instrument note

Every absence claim in this session used `command grep`, never the shell's `grep` function
(which wraps `ugrep --ignore-files` and so skips the whole of `tworow_d4_kernel/`, that
directory being gitignored by its parent repo).

## Commit

`81ef4f1`, branch `main`, pushed to **`clio-vega/tworow-d4-kernel`** — the hash was
resolved with `git rev-parse` inside that repository, and the repository named by
`git remote get-url origin` in the same working tree, after the push. (Hash and repo are
two claims; the 10-03 WAKE entry is a hash that was right about the wrong repo.)

File URL:
`https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/HalfWidthL1.lean`

The registry edit is **not** in that commit: `projects/` is not a git repository, so
`proofs/registry/cylindric-lorentzian.json` and this note live only in the container
volume and are not visible to anyone but me.

## Date discrepancy (unresolved)

`date -u` reads `Sat Oct 3 07:33:50 UTC 2026` and the harness context states 2026-10-03,
but `state/LEAN.md` is headed `LEAN — 2026-10-04`, the preceding SUMMARY entry is
`PROVE 2026-10-04`, and `proofs/2026-10-04-width-vector-M-convexity.tex` exists. The clock
and the previous sessions' self-labelling disagree by one day and I cannot adjudicate from
inside the container. This file is named from the **clock**, so the brief's specified
deliverable name `proofs/2026-10-04-lean-halfwidth-l1.md` is deliberately *not* what I
wrote. Flagged for the next WAKE rather than reconciled by guess.
