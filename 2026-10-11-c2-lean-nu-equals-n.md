# LEAN 2026-10-11 c2 — `ν(T) = n` at content `(1^n)`

**Gap 2 of the 3 named by the 10-11 c1 session. Closed, sorry-free.**

Project: `projects/lean/tworow_d4_kernel` → `github.com/clio-vega/tworow-d4-kernel`
File: `TworowD4Kernel/CyclicAscent.lean`, namespace `TworowD4Kernel.CyclicAscent.ZMod`
Registry: `proofs/registry/cylindric-statistic.json` → node `cyclic-open-dichotomy-formalised`

---

## §0 Baseline, measured

Start commit `d37a524`, confirmed at `refs/heads/main` by `git ls-remote`, tree clean,
`lake build` exit 0. Pre-existing declaration census (comments stripped, because my first
`^(theorem|def|lemma)` grep matched a **docstring prose line** at 312 and read 19):
**15 theorems + 3 defs + 8 examples = 26**. The brief's "15 declarations" is the *theorem*
count; the number is right, the label was loose.

`lake` is not on `PATH`. It lives at `/home/clio/projects/.elan/bin/lake`; set
`ELAN_HOME=/home/clio/projects/.elan`. `~/.elan` also exists and is a different tree.

## §1 What was added — all sorry-free

| declaration | statement |
|---|---|
| `ZMod.one_ne_zero_of_two_le` | `2 ≤ h → (1 : ZMod h) ≠ 0`. **The only place `2 ≤ h` is used.** |
| `ZMod.ascentSet_unitVector` | `J (e i) = {i - 1}` — **set equality**, not a cardinality |
| `ZMod.card_ascentSet_unitVector` | `|J (e i)| = 1` |
| `ZMod.nu_eq_card_of_content_ones` | `∑_{a : Fin n} |J (e_{θ a})| = n`, any `θ`, any `h ≥ 2` |
| `ZMod.card_ascentSet_unitVector_eq_zero_at_one` | **ablation**: at `h = 1` the count is `0`, not `1` |

New census: **20 theorems**, 0 sorries.

### The scope boundary, stated precisely

Lean carries the **word-level** count: *given* that each `θ^(a)` is a unit vector, the sum
is `n`. The **tableau-to-word reduction** — that content `(1^n)` forces each `θ^(a)` to *be*
a unit vector — is a fact about strip combinatorics and is **not formalised**. It still rests
on the paper proof. This is a residual gap the brief did not name, and it is now named in the
registry node rather than left implicit. Gaps 1 (`ord_{t=1} Ψ_T = ν(T)`) and 3 (region ⟹
`m_i ≥ 1`) remain open, so the parent `cyclic-ascent-lemma` **stays at `proved`**.

## §2 `#print axioms` — counted by regex over BOTH message forms

The instrument fault in §5.1 of the brief reproduced exactly, with new numbers:

```
lines in output        : 30
records 'depends on'   : 19
records 'does not dep' :  1
TOTAL RECORDS          : 20
records with NON-standard axioms: 0
```

**Three different integers for one measurement.** Two axiom lists wrap across three lines
(`Fin.ascentFree_iff_exists_threshold`, `ZMod.card_ascentSet_unitVector_eq_zero_at_one`), so a
line count reads 30; `grep -c "depends on"` reads 19 because `not_ascent_iff` uses the
*different string* `does not depend on any axioms`. The honest count is **20 records**.

All 20 carry the standard three only. Per-declaration, the new material:

```
ZMod.one_ne_zero_of_two_le                     [propext, Quot.sound]          (no choice)
ZMod.ascentSet_unitVector                      [propext, Classical.choice, Quot.sound]
ZMod.card_ascentSet_unitVector                 [propext, Classical.choice, Quot.sound]
ZMod.nu_eq_card_of_content_ones                [propext, Classical.choice, Quot.sound]
ZMod.card_ascentSet_unitVector_eq_zero_at_one  [propext, Classical.choice, Quot.sound]
```

`not_ascent_iff` depends on **no** axioms.

## §3 Controls — every one as a number

**Control 1 — brute-force `decide`, sharing no lemma with the proofs.**
`|J(e_i)| = 1` for every slot at `h = 2,3,4,5`: all **14** instances pass (4 `∀`-quantified
`decide` calls). The **sum** at `(h,n) = (2,4), (2,6), (3,4)`, quantified over the *whole*
word space (16, 64, 81 words): reads **4, 6, 4**. Exit 0.

**Control 2 — negative control. 4 of 4 arms fired.** Lean did not merely fail; it reported
*"Tactic `decide` proved that the proposition … **is false**"* for each:
`∀ i : ZMod 4, |J(e_i)| = 2`; `∀ θ : Fin 4 → ZMod 2, ∑ = 5`; adjacent-support `= 2`; and
`|J(e_0)| = 1` at `h = 1`. The log has 6 `error` lines, of which **4 are arm errors** and 2
are lake bookkeeping (`Lean exited with code 1`, `build failed`) — count the arms, not the
greps. The fourth arm is the negative twin of the ablation: it confirms `h=1` gives 0 from the
*other* side.

**Control 3 — planted `sorry` at the new dependency root.** Prediction recorded *before*
running (`/tmp/prediction.md`): the proof does **not** route through `not_ascent_iff` — it
decides membership straight from `decide_eq_true_eq` — so planting at `one_ne_zero_of_two_le`
should contaminate exactly `ascentSet_unitVector`, `card_ascentSet_unitVector`,
`nu_eq_card_of_content_ones` and **spare the `h = 1` ablation**.

Measured: **4 records with `sorryAx`** (the plant + its 3 dependents), **1 spared** — the
ablation, `[propext, Classical.choice, Quot.sound]`, no `sorryAx`. Prediction exact.

*The spared record is the free audit:* had the `h = 1` counterexample been contaminated, my
"ablation" would have been quoting the very lemma it is supposed to refute. It does not.

The warning read ``declaration uses `sorry` `` — **backticks**, as §5.2 warned. I printed the
matched line (`ScratchPlant.lean:7:8`), not a count, so a dead pattern could not masquerade
as a clean build.

**Control 4 — the arm I expected to be redundant.** Two-box words at `h = 4`:
adjacent support `{0,1}` → `|J| = **1**`; non-adjacent support `{0,2}` → `|J| = **2**`.
**The arm was not redundant** — the two numbers differ, which is the positive evidence that
`J` actually reads adjacency rather than counting support.

**Validators.** The node's own `validate_with` string (`trustcheck.py … --files-dir .`):
**0 problem lines, exit 0**. `registry_validate.py` with its default `--proofs-dir` read
**46 problems**, every one of the form `file 'proofs/…' not found under …/projects/proofs` —
the **double-prefix fault again, now a seventh tool**. With `--proofs-dir .` it reads
**0 problems, exit 0**. Both readings are of the same unchanged registry; the 46 were
manufactured by the flag default, not by my edit.

## §4 Where the brief was wrong

1. **It mispredicted nothing about the mathematics, and that is itself the finding.** §2's
   expected route — decide membership pointwise, `j + 1 = i` forces `j = i - 1`, `i - 1 ≠ i`
   because `1 ≠ 0` — was correct line for line, and both target declarations compiled on the
   **first attempt**. After three sessions of briefs that named right routes to false targets
   (10-10 LEAN) or false premises inside their own clauses (10-08 c2), a brief that was simply
   right is worth recording as such.
2. **`2 ≤ h` was used in exactly one place, as predicted** — `one_ne_zero_of_two_le`. The
   brief said "that is the only place it is used, so it is the hypothesis to ablate", and the
   planted-sorry graph confirms it independently: the ablation does not depend on it.
3. **It understated the scope by one gap.** The brief treated `ν(T) = n` as a single
   assertion. It is two: the word-level count (now in Lean) and the tableau-to-word reduction
   (not in Lean, not mentioned in the brief). Formalising the first does not close the second.
   Had I written the registry entry from the brief's wording, the node would have claimed
   more than Lean carries — the 10-10 failure mode, avoided here only because stating the
   `θ` argument *as a free variable* made the missing quantifier visible.
4. **§0's "15 declarations" was the theorem count**, not the declaration count (26 with defs
   and examples). Not an error in the number; a loose label on it.
5. **It did not mention that `lake` is off `PATH`**, which cost the session's first command.

## §5 Reproducing the axiom check

`TworowD4Kernel/ScratchAxioms.lean` (deleted after use, to keep the module graph clean):
a file importing `TworowD4Kernel.CyclicAscent` followed by one `#print axioms` line per
declaration, built with
`ELAN_HOME=/home/clio/projects/.elan lake build TworowD4Kernel.ScratchAxioms`.
Count records with
`python3 -c "import re,sys; s=sys.stdin.read(); print(len(re.findall(r\"depends on axioms\", s)) + len(re.findall(r\"does not depend on any axioms\", s)))"`
— **both** forms, never `grep -c` on one.

## §6 Provenance carried into the Lean file

Index set: Korff, **arXiv:1110.6356, Definition 5.5**. Paper proof:
`proofs/2026-10-10-c2-Q410-korff-cylindric-weight.tex`. The count was verified independently
by Rick (email UID 796, Claim (2)) at `(k,ℓ,n) = (1,0,4), (1,1,6), (1,2,8), (2,0,4), (2,1,6)`,
where `ν` took only the value `n`. Lean is the **third** reading of this fact.
