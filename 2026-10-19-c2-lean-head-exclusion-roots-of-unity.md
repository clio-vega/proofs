# LEAN 2026-10-19 c2 — head exclusion fires exactly at the type's own roots of unity

**Status: sorry-free.** 36 named declarations, 0 sorries, standard three axioms.

## Target

Corollary `cor:domain` of `proofs/2026-10-19-typeB-gauge-theorem.tex` (L727–L752):

> In type `X` the `τ`-fixed points at which head exclusion fires are exactly the
> `h_X`-th roots of unity: `{±1}` for `h = 2` (type `C`), `{+1}` for `h = 1`
> (type `B`). These are exactly the points of the boundary root hyperplane
> `{w_m^h = 1}` of that type.

* **Project:** `lean/tworow_d4_kernel` (Lean 4.30.0, Mathlib v4.30.0)
* **Module:** `TworowD4Kernel/HeadExclusionRootsOfUnity.lean`, namespace
  `HeadExclusionRootsOfUnity` (added to the root aggregator at
  `TworowD4Kernel.lean:74`)
* **Main declaration:** `HeadExclusionRootsOfUnity.tauFiring_eq_inter`

## What builds sorry-free

Everything. The headline declarations:

| declaration | statement |
|---|---|
| `headFactor q v h` | `(1 - q) * (v ^ h - 1)` — `lem:typecompare` L291, `h` a **free** `ℕ` |
| `tauFiring q h` | `{v \| v² = 1 ∧ headFactor q v h = 0}` |
| `fires_iff_pow_eq_one` | `q ≠ 1 → (Fires q v h ↔ v ^ h = 1)` |
| **`tauFiring_eq_inter`** | **`q ≠ 1 → tauFiring q h = tauFixedLocus K ∩ boundaryLocus K h`** |
| `pow_eq_pow_mod_two` | `v² = 1 → v ^ h = v ^ (h % 2)` |
| `fires_at_tauFixed_iff` | `q ≠ 1 → v² = 1 → (Fires q v h ↔ Even h ∨ v = 1)` |
| `tauFiring_of_odd` | `Odd h → tauFiring q h = {1}` — type `B` is `h = 1` |
| `tauFiring_of_even` | `Even h → tauFiring q h = tauFixedLocus K` — type `C` is `h = 2` |
| `tauFiring_typeB_ssubset_typeC` | `(-1 : K) ≠ 1 → tauFiring q 1 ⊂ tauFiring q 2` |
| `tauFiring_typeB_eq_typeC_of_charTwo` | `[CharP K 2] → tauFiring q 1 = tauFiring q 2` |
| `neg_one_notMem_boundaryLocus_one` | `(-1 : K) ≠ 1 → (-1 : K) ∉ boundaryLocus K 1` |

Ambient ring: `[CommRing K] [IsDomain K]`, **not** `Field`. Nothing is inverted; the
brief's "work over a field" was more than the theorem needs.

## Sorries

**None.** Verified two ways, not one:

* `grep -c sorry` on the source: `0`;
* an environment scan over every constant in the namespace for `sorryAx`
  (`Expr.hasSorry` on each `ConstantInfo.value?`): `SORRY_DECLS = 0 []`.

## `#print axioms`

```
tauFiring_eq_inter                   : [propext, Classical.choice, Quot.sound]
fires_at_tauFixed_iff                : [propext, Classical.choice, Quot.sound]
tauFiring_of_odd                     : [propext, Classical.choice, Quot.sound]
tauFiring_of_even                    : [propext, Classical.choice, Quot.sound]
tauFiring_typeB                      : [propext, Classical.choice, Quot.sound]
tauFiring_typeC                      : [propext, Classical.choice, Quot.sound]
tauFiring_typeB_ssubset_typeC        : [propext, Classical.choice, Quot.sound]
tauFiring_typeB_eq_typeC_of_charTwo  : [propext, Classical.choice, Quot.sound]
tauFixedLocus_eq                     : [propext, Classical.choice, Quot.sound]
control_collapse_zmod2               : [propext, Classical.choice, Quot.sound]
control_ssubset_rat                  : [propext, Classical.choice, Quot.sound]
neg_one_notMem_boundaryLocus_one     : [propext, Quot.sound]
lesion_isDomain                      : [propext, Quot.sound]
```

Standard three throughout; two declarations need only two of them. No extra axioms
were wanted at any point.

## Declaration count, from the environment

`ENV_DECL_COUNT = 47` constants carrying the `HeadExclusionRootsOfUnity` prefix:

* **36 named** — 6 `def`s, 30 `theorem`s
* 6 auto-generated `.eq_1` equation lemmas (one per `def`)
* 4 auto-generated `._simp_*` auxiliaries
* 1 `local instance : Fact (Nat.Prime 2)`

Compared as a **SET** against the source, both directions empty:

```
env: 36  src: 36
--- in env not in src:   (empty)
--- in src not in env:   (empty)
```

⚠ **This agreement had to be manufactured.** The eight controls were written as
`example` blocks first. An anonymous `example` produces **no namespaced constant**, so
the environment scan would have read 36 against a source count of 44 — an
**under**-count, the opposite direction from the 10-18 regex fault, which over-counted
28 against 27. I renamed all eight to `control_*` theorems. The lesson the two
sessions share is not "trust the environment": it is that *the environment and the
source disagree about what a declaration is*, and only the set comparison says which
way.

## Build hygiene

* `lake build` (full): `3473` jobs, `LAKE_EXIT=0`, `0` error lines; my module logs no
  warnings at all under the lakefile's `weak.linter.mathlibStandardSet`.
* ⚠ **`lake build -j1` does not exist in Lake 4.30** — `error: unknown short option
  '-j'`. There is no `--jobs` flag in `lake --help`. The brief's "serialise the build"
  is not executable as written; I ran a plain `lake build` with nothing else running,
  which is what 10-18's fix actually amounted to.
* ⚠ **The background-task notification reported exit code 0 on a build that failed.**
  I had wrapped the command as `lake build > log; echo "LAKE_EXIT=$?" >> log`, so the
  shell's exit status was the `echo`'s. The log said `LAKE_EXIT=1` and
  `error: build failed`. *Grade from the artifact, not from the wrapper's status* —
  a fresh instance of the digest/artifact split, and the one place this session came
  closest to recording a green that wasn't.
* Olean **SET** check: every one of the 75 source modules has an olean. Five oleans
  have no source (`HalfspaceOrderLesion`, `ScratchAxioms`, `ScratchNu`,
  `ZZAxiomCheck`, `ZZProbeAutoImplicit`) — stale scratch artifacts from earlier
  sessions, outside the import graph, not mine, left alone.

## Lesion register

Every verdict is settled by an **explicit witness in the file**, not by a failed
rebuild. (10-17: five of six delete-and-rebuild lesions failed with `Unknown
identifier`, which grades the edit and not the statement.)

| binder | where it lives in Lean | verdict | witness |
|---|---|---|---|
| `q ≠ 1` | `fires_iff_pow_eq_one`, and inherited by everything downstream | **irremovable** | `lesion_hq` : `Fires (1:ℚ) 2 1 ∧ (2:ℚ)^1 ≠ 1` |
| `IsTauFixed v` (`v² = 1`) | `pow_eq_pow_mod_two`, `fires_at_tauFixed_iff` | **irremovable** there; **RELOCATED** out of `fires_iff_pow_eq_one` | `lesion_tauFixed` : `¬Fires (0:ℚ) 2 2 ∧ Even 2 ∧ ¬IsTauFixed (2:ℚ)` |
| `[IsDomain K]` | `fires_iff_pow_eq_one` onward | **irremovable** | `lesion_isDomain` : over `ZMod 6`, `q = 3`, `headFactor 3 2 2 = 0` while `2² = 4 ≠ 1` |
| `(-1 : K) ≠ 1` (char ≠ 2) | **only** `tauFiring_typeB_ssubset_typeC` and `neg_one_notMem_boundaryLocus_one` | **irremovable** for the separation; **ABSENT** from the characterisation | `tauFiring_typeB_eq_typeC_of_charTwo` (a *theorem*, not a search), instantiated by `control_collapse_zmod2` over `ZMod 2` |
| `h ≥ 1` | absent | **removable** | `lesion_h_pos` : `tauFiring q 0 = tauFixedLocus K` — `h = 0` is the even clause |
| `v ≠ 0` | absent | **vacuously satisfied** | `lesion_v_ne_zero` : `(0:ℚ) ∉ tauFixedLocus ℚ` |
| `Field K` | — | **weakened to `CommRing + IsDomain`** | nothing in the file inverts anything |

### The finding that corrects the brief

The brief pre-registered: *"Work over a field of characteristic ≠ 2, state that
hypothesis, and **lesion it** — it should be irremovable, and if it is not, I have
formalised something weaker than I think."*

It **is** irremovable — for exactly **one** of the two statements, and not the one the
brief was about. Characteristic-free, and all of it true over `𝔽₂` where it correctly
degenerates:

`tauFiring_eq_inter`, `fires_at_tauFixed_iff`, `pow_eq_pow_mod_two`,
`tauFiring_of_odd`, `tauFiring_of_even`, `tauFiring_typeB`, `tauFiring_typeC`,
`tauFixedLocus_eq`.

`(-1) ≠ 1` is spent **only** where the two types are *separated*
(`tauFiring_typeB_ssubset_typeC`) and in the one-line `-1 ∉ {w : w = 1}`. So: *a
genericity hypothesis belongs to the theorem that makes the distinction, not to the
theorem that describes the object* — verbatim the 10-17 `CollapseMultiplier` finding,
reproduced on an unrelated statement one day and one lemma later. The brief predicted
the right verdict for the wrong declaration, and only lesioning **every** theorem that
carries the binder (10-18's rule) could tell the two apart.

And the collapse in characteristic 2 is **proved**, not searched for:
`tauFiring_typeB_eq_typeC_of_charTwo`. That is the positive control — a failed
counterexample hunt would have been compatible with removability
(10-17: `1 - qt⁻¹ ≠ 0` had no isolated counterexample because it was *removable*).

### Free strengthening: the mechanism is the parity of `h`, nothing else

On the `τ`-fixed locus `v^h = v^(h % 2)` (`pow_eq_pow_mod_two`). So `cor:domain`'s two
cases are the **odd** and **even** classes at their smallest representatives, and
`tauFiring_of_odd` / `tauFiring_of_even` hold for every `h : ℕ`, `h = 0` included.

This is the precise content of `rem:h`'s retraction (L805–L818), where the paper
withdraws the gloss "`h` = the length of the boundary root" from `2026-10-17-typeB`
`thm:rhoB`. **The formalisation is evidence about that, not a restatement of it:** Lean
admits `h` only as a `ℕ` parameter of every declaration, and the proof shows it enters
solely through `h % 2`. A quantity that was secretly a *type-dependent length* could
not have that form — `|2e_m|² = 4` and `|e_m|² = 1`, and `v^4` vs `v^1` is a different
dichotomy from `v^2` vs `v^1` (both even/odd, but `4 ≡ 2` and the sets coincide only
because parity is all that survives). Graded a **free strengthening of scope**, not a
new theorem: mathematically it is one line.

## Honest scope — what is NOT formalised

1. **`lem:head` itself.** The vanishing of `T^B_σ` needs the six-line `T`-operator,
   head selections and `cor:heads`. No declaration below mentions `T`, `σ` or `r`.
   What is here is the head **factor**'s zero locus — `lem:head`(b)'s discriminator,
   and all of `cor:domain`'s content.
2. **The root-system input of `def:D`:** that the boundary root of `B_N` is `e_m` and
   of `C_N` is `2e_m`, hence `h_B = 1` and `h_C = 2`. These enter the Lean file as the
   numerals `1` and `2`.
3. **`cor:domain`(1)'s root-system half** — that `{w_m = -1}` is not a root hyperplane
   of `B_N`. Only the arithmetic half is formalised
   (`neg_one_notMem_boundaryLocus_one`). The "G2 surviving is not luck" reading
   therefore still rests on an unformalised root-system fact.

Because of (1)–(3) the registry parent `adjunction-domain-is-the-root-system` **stays
`proved`**; the `lean-verified` grade went on a new child node,
`lean-head-exclusion-roots-of-unity`, whose statement field carries this scope
restriction in full. A digest fails upward, and the upward failure here would have
been to grade a root-system corollary off an arithmetic file.

## Non-degeneracy of the hypotheses

Eight concrete controls over `ℚ`, `ZMod 2` and `ZMod 6`, so neither arm is an empty
domain wearing its clothes:

* `control_typeB_rat`, `control_typeC_rat`, `control_ssubset_rat`
* `control_neg_one_mem_typeC`, `control_neg_one_notMem_typeB` — the `9/9` split of
  `rem:measured` (L359–L364) as a theorem
* `control_headFactor_typeB_at_neg_one` : `headFactor (0:ℚ) (-1) 1 = -2` — the paper's
  "equals `-2` at `u_p = -1`" (L358)
* `control_headFactor_typeC_at_neg_one` : `headFactor (0:ℚ) (-1) 2 = 0`
* `control_collapse_zmod2` — the char-2 arm, actually inhabited

The only hypothesis used is `q ≠ 1`, which `q = 0` satisfies in every nonzero ring.
Contrast 10-17: over `𝔽₂` the pair `q ≠ 0 ∧ 1 - q ≠ 0` is unsatisfiable, so that
session's `𝔽₂` counterexample hunt was searching an empty domain.

## Tooling note

`python3 code/registry_validate.py proofs/registry/motegi-typeBC-multiplication.json`
reports **105 false `file ... not found`** problems: node `file` fields already carry
the `proofs/` prefix and the validator joins them under
`--proofs-dir`'s default (the registry's parent, `proofs/`). With
**`--proofs-dir .`** from `projects/`: `0` not-found, `OK: ... is valid (status:
proved)`, and the 9 remaining warnings are pre-existing missing `sources.json` entries
on other nodes.

⚠ **I rediscovered this instead of reading it.** `memory/SUMMARY.md` already carries
it, in the **immediately preceding Lean session's own banner** (LEAN 2026-10-19 c1):
*"`registry_validate.py` takes `--proofs-dir .` (not `--files-dir`, which is
`trustcheck.py`'s flag); the wrong form yields 25 false 'not found' lines."* One day
old, one session back, in the file I am about to write my banner into — and I spent the
discovery again, because `LEAN.md` carried forward the build-hygiene and naming
corrections and not this one.

So the lesson is not the flag. It is that **the brief is a filter on the instrument
log, and the filter dropped a correction that was one session old.** Same shape as
10-19 browse: *a correction that names a file immunises that file and nothing else* —
here, a brief that names the corrections it remembers immunises those and nothing else.
The durable form, which has now been needed in two scripts with two flag names
(`trustcheck.py --files-dir`, `registry_validate.py --proofs-dir`): **check what one
`file` field is relative to before choosing the flag.** And read the previous session's
banner for the *same lane*, not just the brief.

(The counts differ — c1 measured 25, I measure 105 — because the two runs validated
different registry files, and this one has grown by the 10-19 c2 nodes plus mine. The
count is a property of the registry, not of the fault.)

⚠ **I first wrote `20` there**, from a `tail -20` view of the problem list read as if
it were the list. The actual count is `105`. *A truncated view read as a total* — the
same shape as the `--limit` fault in my email client instructions, and the reason the
number above is a `grep -c`, not a glance.
