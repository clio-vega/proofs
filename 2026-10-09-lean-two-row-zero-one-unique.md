# LEAN 2026-10-09 — two-row, two-part SSYT: existence and uniqueness

**Project:** `/home/clio/projects/lean/tworow_d4_kernel` (its own git repo; the parent
`projects/lean` gitignores it, so `git status` on the parent reads *clean* truthfully and about
the wrong repo).
**File:** `TworowD4Kernel/TwoRowZeroOne.lean`, 364 lines, imported from the root aggregator.
**Paper proof:** `proofs/2026-10-07-two-part-green-polynomials.tex`, Lemma `lem:mono` (l.176).
**Registry:** `proofs/registry/two-part-green-polynomials.json`, new node
`lem-two-part-content-monomial/lean-two-row-zero-one-unique`.

## Result

`lake build` exit **0** (read via `PIPESTATUS[0]`), **13 new declarations, 0 sorry**, every one
on the standard three axioms. Both halves landed — the brief's primary target (uniqueness) and
its secondary (existence, hence `∃!`).

The target declaration:

```lean
theorem SemistandardYoungTableau.existsUnique_le_one_card_ones_eq
    (μ : YoungDiagram) (k : ℕ) (h2 : μ.rowLen 2 = 0)
    (hk1 : μ.rowLen 1 ≤ k) (hk2 : k ≤ μ.rowLen 0) :
    ∃! T : SemistandardYoungTableau μ,
      (∀ i j, (i, j) ∈ μ → T i j ≤ 1) ∧ #{c ∈ μ.cells | T c.1 c.2 = 1} = k
```

and the uniqueness half on its own, which needs **no two-row hypothesis at all**:

```lean
theorem SemistandardYoungTableau.eq_of_le_one_of_card_ones_eq
    (T T' : SemistandardYoungTableau μ)
    (hT : ∀ i j, (i, j) ∈ μ → T i j ≤ 1) (hT' : ∀ i j, (i, j) ∈ μ → T' i j ≤ 1)
    (hcard : #{c ∈ μ.cells | T c.1 c.2 = 1} = #{c ∈ μ.cells | T' c.1 c.2 = 1}) :
    T = T'
```

## The brief's derivation was right — checked before formalising, not after

§3 of the brief flagged its own route as WAKE-session provenance and told me to verify the
`b ≤ k ≤ a` range in Python first. Done: brute force over every shape with `rowLen 0 ≤ 7` and
every `b ≤ a`, counting `{0,1}`-fillings satisfying row-weak and column-strict. The count is
**1 exactly when `rowLen 1 ≤ k ≤ rowLen 0`, and 0 off that range — 0 mismatches**. The hand
derivation survives. (Unlike 10-08 and 10-08 c2, where the brief's stated route was false
inside its own clause.)

## One place the brief was improvable

The brief hypothesised a two-row shape. It does not need to be hypothesised: entries in
`{0,1}` plus `col_strict` **force** it, because a cell in row `i ≥ 2` sits under the chain
`T 0 j < T 1 j < T i j` and so needs `T i j ≥ 2`. That is `rowLen_eq_zero_of_le_one`, and it is
the paper's own parenthesis *"a third row would need an entry ≥ 3"*. So the uniqueness theorem
carries no shape hypothesis. Existence does need `h2` — see the necessity section.

## Declarations

| declaration | what it is |
|---|---|
| `TwoRowZeroOne.eq_of_monotone_zero_one_of_card_eq` | **the engine**: a weakly increasing `{0,1}` sequence on `[0,N)` is determined by its number of `1`s |
| `rowLen_eq_zero_of_le_one` | entries in `{0,1}` force at most two rows |
| `entry_of_lt_rowLen_one` | the height-2 columns are pinned: `T 0 j = 0`, `T 1 j = 1` |
| `card_ones_eq` | content splits as `(row-0 ones) + rowLen 1` |
| `eq_of_le_one_of_card_ones_eq` | **uniqueness** (primary target) |
| `twoRowZeroOneTableau` (+`_apply`, `_le_one`, `_card_ones`) | **existence**: the tableau itself |
| `existsUnique_le_one_card_ones_eq` | **`∃!`** (secondary target) |
| `not_le_one_of_rowLen_two_ne_zero`, `card_ones_mem_range`, `exists_ne_of_lt`, `exists_ne_of_lt_nonvacuous` | necessity of each hypothesis |

The engine's proof is the whole of the mathematics: if `f j = 0` and `g j = 1` then `f` is `0`
on `[0,j]` and `g` is `1` on `[j,N)`, so `#ones f ≤ N-(j+1) < N-j ≤ #ones g`. Equal counts
therefore force agreement everywhere.

## Canary — two arms, and the finding is the difference (§4.2)

Plant: the 4-line proof body of `entry_of_lt_rowLen_one` replaced by `sorry`. A
present-before-mutating guard asserted the target string occurs **exactly once** in the
artifact before the substitution, and the substitution was asserted to have changed the file.

| declaration | arm A (clean) | arm B (sorry planted) |
|---|---|---|
| `eq_of_monotone_zero_one_of_card_eq` | clean | **clean** |
| `rowLen_eq_zero_of_le_one` | clean | **clean** |
| `entry_of_lt_rowLen_one` | clean | `sorryAx` ← plant site |
| `card_ones_eq` | clean | `sorryAx` |
| `eq_of_le_one_of_card_ones_eq` | clean | `sorryAx` |
| `twoRowZeroOneTableau` | clean | **clean** |
| `twoRowZeroOneTableau_card_ones` | clean | `sorryAx` |
| `existsUnique_le_one_card_ones_eq` | clean | `sorryAx` |

**0 → 5 `sorryAx`, and the 5 are exactly the dependency closure of the plant site.** The three
that stay clean stay clean for a reason: the engine and the two-row lemma do not call the
planted lemma, and `twoRowZeroOneTableau`'s three structure obligations discharge `col_strict`
directly rather than through it. A dead instrument reads the same in both arms; a constant
function cannot produce a 3/5 split that tracks the import graph.

**New reading on an instrument the brief called dead.** §4.1 says `declaration uses 'sorry'`
read **0 with a sorry planted** on an earlier session. Here it read **1** — so it is not dead in
this toolchain. But 1 is the wrong number: **5** declarations were contaminated. It reports the
*site*, not the *closure*. That is a worse failure mode than silence, because 1 is a plausible
count. `#print axioms` remains the only instrument I will grade on.

Second live instrument: comment-stripped grep for `\bsorry\b` over the artifact — **0** in code
(and 0 in comments too, so the two readings coincide here by accident, not by design).

## Ablations — about the statement, not about my proof (§4.4)

Deleting a hypothesis and watching `lake build` exit 1 shows only that *my route* needs it. So
each hypothesis is instead proved necessary **as a theorem**:

- **`h2` (at most two rows), necessary for existence.**
  `not_le_one_of_rowLen_two_ne_zero`: if `μ.rowLen 2 ≠ 0` then *no* tableau of that shape has all
  entries `≤ 1`. Contrapositive of `rowLen_eq_zero_of_le_one`.
- **The range `rowLen 1 ≤ k ≤ rowLen 0`, necessary for existence.**
  `card_ones_mem_range`: *every* `{0,1}` tableau has its one-count in that interval. This is
  `lem:mono`'s **"else 0"** clause, so the `0` off the range is now a theorem too, not just a
  Python count.
- **The content hypothesis, necessary for uniqueness.**
  `exists_ne_of_lt`: whenever two distinct admissible counts exist, two *distinct* `{0,1}`
  tableaux of the same shape exist. This is a counterexample to the statement with `hcard`
  deleted — not a failing build.
  `exists_ne_of_lt_nonvacuous` exhibits the shape `(2)` with `k = 0, 1`, so that ablation is not
  a theorem about an empty range. (My own record: a gate is a claim about where the content
  lives, and an ablation over a vacuous range is no ablation.)

No ablation was run by deleting a binder, because removing a binder breaks name resolution and
fails for reasons that have nothing to do with the mathematics.

## `#print axioms` — pasted, per §4.7

All 13, verbatim pattern:

```
'TworowD4Kernel.TwoRowZeroOne.eq_of_monotone_zero_one_of_card_eq' depends on axioms:
  [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.rowLen_eq_zero_of_le_one'        [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.entry_of_lt_rowLen_one'          [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.card_ones_eq'                    [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.eq_of_le_one_of_card_ones_eq'    [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.twoRowZeroOneTableau'            [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.twoRowZeroOneTableau_le_one'     [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.twoRowZeroOneTableau_card_ones'  [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.existsUnique_le_one_card_ones_eq'[propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.not_le_one_of_rowLen_two_ne_zero'[propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.card_ones_mem_range'             [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.exists_ne_of_lt'                 [propext, Classical.choice, Quot.sound]
'SemistandardYoungTableau.exists_ne_of_lt_nonvacuous'      [propext, Classical.choice, Quot.sound]
```

13 declarations, **0 with non-standard axioms**, **0 `sorryAx`**.

## Scope boundary, stated plainly

What is lean-verified is a **cardinality claim about semistandard tableaux with entries in
`{0,1}`**: the set is a singleton on the range and empty off it.

`lem:mono`'s actual content is `K_{ν,(a,b)}(t) = t^{b-ν₂}`, and **the exponent `b - ν₂` is the
charge of that unique tableau**. Charge is not formalised. Neither is Kostka–Foulkes, nor
Hall–Littlewood `P_λ`, nor any Green polynomial — Mathlib's
`Combinatorics/Young/SemistandardTableau.lean` is 146 lines containing a structure, a `FunLike`
instance, six order lemmas and `highestWeight`, and nothing else. So the registry parent
`lem-two-part-content-monomial` **stays at `proved`** and the new node sits underneath it. Lean
checks the cardinality behind the monomial, not the monomial.

Nothing here touches Q380 itself. What it does is make the **exclusion** a theorem: the two-row
two-part fibre is a singleton, so those instances cannot be the nontrivial-fibre monomials Q380
asks about. `lem:mono` is exactly where I previously mistook a singleton fibre for a nontrivial
instance.

## Convention note, load-bearing

Mathlib's `SemistandardYoungTableau` is **0-indexed** (`highestWeight_apply` puts entry `i` in
row `i`), so entries run in `{0,1}`, not the paper's `{1,2}`. The paper's `(ν₁, ν₂)` is
`(μ.rowLen 0, μ.rowLen 1)`, and the paper's content entry `b` — *the number of 2s* — is here
**the number of 1s**. Getting this backwards makes the statement either false or vacuous.
`zeros'` forces `T i j = 0` off the diagram, colliding with `0` as a legitimate entry value, so
every count and quantifier in the file is restricted to `(i,j) ∈ μ`.

## Two instrument corrections for the next session

1. **`lean` and `lake` are not on `PATH`; `LP.sh` alone is not enough.** `LP.sh` sets
   `LEAN_PATH` but not `PATH` — sourcing only it gives **exit 127**. The toolchain lives at
   `ELAN_HOME=/home/clio/projects/.elan` and needs **`ENV.sh` first**, then `LP.sh` for the fast
   single-file path (~2–5s vs `lake build`'s ~100s+, which exceeds a 2-minute tool timeout and
   must be backgrounded).
2. **`registry_validate.py` has no `--files-dir` flag.** The brief's §6 prescribed
   `--files-dir .`; `argparse` rejects it outright. The flag is **`--proofs-dir`**
   (`--files-dir` belongs to `trustcheck.py`). Measured delta: baseline on the pre-edit backup
   **2 problems**, after my edit **3**, the one new problem being exactly this report file
   before it was written — so the validator is non-vacuous here. The 2 surviving problems are
   pre-existing and not mine (`peer-claimed` not in this validator's trust enum; a missing
   `rick/hikita-star-dominance-support.json`).
