# LEAN 2026-10-20 c1 — `lem:order` + `cor:order` are Lean-verified

**Target.** `lem:order` (L261–L285) and `cor:order` (L286–L305) of
`proofs/2026-10-19-c2-sa-divisors-full-q.tex`.

**Project / module.** `lean/tworow_d4_kernel/TworowD4Kernel/MotegiOrderTailCount.lean`,
namespace `TworowD4Kernel.Motegi` — the *existing* namespace, extended, not forked.
Added to `TworowD4Kernel.lean` root imports (line 42).

**Status: sorry-free.** 0 sorries (environment scan and source grep agree).
Axioms exactly the standard three on all seven main declarations.

---

## 1. What builds

Everything. There are no sorries and no gaps; the paper proof was complete and
survived formalisation unmodified.

| tag | declaration | paper |
|---|---|---|
| — | `tail A k := (A.filter (k ≤ ·)).card` | `A_{≥k}`, L262 |
| — | `getD_le_iff_lt_countP` | the sorted-list fact underneath the whole lemma |
| — | `le_desc_iff_lt_tail` | ditto, transported to `Finset` |
| **T1** | `MLe_iff_tail` | `lem:order`, L261–269 |
| **T2** | `MLe_refl`, `MLe_trans`, `MLe_antisymm` | `cor:order`(1), L299–301 |
| **T2′** | `MLe_antisymm_asc` | same, *second mechanism* (see P3) |
| — | `tail_succ`, `mem_iff_tail_lt` | "the counting function determines the set", L300–301 |
| **T3** | `MLe_erase`, `downset_erase_closed` | `cor:order`(2), L293–294, L302 |
| **T3′** | `MLe_of_subset` | **free strengthening**, not in the paper |
| — | `tail_range`, `tail_sdiff_range` | `(Aᶜ)_{≥k} = (N−k+1) − A_{≥k}`, L303–305 |
| **T4** | `MLe_iff_compl` | `cor:order`(3), L295–297 |
| **T2/T4 bundled** | `MotegiBelow`, `instPartialOrderMotegiBelow`, `motegiCompl`, `motegiCompl_involutive`, `motegiComplOrderIso` | **strengthening**, see P2 |

## 2. Declaration census — from the environment, both directions

```
ENV_DECL_COUNT 40        (non-internal constants with module = MotegiOrderTailCount)
SORRY_COUNT    0
SRC_COUNT      35        (attribute-aware regex over the source)
```

`src \ env` = **∅**.  `env \ src` = **5**, all auto-generated:
`tail.eq_1`, `ascList.eq_1`, `desc.eq_1`, `descList.eq_1`,
`getD_le_iff_lt_countP.match_1_1`.  So 35 + 5 = 40, both directions closed.

⚠ **The source regex under-counted by 2 on its first run.** My first pattern was
anchored `^(theorem|lemma|def|instance|example)` and silently missed
`@[simp] lemma tail_range` and `@[simp] lemma tail_zero` — attribute-prefixed
declarations. Totals were 33 vs 40 and the gap *looked* like exactly the
auto-generated set plus noise; only `comm -13` naming the two missing lemmas
revealed that two of the seven were **real declarations my scan could not see**.
This is the 10-18/10-19 lesson firing a third time, and in the **under**-count
direction: ***only the symmetric difference, computed both ways and read as a
SET of names, says which way the error runs.*** Attribute-aware regex: 35,
symmetric difference exactly the 5 generated constants.

Anonymous `example`s in the file: **0**. All ten controls are named theorems, so
the environment scan can see them (10-19 lesson).

## 3. Axioms

```
MLe_iff_tail          : [propext, Classical.choice, Quot.sound]
MLe_antisymm          : [propext, Classical.choice, Quot.sound]
MLe_antisymm_asc      : [propext, Classical.choice, Quot.sound]
MLe_erase             : [propext, Classical.choice, Quot.sound]
MLe_iff_compl         : [propext, Classical.choice, Quot.sound]
motegiComplOrderIso   : [propext, Classical.choice, Quot.sound]
mem_iff_tail_lt       : [propext, Classical.choice, Quot.sound]
```

Standard three, nothing else, no `sorryAx`.

## 4. The structural finding: both directions of `lem:order` are one lemma

The paper proves `lem:order` by two separate index chases (L272–276 for ⇐,
L278–283 for ⇒). In Lean they are the **same iff applied at two different
points**:

> `getD_le_iff_lt_countP` — for a list sorted **decreasingly**, and `i` in range,
> `k ≤ l[i] ↔ i < l.countP (k ≤ ·)`.
>
> *The `(i+1)`-st largest entry is `≥ k` exactly when more than `i` entries are `≥ k`.*

- (⇒) is this iff at `i = (tail B k) − 1`, both sets, with `≤` transported by
  `MLe_iff_desc`.
- (⇐) is this iff at `k = desc B c`, both sets.

Nothing else about sortedness is used anywhere in the file. The paper's
"$(r-j+1)$-th largest of these" sentence (L275) and its "$a_{s-i+1},\dots,a_s$
are all $\ge k$" sentence (L283) are the two readings of that one biconditional.
The two index chases are not two arguments; they are one argument read forwards
and backwards. This is a presentational finding, not a correction — the paper is
right — but it is the reason the formalisation is 40 declarations and not 400.

The entry point is `MLe_iff_desc` from `MotegiOrder.lean` (10-17), so T1 is a
genuinely falsifiable theorem against a definition written in **Motegi's
increasing enumeration** (arXiv:2610.15344 Thm 2.6(i), L625–L636) by a session
three days ago that had never heard of tail counting. The 10-17 lesson —
***define it the source's way so the convenient reading is a theorem*** — paid.

## 5. Pre-registration — which arm fired

### P1 (T1 holds exactly as the tex states it) → **SPLIT: arm 1 ∧ arm 3**

Both, and they are compatible because they are claims about **different
domains**.

- **arm 1.** The tex's *criterion* is exactly right. `MLe_iff_tail` is the tex's
  displayed biconditional with no repair, no side condition, no extra hypothesis.
- **arm 3.** The tex's `∀ k ∈ [N]` and the Lean `∀ k : ℕ` are **different
  statements, and the tex's is strictly weaker over `Finset ℕ`** —
  `control_L1b_restrict_one_fails` is the witness: `B = {0}`, `A = ∅` satisfies
  every tail inequality for `k ≥ 1` and yet `B ≰ A`.

The tex is **not wrong**: in `[N] = {1,…,N}` no element is `< 1`, so its `k = 1`
instance *is* `tail B 0 ≤ tail A 0`, i.e. the cardinality clause `r ≤ s`. The
paper's own proof says so ("Taking $k=1$ gives $r\le s$", L272). Over `Finset ℕ`
that role relocates to `k = 0`, and the correct quantifier is `∀ k : ℕ`.

⇒ ***The cardinality clause `r ≤ s` is not a second condition; it is the tail
inequality at the smallest index the ambient admits.*** The two forms of
Motegi's order encode it in different places — one as an explicit conjunct, one
as the bottom instance of a quantifier — and that relocation is the content of
`lem:order`. Which index carries it is a **function of the ambient**, so a
statement of `lem:order` that does not name its ambient is incomplete.

### P2 (T4 is an iff in both directions with no extra hypothesis) → **arm 1, strengthened past arm 3**

`MLe_iff_compl` is an iff in both directions, no extra hypothesis beyond
`A, B ⊆ range N` (which is what makes `Aᶜ` mean anything).

Arm 3 asked whether "anti-automorphism" is the wrong word, with the detector
"try to state it as an `OrderIso` onto the order dual and see whether it
type-checks without a hack." **It type-checks.**

```lean
def motegiComplOrderIso (N : ℕ) : MotegiBelow N ≃o (MotegiBelow N)ᵒᵈ
```
sorry-free, standard axioms, where `MotegiBelow N := {A : Finset ℕ // A ⊆ range N}`
carries `MLe` as its `PartialOrder` (`instPartialOrderMotegiBelow`, which is T2
bundled). `map_rel_iff'` is literally `(MLe_iff_compl B.2 A.2).symm` — one term,
no hack. `motegiCompl_involutive` is `Finset.sdiff_sdiff_eq_self A.2`.

⇒ **The paper should say *order-reversing involution* and may say *`OrderIso`
onto the order dual*.** "Anti-automorphism" is defensible but weaker than what
is true: it does not record the involutivity, and it does not record that the
ambient must be **fixed** — `cor:order`(3) writes `Aᶜ` without naming `N`, and
`control_L4a_ambient_needed` shows the ambient is irremovable. This is a
strengthening in my own favour, so per the 10-18 lesson I state it as the Lean
declaration it is rather than as a banner: the theorem is
`motegiComplOrderIso`, and anyone can read its type.

⚠ **`Finset ℕ` cannot carry `MLe` as `≤`** — it already has inclusion. That is
why the bundled form lives on `MotegiBelow N`, not on `Finset ℕ`. The brief
anticipated this and said to name the theorems rather than fight the instance;
the subtype is not a fight, it is the right home, because `cor:order`(3) is only
*stateable* in a fixed ambient anyway. T2 therefore exists in both forms: three
named theorems on `Finset ℕ` and one `PartialOrder` instance on `MotegiBelow N`.

### P3 (antisymmetry needs the counting-function fact) → **arm 2**

**No.** `MLe_antisymm_asc` proves `MLe B A → MLe A B → B = A` from
`MLe`'s definition alone, with no `tail` anywhere:

`MLe B A` and `MLe A B` give `|B| = |A|`, so both cardinality offsets
`M.card − I.card` collapse to `0`, so the two increasing enumerations compare
entrywise in both directions, so `ascList B = ascList A` by
`List.ext_getElem`, so `B = A` by `Finset.sort_eq` + `Finset.val_injective`.

Not arm 3 — it is **not definitionally free** (unlike `thm:dict` on 10-17); it
needs list extensionality and the injectivity of `Finset.val`. But it does not
need the paper's route. Both are in the file: `MLe_antisymm` (paper's route, via
`mem_iff_tail_lt`) and `MLe_antisymm_asc` (sorted encoding). I kept the counting
fact because it is the paper's argument and because `mem_iff_tail_lt` is
reusable — it is the statement that `tail` is a *faithful* invariant, which is
what licenses reasoning about sets through their tail profiles at all.

⇒ ***Antisymmetry is cheap in whichever encoding you are already holding.*** The
paper reaches for the counting function because it has just proved T1 and is
standing in the tail picture; a session standing in the sorted picture never
needs it. The content of `cor:order`(1) is **reflexivity and transitivity** —
those really are "a conjunction of inequalities between integers" (L299) and
really are free. Antisymmetry is the only part with a mechanism, and it has two.

## 6. Lesion register — one named witness per verdict, all in the file

As pre-registered, a **mixed** register: the same binder is removable in one
theorem and irremovable in the adjacent one.

| # | lesion | pre-registered | **verdict** | witness in file |
|---|---|---|---|---|
| L1 | T1: restrict `∀ k` to `k ≥ 2` | FALSE | **FALSE — irremovable** | `control_L1a_restrict_two_fails`: `B={1}`, `A=∅`; all `k ≥ 2` tails are `0 ≤ 0`, yet `¬ MLe {1} ∅` |
| L1b | T1: restrict `∀ k` to `k ≥ 1` (the tex's own range) | *not pre-registered* | **FALSE over `Finset ℕ`** | `control_L1b_restrict_one_fails`: `B={0}`, `A=∅`. New: the tex's quantifier is the lesioned one, and it is safe only because `0 ∉ [N]` |
| L2 | T3: drop `m ∈ A` | REMOVABLE | **REMOVABLE** — confirmed | `control_L2_erase_notMem`: `7 ∉ {1,2}` and `MLe ({1,2}.erase 7) {1,2}`. `MLe_erase` is **stated without the hypothesis** |
| L3 | T1/T3: drop the ambient | REMOVABLE / never needed | **ABSENT** — never stated | `control_L3_no_ambient`: `MLe {1000} {1000, 10⁹}`, no `N` anywhere. `MLe_iff_tail` and `MLe_erase` quantify over all of `Finset ℕ` |
| L4a | T4: drop `A ⊆ range N` | (implied by arm 3) | **IRREMOVABLE** | `control_L4a_ambient_needed`: `MLe {0} {5}` holds, `¬ MLe (range 3 \ {5}) (range 3 \ {0})` — the complemented comparison fails when `A ⊄ range 3` |
| L4b | T4: vary `N` | **arm 3: free parameter** | **arm 3 CONFIRMED** | `control_L4b_ambient_free`: for any two `N, N'` both bounding `A ∪ B`, the two complemented comparisons are **equivalent** |

**L3 vs L4 is the sharp control and it came out as pre-registered.** One binder —
"an ambient `range N` containing the sets" — is *absent* from T1 and T3 and
*irremovable* from T4, with a concrete witness on each side, in adjacent
theorems in one file. Neither verdict is visible from the other: `MLe_iff_tail`
compiles with no `N` in scope, and `MLe_iff_compl` cannot even be **stated**
without one.

Per the 10-19 correction, L2 and L3 were graded by **locating the statement that
carries the hypothesis**, not by a failed build. L2: the paper *does* state
`A ∋ m` (L293, "for every $\bA\ni m$") and *does not use it* — genuinely
removable, not relocated. L3: no paper statement carries an ambient for T1 — the
tex's `k ∈ [N]` is a quantifier range, not a containment hypothesis, and that
range is exactly what L1b shows is doing work. So L1's binder and L3's binder
are the **same `N` playing two roles**, and only one of them is removable.

### A free strengthening the paper does not state

`cor:order`(2) is the one-element case of `MLe_of_subset`: **every** subset is
Motegi-below. `MLe_erase` is a one-line corollary of it, and `MLe_of_subset` is
the shortest route — `tail` is monotone in the set by
`Finset.filter_subset_filter`, so T1 gives it immediately. The paper's
"deleting an element decreases every $\bA_{\ge k}$ by $0$ or $1$" (L302) is
strictly more careful than necessary: the *amount* never matters, only the sign.

## 7. Scope honesty — what is NOT formalised

Only the order theory. The file contains no `w`, no `q`, no field, no root
system, no `Ã`, no transfer operator. In particular **not** formalised:

- `lem:inverse`, `lem:gamma`, `lem:gnz`, `thm:S` — everything the order is *for*.
- The use of `cor:order`(3) that motivated the target: that the support of
  Motegi's inverse formula (`Ã_{Mᶜ,νᶜ} = 0` unless `νᶜ ≤ Mᶜ`) agrees with the
  chain-expansion support of `lem:gamma`. T4 is the bridge; the two mechanisms
  it bridges are not here.
- That Motegi's `≤` is the order appearing in his Theorem 2.6(i) — that is
  `MotegiOrder.lean`'s docstring citation (arXiv:2610.15344 L625–L636), carried
  forward into this file's docstring, and it is a **reading of the source**, not
  a Lean theorem.

⇒ Parent `boundary-divisors-converse-full-q` stays `in-progress`; child
`tail-counting-order` stays `proved`; a **new `lean-verified` child** records
this file. ***Grade what the file proves, not what the file is about.***

## 8. Build discipline — what was actually checked

- **olean, not `$?`.** Deleted
  `.lake/build/lib/lean/TworowD4Kernel/MotegiOrderTailCount.olean`, then
  `lake build`: 3474 jobs, `✔ Built TworowD4Kernel.MotegiOrderTailCount (2.6s)`,
  olean present, 606216 bytes, dated 2026-10-20 07:35. No command was wrapped,
  so no exit status was trusted in either direction.
- **No competing build.** `.lake/*.lock` absent before starting; `pgrep -x lake`
  and `pgrep -x lean` both empty.
- **Liveness detector with a positive control.** `pgrep -x <name>` matches the
  exact process *name* and therefore cannot match my own prompt text, unlike
  `pgrep -f "lake|lean"` (the 10-19 fault). Three arms recorded: idle → `[]`,
  busy → `[894931 894932]`, idle → `[]`. The detector fires.
- **`lake build -j1` was not attempted** — it does not exist in this Lake. No
  serialisation was needed; nothing stalled.
- **Coverage graded from the olean and the environment**, never from the log.
  `lake` logged only the modules that warn (6 of 3474), which says nothing about
  the rest.

## 9. Session faults

1. **Source regex under-counted by 2** (§2). Caught by the symmetric
   difference, which is the only instrument that could have caught it; the
   totals-only reading (33 vs 40) was consistent with a clean result.
2. `Mathlib.Data.List.CountP` does not exist as an importable module in this
   Mathlib (v4.30.0) — the `countP` lemmas live in core/`List.Basic`. One
   probe, fixed by deleting the import.
3. Five first-round type errors, all mechanical (`Finset.Ico` needed
   `Mathlib.Order.Interval.Finset.Nat`; `Finset.card_sdiff` has the
   `∩`-signature here and `Finset.card_sdiff_of_subset` is the subset form;
   `rw [hk]` where `if_pos hk` was meant). None was a mathematical gap.
4. One real near-miss: the forward direction of `MLe_iff_compl` failed `omega`
   because I never introduced `h k` — I had the hypothesis `∀ k, …` in scope and
   the instance I needed was absent. `omega`'s counterexample report named the
   free variable (`d := ↑(tail B k)`), which is what located it. ⇒ **a
   quantified hypothesis in scope is not an instance of it**; `omega` reporting
   a *free* term is the signature.

## 10. Commit

Repository **`clio-vega/tworow-d4-kernel`** (not `clio-vega/lean` — see §11), commit
`002e6390b317e657c3017581c9794e83b42eeaef`, pushed `9fa1028..002e639`. Hash recorded in the registry node
`lean-motegi-order-is-tail-counting`.

---

## 11. Two tooling corrections the brief needs

1. ⚠ **The Lean project is not in `clio-vega/lean`.** The brief said "commit to
   `clio-vega/lean`". That repo's `.gitignore:6` is `tworow_d4_kernel/`, so
   `git status` in `projects/lean` reported a **clean tree** with the new file
   invisible — a silent, confident "nothing to commit". The Lean project is its
   own repository, `clio-vega/tworow-d4-kernel`, rooted at
   `projects/lean/tworow_d4_kernel`. Commit `002e639`, pushed
   `9fa1028..002e639`.
   ⇒ ***A clean `git status` is not evidence that the work is committed; it is
   also the signature of an ignored path.*** The discriminator is
   `git check-ignore -v <file>`, which names the rule and the line.
   Same shape as 10-17's `memory/for-robin/` is-not-a-git-repo: ***canonical ≠
   reachable***, and this time the containing repo *existed* and *was the wrong
   one*, which is harder to see than an absent repo.

2. ⚠ **`trustcheck.py` is unrunnable in this container.** It requires
   `--deployment <descriptor.json>` plus a subcommand
   (`validate|sources|report|ops`); `find /home/clio -name '*.json' -path
   '*deploy*'` returns **nothing**, and there is no `deployments/` directory
   anywhere under `/home/clio`. So the brief's prescription
   (`trustcheck … --files-dir .`) cannot execute — not because of the flag, but
   because its required input does not exist here. The `--files-dir` note in my
   standing instructions is correct *about the flag* and describes an invocation
   shape this copy of the script no longer accepts bare.
   The gate that **did** run is `registry_validate.py`, whose analogous flag is
   **`--proofs-dir`**, not `--files-dir`, and which exhibits the *same*
   double-prefix fault: its default (`parent of the registry's directory` =
   `proofs/`) against `file` fields that already carry the `proofs/` prefix
   gives **110** false `not found` lines; `--proofs-dir .` gives **0** and
   `OK: … is valid`.
   ⇒ ***A remembered flag is attached to a script, not to a task.*** Two
   different validators, one fault shape, two different flag names, and the
   remembered name belonged to the one that cannot run.
