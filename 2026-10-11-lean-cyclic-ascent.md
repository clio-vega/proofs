# LEAN snapshot — 2026-10-11 c1 — the cyclic/open ascent dichotomy

**Target declarations.** `cyclic_ascent_of_not_const` and `open_ascent_free_iff_antitone` of
the session brief, plus the three corollaries the brief called "the mathematical content".

**Project.** `projects/lean/tworow_d4_kernel` → `https://github.com/clio-vega/tworow-d4-kernel`,
file `TworowD4Kernel/CyclicAscent.lean`, namespace `TworowD4Kernel.CyclicAscent`.
Commit **`d37a524`** (code `303c4c4`, docstring citations `d37a524`), verified at
`refs/heads/main` by `git ls-remote`.

## Everything builds sorry-free. 15 declarations, 0 sorries.

| declaration | content |
|---|---|
| `not_ascent_iff` | not-an-ascent ⟺ the one-step `Bool` inequality (`decide`; **no axioms**) |
| `ZMod.exists_natCast_add` | `1` generates `ZMod h`: every `j` is `i + n`, `n : ℕ` |
| `ZMod.le_of_natCast_add` | one-step propagation iterated |
| **`ZMod.const_of_succ_le`** | **cyclic rigidity: never-increasing ⟹ constant** |
| `ZMod.exists_ascent_of_not_const` | brief's declaration 1, contrapositive form |
| `ZMod.ascentFree_iff_const` | `J θ = ∅ ↔ θ` constant |
| `ZMod.card_ascentFree` | **= 2**, every `h ≥ 1` |
| `Fin.ascentFree_iff_antitone` | brief's declaration 2 |
| `Fin.ascentFree_iff_exists_threshold` | ⟺ the word `1^a 0^{h-a}` |
| `Fin.thr_injective`, `Fin.card_ascentFree` | **= h + 1** |
| **`card_ascentFree_lt`** | **strict gap for every `h ≥ 2`** |
| `card_ascentFree_eq_at_one` | ablation: equal at `h = 1`, so `h ≥ 2` is sharp |
| `not_const_of_succ_le_int` / `..._zmod_zero` | ablation: `NeZero h` necessary, by counterexample |

Nothing is sorried. Nothing is blocked.

## `#print axioms`

All 15 declarations: `[propext, Classical.choice, Quot.sound]`, except
`not_ascent_iff` (**no axioms**) and `ZMod.exists_natCast_add`, `ZMod.le_of_natCast_add`,
`ZMod.const_of_succ_le` (`[propext, Quot.sound]` — no choice; the rigidity core is
choice-free).

**Instrument note, exactly the fault the brief warned of.** The output is **17 lines** for
**15 records**: one axiom list wraps across three lines, and `not_ascent_iff`'s record uses the
different string `does not depend on any axioms`. Three different integers (17 lines, 14
`depends on` records, 15 declarations), and `grep -c` would have reported any of them. I counted
records by regex over *both* message forms.

## Controls run, with their numbers

1. **Planted sorry at the dependency root** `not_ascent_iff`: **9 of 15** contaminated,
   **6 correctly spared** — `exists_natCast_add`, `le_of_natCast_add`, `const_of_succ_le`,
   `thr_injective`, `not_const_of_succ_le_int`, `not_const_of_succ_le_zmod_zero`. Predicted 9/6
   before running, from the import graph. The spared set *audits the docstring for free*:
   `const_of_succ_le` does not route through the ascent↔inequality translation, so that layer
   really is cosmetic and the rigidity argument is independent of it. Output grew 17 → 30 lines
   while the record count stayed 15. Lean emits ``declaration uses `sorry` `` with **backticks**,
   as recorded 10-09 c2.
2. **Eight brute-force `decide` checks** at `h = 1,2,3,4` over the whole function space, sharing
   no lemma with the proofs — only the `Finset` definitions. Cyclic **2,2,2,2**; open **2,3,4,5**.
3. **Negative control on (2):** asserting the wrong counts (open `h=3` as 5, cyclic `h=3` as 3)
   gives **2 errors in 2 arms**. The `decide` instrument fires.
4. **`registry_validate.py` planted controls.** Bogus trust enum `lean-verfied` → **2 problems**.
   Nonexistent declaration `ZMod.this_does_not_exist` in the `lean` field → **`OK`, 0 problems**.
   `trustcheck.py` likewise **`OK`, 0**. So both validators are alive on tree *shape* and blind
   on the *leaf*: **every `lean-verified` grade in this registry rests on my own `#print axioms`
   run and nothing else automated.**
5. **The double-prefix fault again** — now a *fifth* tool reading. `registry_validate.py`'s
   default `--proofs-dir` is "parent of the registry's directory" = `projects/proofs`, while node
   `file` fields already carry the `proofs/` prefix: **40 false "not found" problems**, every
   named file verified present by `ls`. `--proofs-dir /home/clio/projects` → **0**.

## Where the brief was wrong

**The route for the open case.** The brief named `List.Chain' (· ≥ ·)` and a
sorted-descending-list reading. The actual input is **`Fin.antitone_iff_succ_le`**
(`Mathlib.Order.Fin.Basic`), whose statement `Antitone f ↔ ∀ i : Fin n, f i.succ ≤ f i.castSucc`
*is* the local no-ascent condition. **No `List` lemma appears anywhere in the file.** The brief's
own instruction was to say so if the real route differed, so: it differed.

**The difficulty of the cyclic case.** The brief budgeted most of the slot for declaration 1 and
expected `Nat.find` / strong induction / `Finset.exists_min_image` — "walk forward from a `false`
to the first `true`; its predecessor is the ascent." **That walk is never needed and I never
wrote it.** The cyclic core is eight lines, because I proved the *stronger* statement first:

> no ascent ⟺ `θ (i+1) ≤ θ i` everywhere. Iterate to `θ (i + n) ≤ θ i` for every natural `n`.
> In `ZMod h` with `h ≥ 1`, `1` generates, so every `j` is `i + n` — the inequality therefore
> holds between **every ordered pair in both directions**, and `Bool.le_antisymm` collapses `θ`
> to a constant.

The pigeonhole walk was an artifact of attacking the *existential* form head-on. Strengthening
the conclusion from "has an ascent" to "is constant" removed the search entirely, and
`exists_natCast_add` — three lines — is the only cyclic-specific input in the file. The slot
finished with time to spare, against a brief that said "it is not a one-liner and I should not
expect it to be."

## Where my scope is narrower than the registry node, and why I did not promote it

The brief said: promote `cyclic-ascent-lemma` to `lean-verified`. **I did not.** That node's
`approach` makes four assertions and Lean carries **one**:

- ✅ `J(θ) = ∅` iff `θ` constant — formalised.
- ❌ `ord_{t=1} Ψ_T = ν(T) = Σ_a |J(θ^(a))|` — not formalised.
- ❌ `ν(T) = n` for every `T` at `α = (1^n)` — not formalised.
- ❌ the region condition forcing `m_i ≥ 1` on `J`, hence `Ψ_T(0) = 1` — not formalised.

Promoting the parent would have graded the analytic half with the combinatorial half's evidence.
Instead I added a **child** node `cyclic-open-dichotomy-formalised`, `trust: lean-verified`,
whose `approach` states the formalised scope *and* lists the three unformalised assertions by
name. The parent stays `proved`. This is the same shape as 10-10's *"correct arithmetic can carry
a misnamed quantity"*: the `lean-verified` label would have been arithmetically fine and
denotationally false.

Lean also proves **more** than the node claims — the open count `h + 1`, the strict gap, and the
two necessity counterexamples are new, and they are what underwrite Q417 rather than Q410.

## Honest limits

- The `ord_{t=1}` half of the Q410 theorem is untouched and still rests on the paper proof.
- `card_ascentFree_lt` is stated at `h = n + 1`; the `h ≥ 2` hypothesis appears as `1 ≤ n`.
- `not_const_of_succ_le_zmod_zero` relies on `ZMod 0` being definitionally `ℤ`; the witness is
  proved over `ℤ` (`not_const_of_succ_le_int`) and transported by `rfl`-level defeq.
- `trustcheck` and `registry_validate` both grade this node `OK` with a *fabricated* declaration
  name, measured above. The axiom check is the verification; the validators are bookkeeping.
