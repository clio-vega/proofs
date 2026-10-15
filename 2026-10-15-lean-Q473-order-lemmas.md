# LEAN 2026-10-15 — Q473: the order lemma and the matching lemma, formalised

**Session:** LEAN slot, 2026-10-15.
**Paper proof:** `proofs/2026-10-15-Q473-GDMW-independence.tex`, `lem:order` (line 208) and
`lem:match` (line 222).
**External result formalised against:** Garbali–de Gier–Mead–Wheeler, *Symmetric functions
from the six-vertex model in half-space*, **arXiv:2312.14348** (source `Half_line_ASEP.tex`;
Conjecture `conj:F-orth`, line 443 — "subsets `κ ⊂ ℕ`", which is why `|κ|` is the *number of
parts* throughout and never the sum).

**Project:** `/home/clio/projects/lean/tworow_d4_kernel` (repo `clio-vega/tworow-d4-kernel`,
parent commit `67b890e`).
**New module:** `TworowD4Kernel/HalfspaceOrder.lean`, namespace `Q473`.

*Why this project and not a sibling:* the 6.2 G of built Mathlib oleans live under
`lean/tworow_d4_kernel/.lake/packages`. A sibling project would have to build Mathlib from
scratch, which does not fit a one-hour slot. The module is a leaf — nothing in the existing
library imports it, and it imports only Mathlib — so it cannot perturb the existing tree.

---

## 1. Result: both targets are sorry-free

| | |
|---|---|
| `lake` binary | `/home/clio/projects/.elan/bin/lake`, Lake 5.0.0-src+d024af0, Lean 4.30.0 |
| `lake build TworowD4Kernel.HalfspaceOrder` | real exit code **0**, 0 errors, 0 warnings |
| **DECL_COUNT** | **37** |
| **AXIOM-RECORD_COUNT** | **37** |
| sorries | **0** |
| `sorryAx` records | **0** |
| non-standard-axiom records | **0** |

The two counts are printed separately and asserted equal inside the audit itself
(`if keep.size != records then logError`). 37 = **22 hand-written declarations** + 15
auto-generated (`SP.mk`, `SP.rec`, `SP.parts`, `SP.chain`, `SP.pos`, `noConfusion*`,
`mk.injEq`, `mk.inj`, `mk.congr_simp`, `mk.sizeOf_spec`, `casesOn`, `recOn`, `ctorIdx`, …).

Records are enumerated **from the environment** (`env.header.moduleData[idx].constNames`),
not by a regex over source lines, so `@[simp]`, wrapped signatures and multi-line
`#print axioms` output cannot hide or split a declaration. Full record list:
`proofs/2026-10-15-lean-Q473-axioms.txt`. Audit source: `proofs/2026-10-15-lean-Q473-audit.lean`.

Axioms of the main declarations (from the record file):

```
Q473.wellFounded_of_finite_lowerSet :: propext, Classical.choice, Quot.sound
Q473.SP.finite_downset             :: propext, Classical.choice, Quot.sound
Q473.SP.instWellFoundedLT          :: propext, Classical.choice, Quot.sound
Q473.SP.not_wellFounded_gt         :: propext, Classical.choice, Quot.sound
Q473.exists_inj_iff                :: propext, Classical.choice, Quot.sound
Q473.card_le_of_inj                :: propext, Classical.choice, Quot.sound
```

Standard three only.

---

## 2. The statements, beside the informal ones

### Target 1 — `lem:order`

> **tex.** `⪯` is a partial order on `P` (strict partitions) with least element `∅`, and every
> principal down-set `P_{⪯κ} := {κ' ∈ P : κ' ⪯ κ}` is finite. Consequently `⪯` is well-founded.
> Here `κ' ⪯ κ :⟺ |κ'| ≤ |κ|` and `κ'_j ≤ κ_j (1 ≤ j ≤ |κ'|)`.

```lean
structure SP where                        -- `P`: strict partitions
  parts : List ℕ
  chain : parts.IsChain (· > ·)
  pos   : ∀ x ∈ parts, 0 < x

def SP.len  (κ : SP) : ℕ     := κ.parts.length        -- `|κ|`, the NUMBER OF PARTS
def SP.part (κ : SP) (j : ℕ) : ℕ := κ.parts.getD j 0  -- `κ_j`, 0-indexed, padded by 0

instance SP.instLE : LE SP :=                          -- `def:order`, verbatim
  ⟨fun κ' κ => κ'.len ≤ κ.len ∧ ∀ j < κ'.len, κ'.part j ≤ κ.part j⟩

instance SP.instPartialOrder : PartialOrder SP          -- partial order
instance SP.instOrderBot     : OrderBot SP              -- least element ∅
theorem  SP.finite_downset (κ : SP) : {κ' : SP | κ' ≤ κ}.Finite
instance SP.instWellFoundedLT : WellFoundedLT SP        -- well-founded
```

### The `§3` contrast — delivered

```lean
def     SP.single (m : ℕ) : SP := ⟨[m + 1], List.isChain_singleton _, by simp⟩
theorem SP.single_lt_single (m : ℕ) : single m < single (m + 1)
theorem SP.not_wellFounded_gt : ¬ WellFounded ((· > ·) : SP → SP → Prop)
```

`(1) ≺ (2) ≺ (3) ≺ ⋯`. This is the half of the Q473 headline that no computation in the
`.tex` verified, and it is the reason the triangularity had to be taken with respect to `⪯`
rather than `⪰`: it is now machine-checked that the choice was forced, not convenient.

### Target 2 — `lem:match`

> **tex.** Let `κ' = (κ'_1 > … > κ'_{n'})` and `c_1,…,c_N ∈ ℕ` with decreasing rearrangement
> `c_[1] ≥ … ≥ c_[N]`. There is an injection `ι : [n'] → [N]` with `κ'_j ≤ c_{ι(j)}` for all `j`
> iff `n' ≤ N` and `κ'_j ≤ c_[j]` for all `j ≤ n'`.

```lean
theorem card_le_of_inj {n N : ℕ} {ι : Fin n → Fin N}
    (hι : Function.Injective ι) : n ≤ N                       -- the `n' ≤ N` half

theorem exists_antitone_rearrangement {N : ℕ} (c : Fin N → ℕ) :
    ∃ σ : Equiv.Perm (Fin N), Antitone (c ∘ σ)                -- `c_[·]` exists

theorem exists_inj_iff {n N : ℕ} {a : Fin n → ℕ} {c d : Fin N → ℕ}
    (ha : Antitone a) (hd : Antitone d)
    (σ : Equiv.Perm (Fin N)) (hσ : ∀ m, d m = c (σ m)) (h : n ≤ N) :
    (∃ ι : Fin n → Fin N, Function.Injective ι ∧ ∀ j, a j ≤ c (ι j))
      ↔ ∀ j : Fin n, a j ≤ d (Fin.castLE h j)
```

`d` is "the decreasing rearrangement", spelled out as *antitone, and equal to `c ∘ σ` for a
permutation `σ`* — so the statement does not commit to one sort function, and
`exists_antitone_rearrangement` supplies the witness. `a` is `κ'`: only **antitonicity** is
assumed, so the Lean statement covers weakly decreasing `a`, a larger domain than the tex's
strictly decreasing `κ'`.

The `n' ≤ N` conjunct of the tex's right-hand side is split out as the hypothesis `h` plus
`card_le_of_inj`, rather than carried inside the iff: with `n' > N` the term
`d (Fin.castLE h j)` cannot be written at all, so the dependency has to be discharged first.
Together the two declarations are the tex's iff.

---

## 3. What Lean forced, that the `.tex` had asserted in a clause

**(a) The tex's `⟹` step is still unproved — and is no longer needed.** The `.tex` proves
`⟹` from *"pointwise domination is preserved by decreasing rearrangement"*, stated in one
clause and not proved there. Lean did not let me assert it, so `exists_inj_iff`'s forward
direction is instead a **counting argument**: if `a j > d_j`, then antitonicity of `d` puts
every index that can host `a j` strictly before `j`, so the hosting set `S` has
`|S| ≤ j`; but `σ⁻¹ ∘ ι` injects the `j+1` indices `j' ≤ j` into `S` (antitonicity of `a`
gives `a j ≤ a j' ≤ c (ι j')`), so `j + 1 ≤ |S| ≤ j`. Mathlib: `Fin.card_Iio`,
`Fin.card_Iic`, `Finset.card_le_card_of_injOn`, `Finset.card_le_card`.
**Honest status: the rearrangement-domination clause of the `.tex` remains unproved. The
lemma no longer rests on it.**

**(b) The length clause of `def:order` is redundant.** `SP.le_iff_forall` proves
`κ' ≤ κ ↔ ∀ j, κ'.part j ≤ κ.part j` — the pointwise inequalities **alone** force
`|κ'| ≤ |κ|`, because the parts are positive, so `κ'.part |κ| ≤ κ.part |κ| = 0` already says
`|κ'| ≤ |κ|`. The `.tex` states both clauses as if independent. (Positivity is doing the work
that the length clause appears to do; `SP.part_eq_zero_iff` is where it is spent.)

**(c) The strict-decrease hypothesis of `lem:order` is never spent — verified by lesion.**
Copy the module, **delete the `chain : parts.IsChain (· > ·)` field** and the two arguments
supplying it (`bot`, `single`), rename the namespace:

```
python3 - <<'PY'      # recipe, reproducible
s=open('TworowD4Kernel/HalfspaceOrder.lean').read()
s=s.replace('namespace Q473','namespace Q473Lesion').replace('Q473.SP','Q473Lesion.SP')
s=s.replace("  chain : parts.IsChain (· > ·)\n","")
s=s.replace("bot := ⟨[], List.isChain_nil, by simp⟩","bot := ⟨[], by simp⟩")
s=s.replace("def single (m : ℕ) : SP := ⟨[m + 1], List.isChain_singleton _, by simp⟩",
            "def single (m : ℕ) : SP := ⟨[m + 1], by simp⟩")
open('TworowD4Kernel/HalfspaceOrderLesion.lean','w').write(s)
PY
sed -i 's/^end Q473$/end Q473Lesion/' TworowD4Kernel/HalfspaceOrderLesion.lean
lake build TworowD4Kernel.HalfspaceOrderLesion     # exit 0, 0 errors
```

**All 22 hand-written declarations still compile with the field gone.** So `lem:order` holds
for *arbitrary lists of positive naturals*, not just strict partitions — the Lean statement is
strictly stronger than the tex's. The tex's own proof *does* use decrease ("since `κ'` is
decreasing every element of `κ'` lies in `[κ_1]`"); the Pi-type count below does not need it.
(The lesion module was deleted after the check; the recipe above reproduces it. First attempt
failed on `end Q473` not being renamed — a mechanical slip in the script, not a lesion result.)

**(d) The one-liner that hides a quantifier has no Mathlib lemma.** The tex ends `lem:order`
with *"a partial order all of whose principal down-sets are finite is well-founded"*. Mathlib
has `Set.Finite.wellFoundedOn` (well-founded **on** a finite set) but nothing global. Proved
here as its own named declaration — the valuable artifact:

```lean
theorem wellFounded_of_finite_lowerSet {α : Type*} {r : α → α → Prop} [IsStrictOrder α r]
    (h : ∀ a, {b | r b a}.Finite) : WellFounded r
```

via the down-set cardinality as a `ℕ`-valued measure: `r b a` makes `D(b) ⊂ D(a)` (⊆ by
transitivity, strict because `b ∈ D(a) \ D(b)` by irreflexivity), so
`Finset.card_lt_card` + `Subrelation.wf` + `InvImage.wf wellFounded_lt`. **The hidden
quantifier is the induction on the cardinality, and it needs transitivity and irreflexivity
— i.e. `[IsStrictOrder]`, not merely "a partial order".**

**(e) The down-set finiteness count changed mechanism.** The tex bounds the down-set using
`κ'_1 ≤ κ_1` and monotonicity. Lean route: inject into `∀ i : Fin |κ|, Fin (κ_i + 1)`, finite
by `Pi.finite`, injective by `SP.eq_of_part_eq`. Needs only the pointwise bound — which is
why (c) holds.

## Mathlib lemmas used

`List.getD_eq_getElem`, `List.getD_eq_default`, `List.getElem_mem`, `List.ext_getElem`,
`List.isChain_nil`, `List.isChain_singleton`, `Pi.finite`, `Finite.of_injective`,
`Set.toFinite`, `Set.Finite.subset`, `Set.Finite.mem_toFinset`, `Finset.card_lt_card`,
`Finset.ssubset_iff_of_subset`, `Finset.card_le_card`, `Finset.card_le_card_of_injOn`,
`Fin.card_Iio`, `Fin.card_Iic`, `Fin.castLE_injective`, `Fintype.card_le_of_injective`,
`Subrelation.wf`, `InvImage.wf`, `wellFounded_lt`, `RelEmbedding.natGT`,
`RelEmbedding.not_wellFounded`, `Tuple.sort`, `Tuple.monotone_sort`, `irrefl_of`.

---

## 4. Controls

**Plant (contamination), and it reached the control.** The *same* `audit` function, run over a
second domain: three declarations in namespace `Q473Plant` — `finite_downset_plant`
(statement of `finite_downset`, proof replaced by `sorry`), `wf_plant` (uses the plant, so it
tests **propagation**, not just direct detection), and `clean_control` (a real proof in the
same domain).

```
PLANT-RECORD 1 :: thm :: Q473Plant.clean_control         :: propext,Classical.choice,Quot.sound
PLANT-RECORD 2 :: thm :: Q473Plant.wf_plant              :: propext,sorryAx,…  <<<SORRYAX
PLANT-RECORD 3 :: thm :: Q473Plant.finite_downset_plant  :: propext,sorryAx    <<<SORRYAX
PLANT-DECL_COUNT 3 / PLANT-RECORD_COUNT 3 / PLANT-SORRYAX_RECORDS 2
```

**2 of 2 fired, direct and propagated; the clean control in the same domain was spared.** So
the audit is reporting contamination, not the domain.

**Replay check (§1.4), and it caught a replay.** `touch`-ing the source and rebuilding gave
`exit 0` with the **olean mtime unchanged** (`22:07:37` before and after) — a trace-hash
match, not an elaboration. Deleting
`.lake/build/lib/lean/TworowD4Kernel/HalfspaceOrder.olean` and rebuilding printed
`✔ [801/801] Built TworowD4Kernel.HalfspaceOrder (2.4s)` and regenerated the olean at
`22:09:40`. The exit code was reported separately from the pipeline (`PIPESTATUS[0]`), not
from a trailing `echo`.

**A scan that silently measured nothing, and why it is not in the evidence.** I first tried to
decide (c) by scanning proof terms for the constant `Q473.SP.chain`. It returned **0 users —
and also 0 users of `Q473.SP.pos`, which the proofs demonstrably do spend.** Structure
projections are `Expr.proj` nodes, not constant applications, so `getUsedConstants` cannot see
them: the scan was blind over its whole domain and its "0" was unreadable, not clean. The
positivity control is what exposed it. Replaced by the lesion in (c), which fails loudly.

---

## 5. Scope — what is NOT carried

* `thm:tri`, `thm:main` — need `K = Q(q,a,Y)` rational-function arithmetic and the zero/pole
  structure of `φ_m`. Untouched.
* `thm:top` — rests on **Hypothesis H**, the unproved analytic interchange (`sec:gaps` item 2).
  Deliberately untouched: formalising a statement whose hypothesis is open would manufacture a
  `lean-verified` grade for something that is not proved.
* `lem:phi` (zeros and poles of `φ_m`) — untouched.
* Nothing here mentions `G_κ`, `F_κ`, or the half-space model. These are the two *order and
  combinatorics* steps, in the abstract, and that is the whole claim.

The registry nodes `evaluation-triangularity` and `topological-independence` stay at `proved`:
they assert much more than these two lemmas. Two **new** nodes are added rather than
overwriting them.

## 6. Instruction-fault report — one new instance, in the FLAGLESS default

I ran `registry_validate.py` in exactly the form the LEAN skill prescribes —
`python3 code/registry_validate.py proofs/registry/<name>.json` from `/home/clio/projects`,
**no flag** — and got **13 of 15 "problems" that were false**, every one of the form

```
root/<node>: file 'proofs/2026-10-15-Q473-GDMW-independence.tex' not found under /home/clio/projects/proofs
```

`--proofs-dir` **defaults to the parent of the registry's directory** (`proofs/registry/..` =
`proofs`), and node `file` fields already carry the `proofs/` prefix, so the default joins to
`proofs/proofs/...`. This is the **same double-prefix fault** recorded for `trustcheck.py`
(`--files-dir proofs`), `emit.py` (`--files-dir registry`) and `registry_validate.py`
(`--proofs-dir .` noted for a different flag) — but here it is in the **flagless default**,
so there is no wrong flag to omit: following the instruction *is* reproducing the fault.
11 of the 13 were on pre-existing nodes, so the fault predates this session and has been
sitting in a registry that reads as 11-problems-dirty under its own documented invocation.

**Correct form, from `/home/clio/projects`:**

```
python3 code/registry_validate.py --proofs-dir . proofs/registry/<name>.json
```

With that: `OK: proofs/registry/halfspace-sixvertex-orthogonality.json is valid
(status: proved)`, exit 0, **zero** problems.

The other 2 of the 15 were real and are fixed: both new nodes were missing the `children` key.
Reported here, not patched into CLAUDE.md or the skill (read-only).
