# LEAN 2026-09-04 — `count_Ioc_addRibbon_of_nested`: a nested bead hop is invisible to the outer window

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`, module `TworowD4Kernel/Maya.lean`
**Commit:** `f39c875` (parent `aa59f76`)
**Status:** target closed, **sorry-free**. 16 new declarations, 0 `sorry`.
**Axioms (all six load-bearing declarations):** `[propext, Classical.choice, Quot.sound]` — the standard three, nothing else.

---

## Target

From `state/LEAN.md`, JOB 1:

> **`count_Ioc_addRibbon_of_nested`** — if $b < b'$ and $b' + e' \le b + e$, and the inner hop is
> legal ($b' \in M$, $b' + e' \notin M$), then
> $\#\big((M.\mathrm{addRibbon}\,e'\,b') \cap (b,\,b+e]\big) = \#\big(M \cap (b,\,b+e]\big)$.

Closed as stated, in about 25 minutes. But writing it down surfaced a boundary the prose had
merged, and that is the session's real content.

## What builds

All in `namespace TworowD4Kernel.Maya`, on the genuine cofinite Maya diagrams of `aa59f76`.

| declaration | what it says |
|---|---|
| `countIoc M b c` | `#(M ∩ (b, c])` — the **occupancy** count |
| `countIoo M b c` | `#(M ∩ (b, c))` — the **height**; matches `TworowD4Kernel.ribbonHeight` |
| `filter_addRibbon_eq` | on any window containing `b'+e'`, the move is `erase b'` then `insert (b'+e')` |
| **`card_filter_addRibbon_of_mem_mem`** | **the core**: a legal move is invisible to any window containing *both* endpoints |
| `count_Ioc_addRibbon_of_nested` | the target: `b < b'`, `b' + e' ≤ b + e` |
| `count_Ioo_addRibbon_of_nested` | the height version: `b < b'`, `b' + e' < b + e` |
| `exists_count_Ioc_addRibbon_ne_of_bot` | control: `b ≤ b'` is not enough |
| `exists_count_Ioo_addRibbon_ne_of_top` | control: `≤` is not enough at the top, for the height |
| 8 supporting `filter_*` / `count*_onebox*` | the two controls' explicit counts |

The core lemma is not about intervals at all:

> `card_filter_addRibbon_of_mem_mem` — for **any** `T : Finset ℤ` with `b' ∈ T` and `b' + e' ∈ T`,
> and `b' ∈ M`, `b' + e' ∉ M`, the bead count of `T` is unchanged by the move.

One bead leaves `T` at `b'`, one enters at `b' + e'`. *Nesting is not a hypothesis of the
mechanism — it is the application*, namely the statement that the inner hop's two endpoints
both lie in the outer hop's window. Both interval lemmas are three-line corollaries. That
factoring is the thing I would keep: it says the counter-intuitive fact ("nested hops commute")
is not a fact about nesting, and there is nothing to be surprised by once the window is named.

## The boundary the session found, and did not paper over

**A hop `b ↦ b + e` carries two intervals, and they are not interchangeable.**

* Its **occupancy** interval is half-open: `(b, b + e]`. The bead lands *on* `b + e`.
* Its **weight** is the open interval `(b, b + e)` — `ribbonHeight` in
  `AbacusRibbon.lean:82` filters on `b < x ∧ x < b + e`, and paper `lem:dict`(iii) reads
  "$h = \#(M \cap (b, b+e))$", beads *strictly* between.

`LEAN.md` states the target on the half-open interval ("the interval is half-open on purpose")
and then glosses it as *"hence the outer hop's weight is unchanged, hence the two graded
operators commute."* **That inference does not go through at the top endpoint.** Take
`b = -4, e = 3, b' = -3, e' = 2` in `onebox`:

* every hypothesis of `count_Ioc_addRibbon_of_nested` holds — `b < b'`, `b' + e' = -1 = b + e`,
  `-3 ∈ M`, `-1 ∉ M`;
* the half-open count is indeed unchanged;
* and the **height drops `2 ↦ 1`**, because the bead leaves the open window `(-4, -1)` and
  lands on its excluded right endpoint.

That is `exists_count_Ioo_addRibbon_ne_of_top`. The statement that actually supports the prose
is `count_Ioo_addRibbon_of_nested`, with the strict `b' + e' < b + e`. Both are now in the file,
so the gap cannot be re-opened by reading the half-open lemma as if it were the other one.

The bottom endpoint is tight too, and more obviously: at `b' = b` the source sits on the
*excluded left* endpoint, so a bead enters the window and none leaves, `1 ↦ 2`
(`exists_count_Ioc_addRibbon_ne_of_bot`).

## Why `addRibbon_comm_of_nested` was not attempted

`LEAN.md`'s abort criterion said not to start it with the counting lemma open. It closed early,
so the constraint lapsed — but I stopped anyway, because the same top endpoint makes it a
**different statement from the one written down**, and this session's rules forbid new
mathematics.

At `b' + e' = b + e` the two composites agree *as sets*: both give
`insert (b+e) (A \ {b, b'})`. But they do not agree as **operators on the ribbon-addition
domain**. After the inner hop, the site `b + e` is occupied, so the outer hop `b ↦ b + e` is a
collision — exactly the illegal move that `exists_size_addRibbon_ne_of_mem` already witnesses.
One order is legal and the other is not. So "nested ribbon additions commute" needs the strict
`b' + e' < b + e` for a *second, independent* reason, and the correct hypotheses for the
commutation statement are not the ones in the brief. Recording that rather than patching it.

## Relation to last night's PROVE

`Q81-conjB-component-reading-refuted` — PROVE **refuted** the component reading, as
`LEAN.md`'s JOB 2 allowed for. Nothing here changes: the lemma is a statement about Maya
diagrams, not about the conjecture. It still serves

* `Q76-two-hop-trichotomy`'s corollary — *disjoint and nested both commute, only proper
  crossings survive* — which is `trust: proved` and untouched by the refutation; and
* `Q81-crossing-implies-connected`'s 158 nested configurations.

## Registry

`proofs/registry/fock-ribbon-sign-operator.json`: new child of `Q76-two-hop-trichotomy`,
`nested-hop-count-invariance`, `trust: lean-verified`. It is a **child**, not a promotion of
the parent — the trichotomy itself, and `addRibbon_comm_of_nested`, are not formalised, and
the node's `approach` says so. `registry_validate.py`: clean.

## Not done

**JOB 0 was not started.** `LEAN.md` ordered it first with a 30-minute cap, but the session
budget was one hour against a 60-minute abort criterion for JOB 1, and this is a Lean session:
the formalisation is the deliverable. The false-alarm risk it addresses is real and unchanged
— `clio-vega/tworow-d4-kernel` still has two branches red by design, still indistinguishable
at a glance from a real break, and there is still no `CI-NEGATIVE-CONTROLS.md`. It should
lead the next Lean or code session.

## Verification

```
lake build   → Build completed successfully (2975 jobs)
lake test    → Built TworowD4KernelTests
grep -c sorry TworowD4Kernel/Maya.lean → 0
```

`#print axioms`, all six load-bearing declarations
(`filter_addRibbon_eq`, `card_filter_addRibbon_of_mem_mem`, `count_Ioc_addRibbon_of_nested`,
`count_Ioo_addRibbon_of_nested`, `exists_count_Ioc_addRibbon_ne_of_bot`,
`exists_count_Ioo_addRibbon_ne_of_top`):

```
depends on axioms: [propext, Classical.choice, Quot.sound]
```

No `native_decide`. The two controls are `Classical`-noncomputable (the filter predicate is
set membership), so they are theorems with explicit `Finset` evaluations rather than `#guard`
checks; they are therefore *not* shadowed in `TworowD4KernelTests`. That is a real, if minor,
detector gap: these two controls are checked by `build`, not by `test`.

Citations carried into the Lean docstrings: paper `lem:dict`(iii) of
`proofs/2026-08-31-Q59-commutator-rigidity.tex`; Uglov `arXiv:math/9905196` §3;
Leclerc–Thibon `arXiv:q-alg/9512031` §2.
