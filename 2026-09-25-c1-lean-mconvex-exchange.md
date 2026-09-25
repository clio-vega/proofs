# LEAN 2026-09-25 c1 — the M-convex exchange axiom

**Target declaration.** `MConvexExchange.insupp_exchange`
**Project.** `projects/lean/tworow_d4_kernel` = `github.com/clio-vega/tworow-d4-kernel`,
commit `6984fff`, Lean 4.30.0 / Mathlib v4.30.0.
**File.** `TworowD4Kernel/MConvexExchange.lean` (new, 340 lines).
**Paper proof.** `proofs/2026-09-20-c1-cylindric-M-convexity.tex`, `lem:tight` and
`prop:perm-mconvex` — the fourth and last pillar of `thm:main`, and the one the Lean
development had never seen.

## Rung reached

**Rung 1**, for the paper's own definition of M-convexity, with no bound on `ℓ` and no
bound on `λ̂` — but with `J` in the **subset form**, not the sorted form. See "what Lean
cannot see" below; the LEAN.md brief's `InSupp` used `Dom (sortDesc α) lhat`, and I did
not formalise the bridge to it.

## What builds sorry-free

Everything. `lake build` green (2990 jobs), **0 sorries** (the single `grep sorry` hit is
the word inside a prose sentence saying there is none), no `native_decide`, no local
axiom.

| declaration | what it is |
|---|---|
| `asum`, `InSupp` | `α(S) = ∑_{c∈S} α c`; `J = {α ∈ ℕ^ℓ : α([ℓ]) = \|λ̂\|, α(S) ≤ Λ_{\|S\|} ∀S ⊆ [ℓ]}` |
| `asum_ex` | `α(S)` after the surgery `α − eᵢ + e_j` |
| `psum_sub_le` | concavity of `r ↦ Λ_r`, from `λ̂` antitone |
| `rank_submodular` | **`ρ(S) = Λ_{\|S\|}` is submodular** — *derived*, not postulated |
| `tight_union` | **`lem:tight`**: the tight family `F_α` is closed under union |
| `insupp_ex_of_no_tight` | the displayed equivalence: the surgery leaves `J` only through a tight `S ∋ j`, `i ∉ S` |
| `insupp_exchange` | **`prop:perm-mconvex`** — the exchange axiom |
| `insupp_exchange_witness_is_genuine` | non-vacuity: `λ̂=(2,1)`, `α=(2,1)`, `β=(1,2)`, forced move `(2,1) ↦ (1,2)` |
| `insupp_exchange_hyp_load_bearing` | negative control: with `α = β` the hypothesis fails and the conclusion is false |

**Sorries: none.** Nothing is assumed. No import of the paper's conclusion — in
particular the proof never says "because it is a permutahedron"; submodularity comes out
of `λ̂` being antitone.

## `#print axioms`

```
'MConvexExchange.rank_submodular'                  [propext, Classical.choice, Quot.sound]
'MConvexExchange.tight_union'                      [propext, Classical.choice, Quot.sound]
'MConvexExchange.insupp_ex_of_no_tight'            [propext, Classical.choice, Quot.sound]
'MConvexExchange.insupp_exchange'                  [propext, Classical.choice, Quot.sound]
'MConvexExchange.insupp_exchange_witness'          [propext, Classical.choice, Quot.sound]
'MConvexExchange.insupp_exchange_witness_is_genuine' [propext, Classical.choice, Quot.sound]
'MConvexExchange.insupp_exchange_hyp_load_bearing' [propext, Classical.choice, Quot.sound]
```

The standard three. `Classical.choice` is genuinely used and I am stating it rather than
passing over it: it enters at `choose S … using hchoice`, picking a tight witness set
`S_j` for each `j ∈ D`. That is the paper's own "pick `S_j ∈ F_α`". `RobinHood.robin_hood_step`
managed choice-free; this one does not, and I did not try to remove it.

## The defect in the paper

**The last step of `prop:perm-mconvex` was wrong, and Lean is what caught it.** The paper
displays

> `β(S*) − α(S*) = ∑_{j∈S*}(β_j−α_j) ≥ ∑_{j∈D}(β_j−α_j) > 0`,
> "the first inequality because every term with `j ∈ S*∖D` is `≤ 0`"

— and that stated reason gives `≤`, not `≥`. The terms outside `D` drag the sum *down*.
`linarith` refused it, and refused it with the right hypotheses in context, which is how I
saw it.

The conclusion survives, and the repair uses only ingredients the proof had already
established. Count on the **complement** `C = [ℓ] ∖ S*`: no index of `D` lies in `C` (since
`D ⊆ S*`), so `β_c ≤ α_c` throughout `C`; and `i ∈ C` with `β_i < α_i` strictly. So
`β(C) < α(C)`, and since the totals agree, `β(S*) > α(S*) = ρ(S*)`, contradicting `β ∈ J`.
That is what the Lean proof runs.

`.tex` corrected in the same session (`2026-09-20-c1-cylindric-M-convexity.tex`, proof of
`prop:perm-mconvex` plus a new `rem:mconvex-erratum`); recompiles clean; previous version
kept at `.tex.bak-0925c1-lean`.

## What Lean cannot see

Two things, both stated in the module docstring so a later reader meets them where they
would copy from:

1. **The sorted form.** `J` is formalised as `α(S) ≤ Λ_{|S|}`, not as
   `sort(α) ⊴ λ̂`. That these describe the same set is the **first sentence** of the
   paper's proof, and it is *not* formalised — it needs `Multiset.sort` and a
   rearrangement argument. So the green build says nothing about it.
2. **Murota's symmetric axiom.** What is proved is the one-sided exchange — which *is*
   the definition of M-convexity quoted in the paper's §"The problem" (lines 95–97) and
   used by WZZ and Brändén–Huh. The LEAN.md brief asserted in bold that "a one-sided
   version is a different, weaker statement and must not be labelled M-convexity"; that
   is overstated — the two are equivalent by Murota–Shioura. But that theorem is neither
   used nor formalised here, so the symmetric form `β + eᵢ − e_j ∈ J` for the **same** `j`
   is not a Lean fact in this development.

## Differential check against Python

`proofs/code-q254-lean/mconvex_exchange_check.py`, written from the definition and not
from the Lean file. Range: **all `λ̂` with `|λ̂| ≤ 9`, `ℓ ≤ 4`.**

```
(EQ)  sorted-form vs subset-form : 0 disagreements / 164 pairs (λ̂, ℓ)
(EX1) one-sided exchange         : 0 failures / 876317 triples (α, β, i)
(EX2) two-sided exchange         : 0 failures / 876317 triples (α, β, i)
```

(EQ) is exactly gap 1 above; (EX2) is exactly gap 2. Each unformalised half has an
instrument pointed at it. Raw output in `proofs/code-q254-lean/check-d9-ell4.out`.

## Registry

`proofs/registry/cylindric-M-convexity.json`: new child
`root/polymatroid-exchange/polymatroid-exchange-subset-form`, `trust: lean-verified`,
`lean: MConvexExchange.insupp_exchange`. The **parent `polymatroid-exchange` stays at
`proved`** — it states `J` in the sorted form, and that bridge is not formalised.
`registry_validate.py`: with the default root it reports 22 problems — 16 pre-existing
`sources index` warnings on other nodes, and 6 spurious `file not found under
/home/clio/projects/proofs`. Passing **`--proofs-dir .`** (registry `file:` paths are
relative to `projects/`, not to `proofs/`) leaves **0 errors and the same 16 pre-existing
warnings**, none of them on my new node. Backup at `.json.bak-0925c1-lean`.

## Tooling note — the brief's `--files-dir` was wrong again

LEAN.md said in bold: `--files-dir proofs`, "verified working 00:29 today". For **this**
registry that is wrong and produces a screen of spurious `file not found`:

```
--files-dir proofs → ... 'proofs/2026-09-20-c1-cylindric-M-convexity.tex' not found under proofs
--files-dir .      → OK: proofs/registry/cylindric-M-convexity.json is valid
```

`proofs` is correct for exactly one registry (the ribbon one the brief was verified
against); `.` is correct here, and `registry_validate.py --proofs-dir` has the identical
root-dir behaviour. This is the third morning in a row a stale `--files-dir` has been
copied forward from an authoritative artifact. **Run both, read the output, don't trust
the brief.**
