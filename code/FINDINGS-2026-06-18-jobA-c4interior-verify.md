# FINDINGS — Job A (2026-06-18 CODE): verify + stress-test the `c=4` interior even-|J*|

**Script:** `jobA_c4_interior_verify.py` · **Output log:** `results_jobA_c4_interior.txt`
**Target:** independently validate the PROVE deliverable
`proofs/2026-06-18-c4-interior-Jstar-even.md` (three-row `(a,b,4)` interior, `J*⊆{0,2}`).
**Ground truth:** Murnaghan–Nakayama `M_j=⟨s_λ,h₁^{2(m−j)}e₂^j⟩` via the abacus χ engine
(`threerow-c3/mn.py`). Every closed form is cross-checked against MN before use.

## Verdict (one line)

**CONFIRMED.** The `c=4` interior theorem holds exactly as PROVE states: `J*⊆{0,2}`, `0∈J*`
always, `|J*|≤2` even, tie `⟺ (a,b)≡(2,2),(1,3) (mod 4)`. **0 violations to `m≤120`** (closed form,
MN-validated to `m≤40`). The one named residual (heavy-free `Δ̂(j)>0`, `j≥8`) is certified `m≤120`,
**0 violations** — still a *certification*, not a uniform proof, so the §6 residual stands.

## What each section certified (all 0 violations)

| § | Claim verified | Range / gate | Result |
|---|---|---|---|
| 1 | **Closed form (Lemma 1)** `M_j=C(N,b−j)(a−b+1)Q₄/[24(a+5−j)∏₁⁴(b+i−j)]` **= MN** | full `m≤40`, both `b`-parities, all interior `j` (**10 290 cases**) + MN spot-checks `m=20,23,26` | **0 mismatch** |
| 2 | **Decomposition** `Q₄=(a−2)(b−3)H+P₈(j)`, `P₈=8!C(j,8)` | symbolic poly-division (`cancel`); `H` matches `c4finite.py`; `H(0)=(a+3)(a+4)(a+5)(b+2)(b+3)(b+4)` | **exact (=0)** |
| 3 | **Prop-2 Kummer** `Δ(j)=j−2s₂(j)+2vC(a+5,j)+2vC(b+4,j)+2[v₂Q₄(j)−v₂Q₄(0)]` vs direct `val(j)−val(0)` | `m≤23` (1 653 ties seen) | **0 mismatch** |
| 4 | **|J*| classification** (high-`m` sweep, closed form) | **`m≤120`, 6 670 shapes** | see below |
| 5 | **Five finite residue divisibilities** `2^{k_j}∣(a+2)_{j−3}(b+1)_{j−3}H(j)`, `k=3,6,6,8,8`, `j=3..7` | **exhaustive** over all residues mod `2^{k_j}`, `a≡b (mod 2)` | **0 exceptions** |
| 6 | **Residual** `Δ̂(j)>0` (heavy-free) for `j≥8` | closed-form sweep `m≤120` | **0 violations** |

## The `(a,b) mod 4` box table (§4, `m≤120`, the data Rick said was only `m≤45`)

```
 (a%4,b%4) -> J* sets observed         tie locus = {(1,3),(2,2)}
   (0,0): {(0,)}      (2,0): {(0,)}     max |J*| observed   = 2
   (0,2): {(0,)}      (2,2): {(0,2)}    0 ∈ J* failures     = 0
   (1,1): {(0,)}      (3,1): {(0,)}     J*⊄{0,2} violations = 0
   (1,3): {(0,2)}     (3,3): {(0,)}     total ties (|J*|=2) = 1 653
```
- **Minimal tie witness:** `(a,b,m)=(6,6,8)`, `J*={0,2}` (matches `(2,2) mod 4`).
- **No generator 4, 6, or 8 ever appears** — confirming PROVE's headline correction to the
  PROVE.md expectation (the tip `P₈(j)` vanishes at `j=0..7`, killing `c=2`'s gen-4 mechanism).
- The min-`Δ` floors match PROVE's table: `minΔ̂(j)=4,5,6,9,8,…` for `j=8,9,10,11,12,…` — all `>0`.

## Honesty notes
- I initially added a "`M_j=f^{(a−j,b−j,4−j)}` conjugate-dim" cross-check; **it is FALSE** for
  `j≥1` (e.g. `(18,10,4)` `j=1`: MN `35 562 372 300` ≠ `fdim(17,9,3)=3 365 041 680`). Dropped it.
  This identity is *not* part of the proof — `M_j` is the e₂-power pairing, not a hook dimension.
  The real cross-check (closed form vs MN, §1) is clean.
- The §6 residual is genuinely *certified, not proved*: for fixed `j` it is a finite divisibility
  `2^{k_j}∣(a+2)_{j−3}(b+1)_{j−3}H(j)`, but `k_j∼2·vfact(j)` grows, so no single finite check
  covers all `j≥8`. Closing it uniformly needs the **`c=4` Number Lemma** = a `v₂H(j)` floor for
  the sextic heavy quotient `H` (PROVE's stated next target). My sweep `m≤120` adds no
  counterexample.

**Net:** the `c=4` interior is proved for `1≤j≤7` and `j≥8` tip-dominated; the deep heavy-dominated
tail is certified to `m≤120` (was `m≤45`). The three-row `(a,b,4)` interior is **complete except for
the one named §6 residual**, exactly as PROVE claims — no overclaim, no hidden gap.
