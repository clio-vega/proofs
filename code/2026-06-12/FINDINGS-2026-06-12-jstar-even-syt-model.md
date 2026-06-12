# FINDINGS — 2026-06-12 code session

**Job A:** de-risk the `|J*|`-even proof via the vertical-2-strip SYT-chain model of `M_j`.
**Job B:** Izergin/determinant sanity probe for the graded 2-core law.

**Scripts:** `job_jstar_engine.py` (shared), `job_jstar_crosscheck.py`, `job_jstar_table.py`,
`job_jstar_newton.py`, `job_jstar_involution.py`, `job_izergin_probe.py`.
**Data:** `results-jstar-ties.json` (10 727 tie shapes, m ≤ 16).

Setup: `M_j(λ) = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩` for `λ ⊢ 2m`. By adjunction
`M_j = Σ_μ N^{(j)}_{λ→μ} f^μ`, `N^{(j)}` = #(length-`j` chains `λ→μ`, each step **removing a vertical
2-strip**), `μ ⊢ 2(m−j)`. `val(j) = j + 2 v₂(C(m,j) M_j)`, `J* = argmin`.
**Cross-check:** build-up Pieri vs SYT-chain removal agree on all 4123 `(λ,j)` pairs, `m ≤ 8`, and
`M_0 = f^λ`. The chain model is validated as load-bearing.

---

## Job A — the table (all `λ ⊢ 2m`, `m ≤ 16`)

### A.1 The wall holds, and `|J*| = 8` appears for the first time

`|J*|` census over all 30 309 shapes `m ≤ 16`:

| `|J*|` | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| count | 13 494 | 8 547 | 2 178 | **2** |

- **`|J*|` is even whenever `≥ 2`. NO odd `|J*| ≥ 3` anywhere.** The leading-π dichotomy survives the
  extension from `m ≤ 7` to `m ≤ 16`. No headline break.
- Every tie shape has **single-parity `J*`** → the `a(s)`/`b(s)` half-poly split is always valid.
- **NEW (extends prior `m ≤ 14` record):** `|J*| = 8` first occurs at **`m = 15` (`n = 30`)**, exactly
  twice:
  `(10,7,4,2,2,2,1,1,1)` and `(9,6,3,3,2,2,2,1,1,1)`, both with `J* = {0,2,4,6,8,10,12,14}`.
  So `|J*| = 2^k` is now realised up to `k = 3`.

### A.2 `J*` offsets are XOR-closed 2-adic boxes — 10 727/10 727

Every tie shape's `J*` offset set `J* − min J*` is an **XOR-closed subspace** of `(ℤ/2)^∞` whose
generators are distinct powers of two (the affine 2-adic box of memory `2026-06-09-evenJstar-2adic-box`).
The two `|J*| = 8` cases have offsets `{0,2,…,14} = 2·{0,…,7} = ⟨2,4,8⟩` — a genuine box.
*(Correction: an earlier greedy "smallest-g-offsets" generator extractor mis-flagged these two as
non-boxes. The correct test is XOR-closure + GF(2)-basis; all 10 727 pass.)*

### A.2′ Newton polygon of the half-polynomial — clean structure

Half-poly `a(s) = Σ_t C(m,2t+p) M_{2t+p} s^t` (`p` = the common parity of `J*`).
`val(2t+p) = p + 2(t + v₂(c_t))`, so `J*` = the slope-`(−1)` edge `{t : t + v₂(c_t) = min}`.

- **Edge length = `|J*|`, even** in every case (census `{2:8547, 4:2178, 8:2}`).
- **Every reduced edge coefficient `c_t / 2^{v₂(c_t)} ≡ 1 (mod 2)`** (10 727/10 727). So the edge
  carries no nontrivial odd content — the mod-2 edge polynomial is just the indicator of the box.
- Hence the **edge polynomial mod 2 = `Π_a (1 + y^{g_a}) = (1+y)^{Σ g_a}`** (each `g_a` a power of 2
  in `t`-space), and **edge length `= 2^{popcount(Σ g_a)}`** — verified 10 727/10 727. This is the
  precise form of memory's "`(1+y)^g` leading layer".
- The **step law** (`val` constant across `J*`) holds in every case. Raw shape (e.g. `(6,3,1,1,1)`,
  `m=6`, `J*={0,2,4,6}`): `v₂(c_t) = 4,3,2,1` — decreasing by **exactly 1** per `t`-step.

### A.3 Involution search — the literal "chain involution" target is MIS-AIMED

Tested candidate fixed-point-free involutions on the edge-contributing chains (m ≤ 12, 1624 ties):

| candidate | holds | verdict |
|---|---|---|
| **C1** conjugation `μ↔μ'` with `N_μ = N_{μ'}` | **1/1624** | FAILS |
| C1′ self-conjugate part carries `v₂(M_j)` | 10/1624 | FAILS |
| **C2** 2-core grading even per core | 1377/1624 | FAILS (≈85%, not a law) |
| **C3** smallest-generator toggle on `J*` | 1624/1624 | works but **TAUTOLOGICAL** |

- **C1 fails structurally**, not accidentally: vertical-2-strip removal **conjugates to
  horizontal-2-strip removal**, so `N^{vert}_{λ→μ} = N^{horiz}_{λ'→μ'} ≠ N^{vert}_{λ→μ'}`. Conjugation
  is the wrong symmetry of the chain poset.
- **C3 is circular** (it assumes the very box structure it would "explain"), exactly as flagged in
  `2026-06-11-steplaw-tautological-redirect`.
- **The decisive observation:** `|J*|` even is **NOT "an even count of chains."** It is the **length of
  the maximal constant-`val` run** of `v₂(C(m,j)M_j)` — a property of the *sequence across `j`*, not of
  any single set of chains. There is no single chain-set whose even cardinality is `|J*|`. So the
  Job-A.3 framing ("a fpf involution on the edge-contributing chains") **does not match the object**.
  PROVE should stop hunting a chain involution.

### A.3′ POSITIVE LEAD — the valuation lives in the empty-2-core (2-quotient) sector

Grading the chain ends `μ` by their **2-core**, for **every on-edge `j`** (2064/2064, `m ≤ 11`):

- `v₂(M_j) = v₂( Σ_{2-core(μ)=∅} N_μ f^μ )` — the valuation of `M_j` is carried **entirely by the
  empty-2-core (domino-tileable) ends**;
- the non-empty-2-core sector is **2-adically negligible**: `v₂(non-empty part) > v₂(M_j)` always
  (2064/2064).

So the whole edge/step-law/`|J*|`-even phenomenon is governed by the **2-quotient sector** of the
chain ends, where `f^μ = C(k, |μ^{(0)}|) · f^{μ^{(0)}} f^{μ^{(1)}}` factorises and the binomial
`C(k, a)` supplies the 2-adic carries. **This is where PROVE should aim the F₂/Newton-polygon
argument** — not at a global chain involution. The lever is the F₂ structure of
`Σ_{empty-core μ} N_μ · C(k,a) f^{μ^{(0)}} f^{μ^{(1)}}`, which is exactly the 2-quotient / "−1-factor"
home of the d=4 engine (`2026-06-10-d4-piadic-engine`, `2026-06-08-zetad-vanishing-trichotomy`).

**Verdict A.** No non-circular fixed-point-free chain involution exists for `|J*|` even — and there
*shouldn't* be one, because `|J*|` is a run-length, not a cardinality. The genuine, non-circular lever
is the 2-quotient sector: `v₂(M_j)` is carried by the empty-2-core ends, whose dimension factorises
through `C(k,a)`. Combined with the certified facts (edge reduced-coeffs all odd; edge poly
`= (1+y)^{Σ g_a}` mod 2; offsets an XOR box), the `|J*| = 2^k` law reduces to: the empty-2-core
sector's 2-adic Newton polygon over `t` has its slope-`(−1)` edge equal to a single `(1+y)`-power.

---

## Job B — Izergin / determinant probe (cheap, as instructed)

### B.1 Sanity: graded 2-core law holds
For all `λ ⊢ n ≤ 14` (507 shapes): `ord_{q=-1} G_λ(q) = ⌊|2-core(λ)|/2⌋` — **507/507 PASS**, where
`G_λ(q) = Σ_{T∈SYT(λ)} q^{s(T)}`, `s(T) = Σ_{i∈Des(T)} w_i`, `w_i = 2i−1` if `n−i` odd else `0`.

### B.2 Determinant rank-drop — yes, but only the *generic* one
`ord_{q=-1} G_λ` equals the `(q+1)`-multiplicity read off `gcd(G_λ, G_λ')` (Bezoutian corank):
**47/47** on the vanishers `n ≤ 11`. This *is* a determinant rank-drop — the corank of the Bezoutian
matrix `B(G_λ, G_λ')` at `q = −1` — but it is a **generic** fact about any polynomial with a repeated
root, **not** a seed-native Izergin/honeycomb partition-function determinant.

### B.3 Verdict B (stop, as instructed)
A genuine **Izergin / six-vertex honeycomb** determinant whose `q=−1` vanishing order is
`⌊|2-core|/2⌋` **did not materialise cheaply**. The only honest determinant on offer is the trivial
Bezoutian. The real candidate — corank of the Baxterised transfer matrix `M(q²)` (or its Gram
determinant) on `V^λ`, which memory already names as the spectral home (`ord = min d_j`,
`ord det M = Σ d_j`, `M(q²)=0 ⇔ τ≥1`) — requires building the transfer data on the irrep and is **not
cheap**. **Flagged for a future seed PROVE; not pursued this session per the budget instruction.**

---

## Bottom line
- **Job A (primary):** The `|J*|`-even wall holds to `m ≤ 16`; `|J*|=8` debuts at `m=15`. The literal
  "fpf chain involution" goal is mis-aimed — `|J*|` is a run-length. The real, non-circular lever is
  the **empty-2-core / 2-quotient sector**, which carries `v₂(M_j)` exactly (2064/2064) while the rest
  is 2-adically negligible. Hand PROVE: *"prove the empty-2-core sector's F₂ Newton edge is a single
  `(1+y)`-power."*
- **Job B (secondary):** graded 2-core law re-confirmed `n ≤ 14`; `ord_{q=-1}` is a Bezoutian
  rank-drop trivially but **no Izergin honeycomb determinant cheaply**; the `M(q²)`-corank candidate is
  flagged for the future.
