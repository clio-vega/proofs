# The leading π-coefficient and the odd-|J*| dichotomy for the d=4 fiber

**Clio — prove session 2026-06-11.**
Target (PROVE.md): the `e₂ mod 2` layer / `g≥1` step toward **even-|J\*|** for the d=4 fiber
`G_λ(i)=0 ⟺ λ=(2,2)`.

**Provenance / honesty.** The leading-coefficient lemma `C ≡ |J*| (mod π)` and the dichotomy
`|J*| odd ⟹ G_λ(i)≠0` were **already obtained on 2026-06-09**
(`proofs/2026-06-09-d4-fiber-reformulations.md`: "leading (1+i)-coeff `w=Σ_{j∈J*}o_j i^{(3μ−j)/2}
≡|J*| mod π`"). I reproduce it here (§1) with a fully self-contained 6-line proof because the
later cycles (06-10, 06-11) and the current PROVE.md drifted *back* to a different, and as it turns
out **incorrect**, framing. The genuinely new contributions of this session are §2 and §3.

This session does three things:

1. **Re-proves cleanly** (Theorem A, §1; *not new* — due to 2026-06-09): the leading π-adic
   coefficient of `G_λ(i)` is `≡ |J*| (mod π)`, giving **`|J*|` odd ⟹ `G_λ(i) ≠ 0`**, pruning
   ~72% of shapes. Included for self-containedness and because it is the correct backbone.
2. **Corrects the PROVE.md target** (Finding B, §2; *new*): the literal statement "`g≥1` on
   every tie (= `χ^λ(2^m)` even)" is **false** — the staircases `(3,2,1)` and `(4,3,2,1)` are
   counterexamples. The correct hypothesis is `|J*|≥2`. The 06-10/06-11 cycles silently conflated
   "`χ`-even" with "genuine tie".
3. **Reduces even-`|J*|` to an honest 2-adic Newton-polygon statement** (Reframing C, §3; *new*):
   even-`|J*|` ⟺ even minimizer-count of the integer half-polynomial `a(s)`/`b(s)`. Re-proves the
   same-parity fact by a transparent conjugate-pairing mechanism, and shows that mechanism
   **cannot** close even-`|J*|`. States the precise remaining gap.

Everything below is verified numerically for all `λ ⊢ 2m`, `m ≤ 8` (partitions of 16); the proof of
Theorem A is general.

---

## 0. Setup and notation

Fix `m ≥ 1` and `λ ⊢ 2m`. Work in the ring of symmetric functions with the Hall inner product
`⟨·,·⟩`; `p_k, h_k, e_k, s_λ` are power-sum, complete, elementary, Schur functions. Let `π = 1+i`,
a prime of `ℤ[i]` with `π² = 2i` (so `2 = -i π²`), and `v_π` the associated valuation
(`v_π(2) = 2`, `v_π(x) = 2v₂(x)` for `x ∈ ℤ`). The fiber value is

```
G_λ(i) = ⟨s_λ, (p₂ + π e₂)^m⟩  ∈ ℤ[i].
```

Using `p₁² = p₂ + 2e₂` and `i π = -(1-i) =: -π̄`, one has the equivalent expansion
(PROVE.md, "engine", form 3):

```
G_λ(i) = Σ_{j=0}^m  C(m,j) i^j π^j M_j,      M_j := ⟨s_λ, p₁^{2(m-j)} e₂^j⟩  ∈ ℤ_{≥0},
```

with `M₀ = f^λ` (the dimension) and `M_m = ⟨s_λ, e₂^m⟩ = K_{λ',(2^m)}`. Positivity `M_j ≥ 0`
holds because `p₁^{2(m-j)} e₂^j = h₁^{2(m-j)} e₂^j` is a product of Schur-positive functions.

Because `i^j` is a unit, the π-adic valuation of the `j`-th term is

```
val(j) := v_π( C(m,j) i^j π^j M_j ) = j + 2 v₂( C(m,j) M_j ),
```

(defined for `j` with `M_j ≠ 0`). Put

```
V := min_j val(j),     J* := { j : M_j ≠ 0, val(j) = V }.
```

`J*` is the set of lower vertices/lattice points of the π-adic Newton polygon of `G_λ(i)` on the
supporting line of slope `-1`.

---

## 1. Theorem A — the leading π-coefficient

> **Theorem A.** Write, for each `j ∈ J*`, `C(m,j) M_j = 2^{a_j} u_j` with `u_j` odd and
> `a_j = v₂(C(m,j) M_j)` (so `val(j) = j + 2a_j = V`). Then
> ```
> G_λ(i) = π^V · C  +  (terms of v_π > V),    C := Σ_{j∈J*} u_j · i^{\,j - a_j} ∈ ℤ[i],
> ```
> and
> ```
> C ≡ |J*|   (mod π).
> ```
> Consequently `v_π(G_λ(i)) = V` **iff** `|J*|` is odd, and in that case `G_λ(i) ≠ 0`.

**Proof.** Fix `j ∈ J*`. Since `2 = -i π²` we have `2^{a_j} = (-i)^{a_j} π^{2a_j}`. Therefore the
`j`-th term of `G_λ(i)` equals, *exactly*,

```
C(m,j) M_j · i^j π^j = 2^{a_j} u_j · i^j π^j
                     = u_j · i^j (-i)^{a_j} π^{2a_j + j}
                     = u_j · i^{\,j - a_j} · π^V,
```

using `-i = i^{-1}` and `2a_j + j = val(j) = V`. For `j ∉ J*` (with `M_j ≠ 0`) the same computation
gives a term of valuation `val(j) > V`, i.e. divisible by `π^{V+1}`. Summing,

```
G_λ(i) = π^V Σ_{j∈J*} u_j i^{\,j-a_j}  +  π^{V+1}(·)  =  π^V C  +  (v_π > V).
```

Now reduce `C` mod `π`. In `ℤ[i]/(π) ≅ 𝔽₂`: every odd integer `u_j ≡ 1`, and `i ≡ 1` because
`i - 1 = iπ` is divisible by `π`. Hence each summand `u_j i^{\,j-a_j} ≡ 1`, and

```
C ≡ Σ_{j∈J*} 1 = |J*|   (mod π).
```

Finally `v_π(G_λ(i)) = V ⟺ π ∤ C ⟺ C ≢ 0 (mod π) ⟺ |J*| ≢ 0 (mod π) ⟺ |J*| odd`
(for an integer `n`, `π | n ⟺ 2 | n`). If `|J*|` is odd then `v_π(G_λ(i)) = V < ∞`, so
`G_λ(i) ≠ 0`. ∎

> **Corollary A1 (pruning).** If `G_λ(i) = 0` then `|J*|` is even. Equivalently, every shape with
> `|J*|` odd is a non-vanisher. For `m ≤ 7`, `212` of `294` shapes have `|J*|` odd:
> **72.1%** (67% for `m ≤ 8`) are settled by this one line. `(2,2)`, the unique vanisher, has
> `|J*| = 2`, even.

**Numerical check.** For all `λ ⊢ 2m`, `m ≤ 8`: the predicted equivalence
"`C` cancels mod π ⟺ `|J*|` even" holds with **0 violations**, and "`|J*|` odd ⟹ `v_π(G)=V`"
holds with **0 violations** (`code/job_jstar.py`, `code/job_dichotomy.py`). Observed
`|J*| ∈ {1,2,4}` throughout.

---

## 2. Finding B — the literal `g≥1` target is false; the right hypothesis is `|J*|≥2`

The PROVE.md target is phrased via the **exact lift**
`Φ_λ(z) = ⟨s_λ,(p₁²+z p₂)^m⟩`. With `w = 1+z` and `p₁²+z p₂ = (1+z)p₂ + 2e₂`,

```
Φ̂(w) = ⟨s_λ, (w p₂ + 2e₂)^m⟩ = Σ_{r=0}^m C(m,r) 2^r R_r w^{m-r},   R_r := ⟨s_λ, p₂^{m-r} e₂^r⟩,
```

and `g` is the order of `(1+z)=w` in `Φ̂/2^e mod 2`, `e = v₂(content Φ)`. The claim was
`g ≥ 1` (`(1+z) | Φ̂/2^e`) for every "tie", with "tie" read as `χ^λ(2^m)=R₀` even.

> **Finding B.** This is false. For the **staircases**
> ```
> λ = (3,2,1),  m=3:   R = (0,0,0,2),   Φ̂(w) = 2³·2 = 16,   so Φ_λ(z) ≡ 16  is constant, g = 0.
> λ = (4,3,2,1), m=5:  R = (0,…,0,24),  Φ̂(w) = 2⁵·24,        Φ_λ(z) constant, g = 0.
> ```
> Both have `χ^λ(2^m) = R₀ = 0` (even), so they *are* in the `χ`-even class, yet `g = 0`.

The reason is structural: for these self-conjugate staircases **only `R_m ≠ 0`**, so
`G_λ(i) = R_m π^m` is a *single* π-term — i.e. `|J*| = 1`. They are not genuine ties at all.

**Conclusion.** The correct hypothesis throughout this circle of ideas is **`|J*| ≥ 2`**, not
"`χ^λ(2^m)` even". (The memory entry's "`g≥2` on 1624/1624 ties" is consistent only under the
reading "tie = `|J*|≥2`"; the `χ`-even reading is refuted by the two staircases above.)

This also pins *why* `g` is not the right invariant for even-`|J*|`: `Φ` lives in the **coarse**
(untilted) `w`-polygon, while `J*` lives in the **sharp** `M_j`-polygon; the two are related by the
2-adic Möbius tilt `Φ(z) = (1+z)^m P(-2z/(1+z))`, where `P(x) := Σ_j C(m,j) M_j x^j` and
`G_λ(i) = P(-π̄)`. Coarse `g≥1` neither implies nor is implied by `|J*|` even.

---

## 3. Reframing C — even-`|J*|` as an honest 2-adic Newton-polygon statement

`Theorem A` reduces the *vanishing* problem to even-`|J*|` shapes but does **not** prove
`|J*|` is even (the congruence `C ≡ |J*|` is symmetric in the parity of `|J*|`). Here is the
cleanest reformulation I obtained, together with what it does and does not settle.

### 3.1 Same parity of `J*` (H1), via conjugate pairing

`val(j) = j + 2v₂(…) ≡ j (mod 2)`, so all `j ∈ J*` share the parity of `V`; gaps in `J*` are even.

A transparent *mechanism*: `P(x) = Σ_j C(m,j) M_j x^j ∈ ℤ[x]` has **real** coefficients, and `J*`
is the slope-`(-1)` edge of its π-adic Newton polygon, whose horizontal length `j_max − j_min`
equals the number of roots `ρ` of `P` with `v_π(ρ) = 1`. A real (rational) number has *even*
π-valuation, so no real root has `v_π = 1`; hence the `v_π=1` roots occur in genuine
complex-conjugate pairs and their count `j_max − j_min` is **even**. (This recovers, conceptually,
the same-parity fact.)

### 3.2 The reduction to half-polynomials

Split `G_λ(i) = A + (-π̄) B` by parity of `j`, using `(-π̄)² = (1-i)² = -2i`:

```
A = Σ_t C(m,2t)   M_{2t}   (-2i)^t = a(-2i),   a(s) := Σ_t C(m,2t)   M_{2t}   s^t ∈ ℤ[s],
B = Σ_t C(m,2t+1) M_{2t+1} (-2i)^t = b(-2i),   b(s) := Σ_t C(m,2t+1) M_{2t+1} s^t ∈ ℤ[s].
```

If `V` is **even**, the minimum valuation is attained inside `A` (odd-`j` terms lie strictly above),
and since `v_π(α_t (-2i)^t) = 2(v₂(α_t) + t)` with `α_t = C(m,2t)M_{2t}`,

```
|J*| = #{ t : v₂(α_t) + t = min_{t'} (v₂(α_{t'}) + t') }
     = (number of on-line lattice points of the slope-(-1) edge of the HONEST 2-adic
        Newton polygon of a(s) over ℤ₂).
```

If `V` is **odd**, the identical statement holds with `b(s)` in place of `a(s)`.

> **Reframing C.** `|J*|` even (when `≥2`) is equivalent to: *the 2-adic Newton-polygon minimizer
> set of the integer half-polynomial `a(s)` (V even) or `b(s)` (V odd) has even cardinality.*

**Numerical check.** This count equals `|J*|` with **0 mismatches** for all `λ ⊢ 2m`, `m ≤ 7`
(`code/job_reframe.py`). Moreover the minimizer set is always a **2-adic box**:
in the `j`-coordinate `J* − min(J*) ∈ { {0,2}, {0,4}, {0,2,4,6} }`, i.e.
`{ Σ_{a∈S} 2^a : S ⊆ {a₁<…<a_k} }` with all `a_i ≥ 1`, so `|J*| = 2^k`. Never `{0,2,4}`.

### 3.3 Why this is genuinely harder than conjugate-pairing

Over `ℤ₂` there is **no** complex conjugation: roots of `a(s)` with `v₂ = 1` need not pair up
(a rational root *can* have `v₂ = 1`, e.g. `s = 2`). So the clean pairing argument of §3.1 does
**not** transfer to the half-polynomial, and the even/box structure of the `a(s)`-minimizers is a
real arithmetic fact about the Kostka-type numbers `M_{2t} = ⟨s_λ, p₁^{2(m-2t)} e₂^{2t}⟩`, not a
soft consequence of reality. This is the same "unbounded-depth 2-adic" wall recorded previously,
now stated in its sharpest, ℤ-only form.

---

## 4. Status

| Statement | Status |
|---|---|
| `C = Σ_{j∈J*} u_j i^{j-a_j}`, `C ≡ |J*| (mod π)` | **Proved** (Thm A) |
| `|J*|` odd ⟹ `v_π(G)=V` ⟹ `G_λ(i) ≠ 0` | **Proved** (Thm A); 0 c.ex. m≤8 |
| `G_λ(i)=0 ⟹ |J*|` even (prunes 72% of shapes) | **Proved** (Cor A1) |
| literal "`g≥1` on `χ`-even ties" | **False** (Finding B: staircases) |
| correct hypothesis is `|J*|≥2` | established |
| same parity of `J*` (H1) | **Proved** (conjugate pairing, §3.1) |
| even-`|J*|` ⟺ even 2-adic-minimizer count of `a(s)`/`b(s)` | **Proved equivalence** (§3.2); verified m≤7 |
| `|J*|` even (≥2) / box `|J*|=2^k` | **OPEN** — the remaining wall (§3.3) |

## 5. Honest gap (for the dream cycle / Robin)

The one missing brick is: **the slope-`(-1)` edge of the 2-adic Newton polygon of
`a(s)=Σ_t C(m,2t)M_{2t}s^t` (resp. `b(s)`) carries an even number of on-line lattice points,
in fact a 2-adic box `{Σ_{a∈S}2^a}`.** Theorem A shows this is exactly the obstruction to a
shape being a non-vanisher by the "odd-leading-coefficient" route; everything else in the
even-`|J*|` story is now either proved or a verified equivalence.

A pairing/involution argument cannot come from reality (no `ℚ₂` conjugation, §3.3); it must use the
arithmetic of the `M_{2t}`.

**A probe that narrows the search (this session).** For the full-box shapes `J* = {0,2,4,6}` (e.g.
`(6,3,1,1,1)`, `(4,4,2,2)` at `m=6`) the half-poly coefficients `α_t = C(m,2t)M_{2t}` have
`v₂ = (4,3,2,1)` — *exactly linear*, so the slope-`(-1)` edge is a single contiguous run and
`|J*|=4`. The edge-polynomial is then an **irreducible cubic over ℤ** with all three roots at
`v₂=1` (e.g. `5s³+870s²+5460s+1848`, content `2`). So the box/even structure is **not** a rational
factorization of `a(s)`: the "factor through `(1+2^{a}s)`" handle is dead (such a factor would give
a root of *negative* `v₂`, whereas every edge root has `v₂=+1`). Note also the two flavours of
`|J*|=2`: `{0,2}` (edge length 1) vs `{0,4}` (edge length 2 but interior `t=1` lifted off the line),
so `|J*|` is genuinely *not* a function of the edge length — the even-ness lives in *which* interior
lattice points return to the supporting line. The remaining viable handle is the SYT model
`M_{2t} = Σ_μ (#vertical-2-strip chains λ→μ of length 2t) f^μ` read 2-adically, i.e. an honest
combinatorial account of `v₂(α_t)` good enough to forbid an odd number of on-line returns.

## Verification

- `code/job_jstar.py` — `|J*|`, `val`, leading `C`, `v_π(C)`; `|J*|∈{1,2,4}`, no odd ≥3,
  0 violations of `C cancels ⟺ |J*| even`, `m ≤ 8`.
- `code/job_dichotomy.py` — exact `G_λ(i)`; `|J*|` odd ⟹ `v_π(G)=V` (0 c.ex.); 72.1% pruned (m≤7).
- `code/job_gtest.py` — coarse `g` via `Φ̂`; finds the two staircase `g=0` counterexamples.
- `code/job_reframe.py` — half-polynomial 2-adic count `= |J*|` (0 mismatches m≤7); box patterns.
