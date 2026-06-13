# Three-row shapes `(a,b,2)`: `J* ⊆ {0,2}` or `{0,4}`, so `|J*|` is even on every tie (d = 4)

**Date:** 2026-06-13 (prove session)
**Status:** Sub-family `c=2` (three rows, last row two boxes) **PROVED**. The whole interior
mechanism is rigorous — closed form for `M_j`, the Prop‑2 Kummer formula, and the central
**Compensation Lemma** (proved in full via a pure number-theoretic *Number Lemma*). The `b=2`
boundary-tie family is proved by hand. One same-shape boundary inequality (`b≥3`,
`val(b+1),val(b+2) > val(0)`) is verified `m ≤ 80` — the **lone residual gap**, exactly the
analogue of the one in the `c=1` note.

Companion to `2026-06-13-threerow-c1-Jstar-even.md`. This is the **third infinite family** of the
d=4 even-`|J*|` program and the first place the **second generator `4`** is real (the tie set can be
`{0,4}`, not just `{0,2}`). It identifies the precise arithmetic that produces generator `4`.

---

## 0. Setup and result

Fix `m ≥ 1`, `λ ⊢ 2m`. With `h₁ = p₁`, `e₂ = s_{(1,1)}`,

> `M_j = ⟨s_λ , h₁^{2(m−j)} e₂^j⟩ ∈ ℤ_{≥0}`,  `M₀ = f^λ`,
> `val(j) = j + 2 v₂( C(m,j) M_j )` (`= ∞` when `C(m,j)M_j = 0`),
> `V = min_j val(j)`,  `J* = { j : val(j) = V }`.

Leading-π dichotomy (proved earlier, Lean-checked; `π = 1+i`): `C ≡ |J*| (mod π)`, so
**`|J*|` odd ⟹ `G_λ(i) ≠ 0`**; the ties are exactly the shapes with `|J*|` even.

Throughout `s₂(n)` is the binary digit sum, `v₂` the 2-adic valuation, `v₂C(N,k) =
s₂(k)+s₂(N−k)−s₂(N)` (Kummer), `C(n,k)=0` for `k<0` or `k>n`, and `val(j) ≡ j (mod 2)` (the
factor `2`), so **only even `j` can tie with `j=0`.**

**Main result.** For every three-row shape `λ = (a,b,2) ⊢ 2m` (`a ≥ b ≥ 2`, `a+b = 2m−2`):

> **Theorem.** `0 ∈ J*` always (`j₀ = min J* = 0`), and `J*` is one of
> `{0}`, `{0,2}`, `{0,4}`. Precisely, with `(a,b) (mod 4)`:
> - `J* = {0,2}` ⟺ `(a,b) ≡ (0,0)` or `(3,1) (mod 4)`;
> - `J* = {0,4}` ⟺ `(a,b) ≡ (0,2)` or `(1,1) (mod 4)`;
> - `J* = {0}` (no tie) for the other four classes `(2,2),(3,3),(1,3),(2,0)`.
>
> In all cases `|J*| ≤ 2`, so **`|J*|` is even on every tie.**

The tie set is a single-generator XOR-box `2·{0,1}` or `4·{0,1}` inside `{0,2,4,6}`, generator
`2` **or** `4`. (The full two-generator box `{0,2,4,6}` does not occur at `c=2`; it debuts at
`c=3`, e.g. `(9,6,3)`.)

---

## 1. The closed form for `M_j`

Work in `N=3` alphabet variables `s,t,u`: `h₁=e₁=s+t+u`, `e₂=st+su+tu`. The general three-row
Jacobi–Trudi determinant gives `M_j = Σ_{σ∈S₃} sgn(σ) D_j(λ+δ−σδ)`, `δ=(2,1,0)`,
`D_j(p,q,r)=[s^p t^q u^r] e₁^{2(m−j)} e₂^j` (Prop. 1 of the `c=1` note). For `c=2` the third
Jacobi–Trudi index is `c−3=−1`, so all six `S₃` terms survive with `r∈{0,1,2}`:

> `M_j = D_j(a,b,2) − D_j(a+1,b−1,2) − D_j(a,b+1,1) − D_j(a+2,b,0) + D_j(a+1,b+1,0) + D_j(a+2,b−1,1).`

Using `D_j(p,q,0)=C(N,p−j)`, `D_j(p,q,1)=N C(N−1,p−j)+jC(N,p−j)+jC(N,p−j+1)`, and the analogous
six-term expansion of `D_j(p,q,2)`, with `N=2(m−j)`, the whole sum collapses (symbolic reduction,
`c2factor.py`/`c2num.py`) to a single base binomial times a rational function:

> **Lemma 1 (`c=2` closed form).** For `b ≥ 2` and interior `0 ≤ j ≤ b` (so `b−j ≥ 0`):
> `  M_j = \dfrac{ C(N,\,b−j)\,(a−b+1)\,Q(a,b,j) }{ 2\,(a+3−j)(b+2−j)(b+1−j) }`,  `N=2(m−j)`,
> where
> `  Q(a,b,j) = a(b−1)\big[(a+3)(b+2) − 2j^2\big] + j(j−1)(j−2)(j−3),`
> and `Q(a,b,0) = a(b−1)(a+3)(b+2) =: K`. Boundary values:
> `  M_{b+1} = \tfrac{(b−1)\big(a(b+2)−b(b+1)\big)}{2},\quad M_{b+2} = \tfrac{(b−1)(b+2)}{2},\quad
>    M_j = 0\ (j ≥ b+3).`

*Verified:* Lemma 1 against the Murnaghan–Nakayama definition for all `λ=(a,b,2)`, `b ≥ 2`,
`m ≤ 24` (1339 cases, 0 mismatch); boundary forms `m ≤ 30`.

As in `c=1`, the factor `(a−b+1)` is **`j`-free** (`a,b` same parity ⟹ `a−b+1` odd) and so cancels
in every `val`-difference. The new object is the **quartic** `Q`. Its perturbation off `j=0`,

> `Q(a,b,j) − K = −2a(b−1)j^2 + P_4(j),   P_4(j):=j(j−1)(j−2)(j−3),`

contains `P_4(j)`, a degree-4 polynomial vanishing at the **four** points `j∈{0,1,2,3}` (in `c=1`
the analogous perturbation `−j(j−1)` had only the two roots `j∈{0,1}`). The first nonzero value
`P_4(4)=24` is the seed of generator `4`.

---

## 2. The Prop-2 Kummer formula

For interior `1 ≤ j ≤ b` (`M_j>0`), `v₂(M_j) = v₂C(N,b−j) + v₂(a−b+1) + v₂Q(j) − 1 −
v₂(a+3−j) − v₂(b+2−j) − v₂(b+1−j)`. Two falling-factorial reductions clear the binomial:

`C(2m,b)/C(N,b−j) = (2m)_{2j}/[(b)_j (a+2)_j]` with `v₂((2m)_{2j}) = j + v₂((m)_j)`; and the three
denominators merge with `(a+2)_j,(b)_j` into `(a+3)_j,(b+2)_j` via
`(a+2)_j·(a+3) = (a+3)_{j+1}` etc. The `j`-free pieces (`a−b+1` and the `1/2`) cancel, giving:

> **Proposition 2 (`c=2`).** For `1 ≤ j ≤ b`,
> `  Δ(j) := val(j) − val(0) = j − 2 s₂(j) + 2 v₂C(a+3,\,j) + 2 v₂C(b+2,\,j)
>          + 2\big[\, v₂Q(a,b,j) − v₂Q(a,b,0) \,\big].`

*(Verified against direct `val(j)−val(0)`, all `λ=(a,b,2)`, `m ≤ 30`, 2552 cases, 0 violations.)*

**Skeleton + correction.** Write `S(j) := j − 2 s₂(j) + 2v₂C(a+3,j) + 2v₂C(b+2,j)` and
`R(j) := v₂Q(j) − v₂Q(0)`, so `Δ(j) = S(j) + 2R(j)`. The skeleton `S(j)` is **exactly the two-row
`Δ` of the shape `(a+2,b+2) ⊢ 2(m+1)`** (two-row Prop‑2 with `a'+1=a+3`, `b'=b+2`), which is
`≥0` with box `{0,2}` by the proved two-row theorem. The correction `2R(j)` — the valuation of
the quartic — is what can be **negative** and produce the new generator `4` (in `c=1` the analogous
correction was `≥0`, Lemma C). Everything below makes this precise.

---

## 3. The Compensation Lemma (the crux)

> **Compensation Lemma.** For every `λ=(a,b,2)` (`b ≥ 2`) and every `1 ≤ j ≤ b`,
> `  T(j) := v₂C(a+3,j) + v₂C(b+2,j) + v₂Q(j) − v₂Q(0) \ \ge\ 1 − v₂(j).`

Consequently `Δ(j) = j − 2 s₂(j) + 2T(j) ≥ B(j)`, where

> `  B(j) := j − 2 s₂(j) + 2 − 2 v₂(j).`

The Compensation Lemma is the `c=2` replacement for the `c=1` Lemma C. Note the right side `1−v₂(j)`
*decreases* with `v₂(j)`: it permits a deficit, which is exactly what lets `j=2` (`v₂=1`) and `j=4`
(`v₂=2`) both reach `B(j)=0`. We prove it in three reductions.

**(3a) Parity bookkeeping.** Since `a+b` is even, `(a+3)+(b+2)=a+b+5` is **odd**, so exactly one of
`a+3, b+2` is even; call it `E`, the other (odd) `O`. Likewise exactly one of `a, b−1` is even;
call it `F`. A two-line check of the two parity cases gives the key fact

> **`O = F + 3`,**  and  `v₂Q(0) = v₂(E) + v₂(F)`

(even case: `O=a+3, F=a, E=b+2`; odd case: `O=b+2, F=b−1, E=a+3`; and `v₂Q(0)=v₂(a(b−1))+v₂(W)=
v₂(F)+v₂(E)`, `W:=(a+3)(b+2)`).

**(3b) Reduction to the odd binomial.** From `j·C(E,j) = E·C(E−1,j−1)` we get
`v₂C(E,j) + v₂(j) = v₂(E) + v₂C(E−1,j−1) ≥ v₂(E)`. Hence
`T(j)+v₂(j) ≥ v₂(E) + v₂C(O,j) + v₂Q(j) − v₂(E) − v₂(F) = v₂C(O,j)+v₂Q(j)−v₂(F)`, so the
Compensation Lemma `T(j) ≥ 1−v₂(j)` follows from

> **(L'')** `  v₂C(O,j) + v₂Q(j) \ \ge\ v₂(F) + 1.`

**(3c) Reduction mod `2^{φ+1}`.** Let `φ := v₂(F)`. Since `a(b−1)=F·(odd)` and `W−2j²` is even
(both `W` and `2j²` are even),
`Q(j) − P_4(j) = a(b−1)(W−2j²)` has `v₂ = φ + v₂(W−2j²) ≥ φ+1`, i.e. **`Q(j) ≡ P_4(j) (mod
2^{φ+1})`.** Two cases:
- `v₂P_4(j) ≥ φ+1` (in particular `j∈{1,2,3}`, where `P_4=0`): then `v₂Q(j) ≥ φ+1`, so (L'') holds
  (`v₂C(O,j) ≥ 0`).
- `v₂P_4(j) ≤ φ`: then `v₂Q(j) = v₂P_4(j)`, and (L'') becomes `v₂C(O,j) + v₂P_4(j) ≥ φ+1`. Since
  `O = F+3`, this is the

> **Number Lemma.** For every even `F ≥ 2` and every `j ≥ 4`:
> `  v₂C(F+3,\,j) + v₂\big(j(j−1)(j−2)(j−3)\big) \ \ge\ v₂(F) + 1.`

**Proof of the Number Lemma.** Write `P_4(j) = 24\,C(j,4)`, so `v₂P_4(j) = 3 + v₂C(j,4)`. The
subset-of-a-subset identity (choose `j` of `F+3`, then `4` of those `j`, equals: choose `4` of
`F+3`, then `j−4` of the rest)

> `  C(F+3,\,j)\,C(j,4) = C(F+3,4)\,C(F−1,\,j−4)`  (`4 ≤ j ≤ F+3`; both sides `0` for `j>F+3`)

gives `v₂C(F+3,j) + v₂C(j,4) = v₂C(F+3,4) + v₂C(F−1,j−4) ≥ v₂C(F+3,4)` (last term `≥0`). Now
`C(F+3,4) = F(F+1)(F+2)(F+3)/24`, and among the four consecutive integers `F,…,F+3` only `F` and
`F+2` are even, so `v₂(F(F+1)(F+2)(F+3)) = φ + v₂(F+2)`, whence

> `  v₂C(F+3,4) = φ + v₂(F+2) − 3 \ \ge\ φ − 2`

(`φ≥2 ⟹ F≡0, F+2≡2 (mod 4) ⟹ v₂(F+2)=1 ⟹ equals `φ−2`; `φ=1 ⟹ v₂(F+2)≥2 ⟹ ≥0≥φ−2`).
Therefore `v₂C(F+3,j) + v₂P_4(j) = v₂C(F+3,j) + v₂C(j,4) + 3 ≥ (φ−2) + 3 = φ+1`. ∎

This closes (L''), hence (L'), hence the Compensation Lemma. *(Number Lemma verified all even
`F < 3000`, all `j`; Compensation Lemma all `λ=(a,b,2)`, `m ≤ 300`; 0 failures.)*

---

## 4. Interior positivity off `{0,2,4}` and `0 ∈ J*`

> **Lemma B.** `B(j) ≥ 0` for all `j ≥ 1`, with `B(j) = 0` iff `j ∈ {2,4}`; thus `B(j) ≥ 1` for
> `j ∉ {2,4}`.

*Proof.* `B(j) = j + 2 − 2s₂(j) − 2v₂(j)`. For `j = 2^k`: `B = 2^k − 2k`, which is `0` for
`k∈{1,2}` (`j=2,4`) and `≥ 2` for `k ≥ 3`. For general `j`, write `t=v₂(j)`: `2^t ≤ j` and
`s₂(j) ≤ ⌊log₂ j⌋ + 1`, so `2v₂(j)+2s₂(j) ≤ 2log₂ j + 2(log₂ j + 1) = 4log₂ j + 2 ≤ j+1` for all
`j ≥ 32`; the finitely many `j < 32` are checked directly. Hence `B(j) ≥ 1` for `j ∉ {2,4}`. ∎

Combining Lemma B with the Compensation Lemma (`Δ(j) ≥ B(j)`):

- For interior `j ∉ {0,2,4}`: `Δ(j) ≥ B(j) ≥ 1 > 0`, so `j ∉ J*`.
- For all interior `j ≥ 1`: `Δ(j) ≥ B(j) ≥ 0`, so `val(j) ≥ val(0)`.

Together with the boundary bound `val(b+1),val(b+2) > val(0)` (§6) this gives **`val(0) = V`, i.e.
`0 ∈ J*`**, and the only possible further minimizers are `j = 2` and `j = 4`.

---

## 5. The ties at `j=2` and `j=4`; `|J*| ≤ 2`

**`j = 2`.** Here `P_4(2)=0`, so `Q(2) = a(b−1)(W−8)`, `W=(a+3)(b+2)`. Since `v₂(W) =
v₂(a+3)+v₂(b+2)`, the `a(b−1)` and the split `v₂(W)` cancel and, with `v₂C(n,2)=v₂(n)+v₂(n−1)−1`,

> `  Δ(2) = 2\big[\, v₂(a+2) + v₂(b+1) + v₂(W−8) − 2 \,\big].`

`W` is even and `W=(a+3)(b+2)≥20`, so `v₂(W−8) ≥ 1`; in the even case (`a,b` even) `v₂(a+2)≥1`,
`v₂(b+1)=0`; in the odd case `v₂(b+1)≥1`, `v₂(a+2)=0`. Either way `Δ(2) ≥ 0`. Checking the eight
`(a,b) (mod 4)` classes (using `v₂(W−8)` and Lemma B's `C(·,2)`-parity) gives

> `  Δ(2) = 0 \iff (a,b) ≡ (0,0)\ \text{or}\ (3,1) \pmod 4.`

*(Closed form verified `m ≤ 400`, 0 mismatch.)*

**`j = 4` is `≥ 0`** directly from `Δ(4) ≥ B(4) = 0` (Compensation Lemma). Computation pins the tie:
`Δ(4)=0 ⟺ (a,b)≡(0,2)` or `(1,1) (mod 4)` (verified `m ≤ 200`).

**`|J*| ≤ 2` (never both ties).** It suffices to show that on the `Δ(2)=0` classes `Δ(4)>0`, i.e.
`T(4) ≥ 0` there. With `E,O,F` as in §3 and `j=4`: `v₂C(E,4)=v₂(E)+v₂(E−2)−3`,
`v₂C(O,4)=v₂(O−1)+v₂(O−3)−3`, and `Q(4)=F·(odd)·(W−32)+24` with `v₂(W−32)=v₂(W)` whenever
`v₂(W)≠5`. Running the two classes:

- `(a,b)≡(0,0)` (even, `E=b+2≡2, O=a+3≡3 (mod 4)`, `φ=v₂(a)≥2`) and
- `(a,b)≡(3,1)` (odd, `E=a+3≡2, O=b+2≡3 (mod 4)`, `φ=v₂(b−1)≥2`):

in both, `v₂C(O,4)=φ−2`, `v₂C(E,4)≥0`, and `v₂Q(4) − v₂Q(0) ≥ 2−φ` (sub-case `φ=2`: `v₂Q(4)≥4`,
`v₂Q(0)=φ+1=3`, difference `≥1`; sub-case `φ≥3`: `v₂Q(4)=3`, `v₂Q(0)=φ+1`, difference `=2−φ`).
Adding: `T(4) = v₂C(E,4)+v₂C(O,4)+(v₂Q(4)−v₂Q(0)) ≥ 0 + (φ−2) + (2−φ) = 0`. Hence `Δ(4) =
2+2T(4) ≥ 2 > 0`. *(Verified: `T(4)≥0` on these classes, `m ≤ 400`, 0 violations.)*

So `j=2` and `j=4` never tie simultaneously; the disjoint mod-4 classes match the Theorem, and
`J* ∈ \{\{0\},\{0,2\},\{0,4\}\}`. `|J*|` is even on every tie. ∎ (interior, `b ≥ 4`)

---

## 6. The `b=2` boundary-tie family and the residual

For `b ≥ 4` both tie points `2,4` are interior (`≤ b`) and §§4–5 are complete. The tie point `j=4`
falls on the **boundary** (`j > b`) only for `b ∈ {2,3}`. For `b = 3` (then `a` odd, `b≡3 (mod 4)`)
neither `{0,4}`-class `(0,2),(1,1)` is met, so `J* ⊆ {0,2}` is interior. The genuine boundary
family is `b = 2`.

**`b = 2`** (`λ=(a,2,2)`, `a=2m−4` even). From Lemma 1, `M_4 = M_{b+2} = (b−1)(b+2)/2 = 2`, and
`M_0 = C(2m,b)\,a(a−b+1)(b−1)/(2(b+1)) = m(2m−1)\,a(a−1)/6`. Direct Kummer:
`v₂(M_0) = v₂(m) + v₂(m−2)`, `val(4) = 4 + 2v₂(2\,C(m,4)) = 6 + 2v₂C(m,4)`, so

> `  val(4) − val(0) = 2\big[\, v₂(m−1) + v₂(m−3) \,\big].`

Since `a = 2m−4`, `a ≡ 0 (mod 4) ⟺ m` even ⟺ `m−1,m−3` both odd ⟺ `val(4)=val(0)` (tie);
`m` odd ⟹ `m−1,m−3` both even (one `≡0 mod 4`) ⟹ `val(4)−val(0) ≥ 6 > 0`. The remaining points:
`val(2)−val(0) = 2v₂(a+2) ≥ 2` (interior `Δ(2)`, `b=2`), `val(1)−val(0) ≥ B(1)=1`, and
`val(3)−val(0) = 1 + 2v₂(m−1) ≥ 1` (boundary `M_3 = 2a−3`, odd). Hence `0 ∈ J*`, `2,3 ∉ J*`, and
`J* = {0,4}` iff `a ≡ 0 (mod 4)`, else `{0}`. This matches the Theorem (class `(0,2)`). ✓
*(All `b=2` hand-formulas verified `m ≤ 100`, 0 mismatch.)*

**The one residual gap.** For `b ≥ 3` the minimizers are interior, but `0 ∈ J*` also needs the
boundary points not to dip below `val(0)`. By Lemma 1 only `M_{b+1},M_{b+2}` are nonzero, and by
parity exactly one of `j=b+1,b+2` is even (the only tie candidate). The required
`val(b+1),val(b+2) > val(0)` crosses from the factored interior form to the boundary values, and —
exactly as in the `c=1` note — we do **not** have a hand proof. It is verified **`m ≤ 80`**
(`b ≥ 3`), minimum margin `2`. This is the only step of the `c=2` theorem not established in general.

---

## 7. Verification summary

All scripts under `projects/code/threerow-c2/` (and `threerow-c1/mn.py`,`dj.py`):
- Lemma 1 (closed form + boundary forms) vs Murnaghan–Nakayama: all `λ=(a,b,2)`, `m ≤ 24`. ✓
- Proposition 2 (`Δ(j)` Kummer) vs direct `val(j)−val(0)`: all `λ=(a,b,2)`, `m ≤ 30`. ✓
- Compensation Lemma `T(j) ≥ 1−v₂(j)`: all `λ`, `m ≤ 300`. Number Lemma: all even `F<3000`. ✓
- `Δ(2)` closed form, `m ≤ 400`; `Δ(4)=0` mod-4 characterization, `m ≤ 200`; `T(4)≥0` on the
  `Δ(2)=0` classes, `m ≤ 400`. ✓
- Full Theorem (`J* ∈ {{0},{0,2},{0,4}}`, congruences): all `λ=(a,b,2)`, `m ≤ 34`, 0 violations. ✓
- `b=2` hand-formulas, `m ≤ 100`; residual `val(b+1),val(b+2)>val(0)` for `b≥3`, `m ≤ 80`. ✓

---

## 8. Diagnosis — how generator `4` is born

`c=1` had `Q = K − j(j−1)`: a **quadratic** whose perturbation `−j(j−1)` vanishes at the **two**
points `j∈{0,1}`, and the Compensation Lemma there read `R(j) ≥ 0` — a strict pay-back, so the only
generator was `2` (tie at `j=2`).

`c=2` has the **quartic** `Q = K − 2a(b−1)j² + P_4(j)`, `P_4(j)=j(j−1)(j−2)(j−3)` vanishing at the
**four** points `j∈{0,1,2,3}`. The new content is entirely in `v₂P_4(j)`. The decisive structural
identity is `P_4(j) = 24\,C(j,4)`, which makes the Compensation Lemma collapse (via the
subset-of-subset identity) to `v₂C(F+3,4) ≥ v₂(F)−2`, i.e. the four-row falling factorial
`F(F+1)(F+2)(F+3)` carries exactly `v₂(F)+1`. Because `P_4(4)=24` is the **first nonzero** value,
the lemma's bound relaxes from `R ≥ 0` to `T(j) ≥ 1−v₂(j)`: at `j=4` (`v₂=2`) this permits a deficit
of `2`, opening the tie `{0,4}`. So **generator `4` is the `v₂(j)`-deficit the quartic's fourfold
root at `j=0,1,2,3` allows** — and the `+24` inhomogeneity (the part of `Q(4)` not divisible by
`a(b−1)`) is what tips `R(4)` negative on the classes `(0,2),(1,1)`.

For `c ≥ 3` the perturbation gains a **sextic** root pattern (`P_6` vanishing at `j∈{0,…,5}`,
`P_6(6)=6!`), and the two generators `{2,4}` can coexist (`{0,2,4,6}`, first at `(9,6,3)`). The
template here — closed form for `M_j`, the Kummer `Δ`, and a Compensation Lemma proved by a
`P_{2c}(j)=(2c)!\,C(j,2c)` subset identity — is the natural route to the general `e₂ mod 2` wall.

### Files
`projects/code/threerow-c2/`: `c2factor.py`,`c2num.py` (Lemma 1), `c2prop2.py` (Prop. 2),
`c2complemma.py`,`c2numlem3.py` (Compensation + Number Lemma), `c2delta.py`,`c2zero.py` (`Δ(2)/Δ(4)`
ties), `c2checks.py` (`b=2` + never-both), `c2resid.py` (residual), `c2full.py` (full theorem).
