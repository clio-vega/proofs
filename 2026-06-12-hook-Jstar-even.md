# The leading‑π layer cancels on every hook tie: `|J*| ≤ 2` for hooks (d = 4)

**Date:** 2026‑06‑12 (prove session)
**Status:** Hook family **PROVED** (complete). General λ: reframed, sharpened, verified `m ≤ 10`, gap precisely located.

---

## 0. The problem

Fix `m ≥ 1` and `λ ⊢ 2m`. With `h₁ = p₁`, `e₂ = s_{(1,1)}`, set

> `M_j = ⟨s_λ , h₁^{2(m−j)} e₂^j⟩ ∈ ℤ_{≥0}`,  `M₀ = f^λ`,
> `val(j) = j + 2 v₂( C(m,j) M_j )` (and `val(j) = ∞` when `C(m,j)M_j = 0`),
> `V = min_j val(j)`,  `J* = { j : val(j) = V }`.

The leading‑π dichotomy (proved earlier, 06‑09/06‑11) gives, with `π = 1+i`,
`G_λ(i) = Σ_j C(m,j) M_j (i−1)^j` and `C ≡ |J*| (mod π)`, so **`|J*| odd ⟹ G_λ(i) ≠ 0`.**
The residual ties are exactly the shapes with `|J*|` even. The target:

> **Theorem (to prove, general).** For every `λ ⊢ 2m` with `|J*| ≥ 2`, `|J*|` is even.

**What this note proves.** The theorem for **every hook** `λ = (a,1^b) ⊢ 2m`; in fact the sharp
statement `|J*| ∈ {1,2}` for hooks, with the only possible tie being `J* = {0,2}`. Along the way I
record several rigorous reductions that sharpen the general program and pin down exactly where the
general case is still open.

Throughout, `s₂(n)` is the binary digit sum, `v₂` the 2‑adic valuation, and Kummer's theorem
`v₂ C(N,k) = ` (number of carries when adding `k + (N−k)` in base 2) `= s₂(k)+s₂(N−k)−s₂(N)`.

---

## 1. Reductions that hold for all λ (rigorous unless flagged)

**1.1 Single parity (H1, trivial reparametrisation).** `val(j) ≡ j (mod 2)`, so all `j ∈ J*` share
one parity `ε = V mod 2`. Writing `j = 2t+ε`, `J*` is the on‑line set of the half‑polynomial
`Σ_t C(m,2t+ε) M_{2t+ε} s^t`.

**1.2 Eisenstein evaluation.** `G_λ(i) = P_λ(α)` where `P_λ(x) = Σ_j C(m,j) M_j x^j ∈ ℤ[x]` and
`α = i−1` is a root of the **Eisenstein** polynomial `x²+2x+2` at 2. Hence `ℚ₂(α)=ℚ₂(i)` is ramified,
`v_π(α)=1`, and `val(j) = v_π( C(m,j)M_j α^j )`. `V` and `J*` are the supporting line of slope `−1/2`
on the 2‑adic Newton polygon of `P_λ`, restricted to the on‑line coefficients.

**1.3 The evenness is *not* formal.** The only thing ramification gives is that the slope‑`(−1/2)`
**edge length** `j_max − j_min` is even — but that is just §1.1 (same parity), so it is content‑free.
`|J*|` counts the coefficients **on** the line, which can skip interior lattice points (e.g.
`λ=(3,3,1,1)`, `m=4`: `J* = {0,4}` but `j=2` is strictly above the line). Concretely, `|J*|` is the
number of nonzero terms of the **residual polynomial** `e(y) = Σ_{t∈T} y^t ∈ 𝔽₂[y]` (`T = J*` in
half‑coordinates), and

> `|J*|` even  ⟺  `e(1) = 0` in `𝔽₂`  ⟺  `(1+y) | e(y)`.

But `e(1) = |J*| mod 2` is exactly the leading‑π residue, so **this reformulation is empty by itself**:
evenness must come from the *arithmetic of the `M_j`*, never from Newton‑polygon formalism alone.

**1.4 The binomial skeleton is trivial.** Let `φ(j) = j + 2 v₂ C(m,j)`. Then `argmin_j φ = {0}` for
every `m` (since `φ(j) = 0 ⟺ j = 0`). *Consequence:* binomials never manufacture a tie; **every tie is
created entirely by `v₂(M_j)`.** This is the precise form of the "4‑core obstruction" and explains why
all earlier purely‑binomial attacks could not see the phenomenon. (Verified, and immediate from
Kummer.)

**1.5 The exact lift and conjugation.** With `Φ_λ(z) = ⟨s_λ,(p₁²+z p₂)^m⟩`:
`Φ_λ(z) = Σ_j C(m,j) M_j (−2z)^j (1+z)^{m−j}`, `G_λ(i) = Φ_λ(−i)/(1−i)^m`, and
`Φ_{λ'}(z) = Φ_λ(−z)` (because `ω(p₁²+z p₂) = p₁² − z p₂`), recovering `G_{λ'} = i^m \overline{G_λ}`.

**1.6 Sharp conjecture (verified `m ≤ 10`, 0 exceptions).** `e(y)` is always a **power of `(1+y)`**:
`e(y) = y^{t_min}(1+y)^k`. Equivalently `T` is a 2‑adic box (translate of the submask set of some `k`),
so `|J*| = 2^{popcount(k)}` — always a **power of 2**, hence even when `≥ 2`. The proof of even (let
alone power‑of‑2) for general λ remains open; see §5.

---

## 2. The hook closed form

> **Lemma 1.** For the hook `λ = (a,1^b)`, `n = 2m = a+b` (`a ≥ 1`, `b ≥ 0`) and `0 ≤ j ≤ m`,
> `  M_j = C(2m−1−j, a−1).`
> Equivalently, the generating identity
> `  Σ_j C(m,j) M_j x^j = ⟨s_λ,(h₁²+x e₂)^m⟩ = [z^{a−1}] (1+z)^{m−1}(1+x+z)^m.`

**Proof.** Since `h₁² = s₂+s_{11}` and `e₂ = s_{11}`, `h₁²+x e₂ = h₂ + (1+x) e₂` (using `s₂=h₂`,
`s_{11}=e₂`). Put `c = 1+x` and `F = (h₂ + c e₂)^m`, a symmetric function of degree `2m`; then
`Σ_j C(m,j) M_j x^j = ⟨s_λ, F⟩`.

*Step 1 — evaluate `⟨F, h_a e_q⟩` by the reproducing kernel.* For the Hall inner product,
`⟨F, H(s)E(t)⟩ = Σ_{a,q} ⟨F, h_a e_q⟩ s^a t^q`, and `H(s)E(t) = exp(Σ_k p_k(s^k − (−t)^k)/k)` is the
Cauchy kernel for the virtual alphabet with `p_1 = s+t`, `p_2 = s²−t²`. Evaluating `F` there:
`h₂ = (p₁²+p₂)/2 = s(s+t)`, `e₂ = (p₁²−p₂)/2 = t(s+t)`, hence
`h₂ + c e₂ = (s+t)(s+ct)`, so `F = (s+t)^m (s+ct)^m =: Π(s,t)`. Therefore
`  ⟨F, h_a e_q⟩ = [s^a t^q] (s+t)^m (s+ct)^m.`

*Step 2 — hook coefficient by dual Jacobi–Trudi.* `s_{(a,1^b)} = Σ_{i=0}^b (−1)^i h_{a+i} e_{b−i}`, so
`⟨s_λ, F⟩ = Σ_{i=0}^b (−1)^i [s^{a+i} t^{b−i}] Π`. As `Π` is homogeneous of degree `2m`, set `t=1`:
`[s^p t^{2m−p}]Π = [s^p] Q(s)`, `Q(s) = (1+s)^m (c+s)^m`. With `p = a+i` this gives
`  A_a := ⟨s_{(a,1^b)}, F⟩ = Σ_{p≥a} (−1)^{p−a} [s^p] Q(s).`

*Step 3 — telescoping.* From the definition, `A_a + A_{a+1} = [s^a] Q(s) = [s^a](1+s)^m(c+s)^m`. The
claimed value `B_a := [z^{a−1}](1+z)^{m−1}(c+z)^m` satisfies the same recursion:
`B_a + B_{a+1} = [z^a](1+z)·(1+z)^{m−1}(c+z)^m = [z^a]Q`. Both vanish for `a > 2m` and
`A_{2m} = B_{2m} = 1`. By downward induction `A_a = B_a` for all `a`. Extracting `[x^j]` from
`B_a = [z^{a−1}](1+z)^{m−1}(1+x+z)^m` gives `C(m,j) M_j = C(m,j) C(2m−1−j,a−1)`, i.e.
`M_j = C(2m−1−j, a−1)`. ∎

**Corollary (bonus).** Setting `x = i−1` (so `c = i`) in Lemma 1:
`  G_λ(i) = [z^{a−1}] (1+z)^{m−1} (z+i)^m`  for every hook.
(Verified `≠ 0` for all hooks `m < 12`, consistent with "(2,2) is the unique vanisher" — `(2,2)` is
not a hook.)

---

## 3. Reduction of the hook tie to a Kummer inequality

Fix a hook. Lemma 1 + `M₀ = C(2m−1,a−1)` give, for any `j` with `C(m,j)M_j ≠ 0` (so `j ≤ min(m,b)`):

> **Proposition 2.** `g(j) := val(j) − val(0) = j + 2[ v₂C(m,j) + v₂C(b,j) − v₂C(2m−1,j) ].`

**Proof.** `val(j)−val(0) = j + 2[ v₂C(m,j) + v₂C(2m−1−j,r) − v₂C(2m−1,r)]` with `r = a−1`. Falling‑
factorial cancellation:
`C(2m−1,r)/C(2m−1−j,r) = (2m−1)!(b−j)! / (b!(2m−1−j)!) = C(2m−1,j)/C(b,j)` (using `2m−1−r=b`,
`2m−1−j−r=b−j`). Hence `v₂C(2m−1,r) − v₂C(2m−1−j,r) = v₂C(2m−1,j) − v₂C(b,j)`. Substitute. ∎

Write `c₁ = v₂C(m,j)`, `c₂ = v₂C(b,j)`, `c₃ = v₂C(2m−1,j)`, so `g(j) = j − 2(c₃ − c₁ − c₂)`.

**A uniform simplification of `c₃`.** `2m−1` is odd. For `j = 2t` (even), the bit‑0 column `0+1`
produces no carry, so the carries of `2t + (2m−1−2t)` equal those of `t + (m−1−t)`:
`  c₃ = v₂C(2m−1,2t) = v₂C(m−1,t).`
For `j = 2s+1` (odd), again bit‑0 gives `1+0` with no carry, so
`  c₃ = v₂C(2m−1,2s+1) = v₂C(m−1,s).`
In both cases `c₃ = v₂ C(m−1, ⌊j/2⌋)`.

Because `g(j) ≡ j (mod 2)`, **for odd `j`, `g(j)` is odd**, so `g(j) ≥ 1` is equivalent to `g(j) > 0`;
no odd `j` can tie with `j = 0`. The whole theorem therefore comes down to controlling `c₃ − c₁`.

> **Lemma K (even case).** For `1 ≤ t ≤ m/2`:  `v₂C(m−1,t) − v₂C(m,2t) ≤ s₂(t) ≤ t`.

**Proof.** Two standard identities:
`C(m−1,t) = C(m,t)·(m−t)/m`, so `v₂C(m−1,t) = v₂C(m,t) + v₂(m−t) − v₂(m)`;
`C(m,2t)C(2t,t) = C(m,t)C(m−t,t)` with `v₂C(2t,t) = s₂(t)` (Kummer: carries of `t+t`), so
`v₂C(m,2t) = v₂C(m,t) + v₂C(m−t,t) − s₂(t)`. Subtracting,
`  c₃ − c₁ = v₂C(m−1,t) − v₂C(m,2t) = s₂(t) + [ v₂(m−t) − v₂(m) − v₂C(m−t,t) ].`
It remains to show the bracket `≤ 0`, i.e. `v₂(m−t) − v₂(m) ≤ v₂C(m−t,t)`.

- If `v₂(m−t) ≤ v₂(m)`: LHS `≤ 0 ≤ v₂C(m−t,t)`. ✓
- If `v₂(m−t) > v₂(m)`: from `m = (m−t)+t` with unequal valuations, `v₂(t) = v₂(m)`. Then
  `C(m−t,t) = (m−t)/t · C(m−t−1,t−1)` gives `v₂C(m−t,t) = v₂(m−t) − v₂(t) + v₂C(m−t−1,t−1)
  ≥ v₂(m−t) − v₂(m)` = LHS. ✓

So `c₃ − c₁ ≤ s₂(t) ≤ t`. ∎

> **Lemma K′ (odd case).** For `j = 2s+1 ≤ m` (`s ≥ 0`):  `v₂C(m−1,s) − v₂C(m,2s+1) ≤ s`.

**Proof.** As above, `v₂C(m−1,s) = v₂C(m,s) + v₂(m−s) − v₂(m)`, and
`C(m,2s+1)C(2s+1,s) = C(m,s)C(m−s,s+1)` with `v₂C(2s+1,s) = s₂(s+1) − 1`
(Kummer: carries of `s + (s+1)`, using `s₂(2s+1) = s₂(s)+1`). Hence
`v₂C(m,2s+1) = v₂C(m,s) + v₂C(m−s,s+1) − s₂(s+1) + 1`. Writing
`C(m−s,s+1) = (m−s)/(s+1)·C(m−s−1,s)`,
`  c₃ − c₁ = v₂(s+1) + s₂(s+1) − 1 − v₂(m) − v₂C(m−s−1,s) ≤ v₂(s+1) + s₂(s+1) − 1.`
For any `w ≥ 1`, `v₂(w) + s₂(w) ≤ w` (small `w` directly; in general both terms are `≤ log₂w + O(1)`),
so with `w = s+1`, `c₃ − c₁ ≤ (s+1) − 1 = s`. ∎

---

## 4. The hook theorem

> **Theorem.** For every hook `λ = (a,1^b) ⊢ 2m` (`b ≥ 0`): `0 ∈ J*` and `J* ⊆ {0,2}`. In particular
> `|J*| ∈ {1,2}`, so **`|J*|` is even whenever `|J*| ≥ 2`** (and the leading‑π layer cancels on every
> hook tie). The unique possible tie is `J* = {0,2}`, occurring iff `m` is odd and `b ≡ 2 (mod 4)`
> or `b ≡ 3 (mod 4)`.

**Proof.** It suffices to show `g(j) ≥ 0` for all `j` (so `val(0)` is the minimum and `0 ∈ J*`), with
`g(j) > 0` for every `j ∉ {0,2}`. We use `g(j) = j − 2(c₃ − c₁ − c₂)` and `c₂ = v₂C(b,j) ≥ 0`.

*Odd `j`.* By Lemma K′, `c₃ − c₁ ≤ s = (j−1)/2`, so `c₃ − c₁ − c₂ ≤ (j−1)/2` and
`g(j) ≥ j − (j−1) = 1 > 0`.

*Even `j = 2t`, `t ≥ 2`.* By Lemma K, `c₃ − c₁ ≤ s₂(t)`, and `s₂(t) ≤ t−1` for `t ≥ 2`. Hence
`c₃ − c₁ − c₂ ≤ t−1` and `g(j) ≥ 2t − 2(t−1) = 2 > 0`.

*`j = 2` (`t = 1`).* Lemma K gives `c₃ − c₁ ≤ s₂(1) = 1`. Explicitly `c₃ − c₁ = 1 − v₂(m)` and
`c₂ = v₂(b)+v₂(b−1)−1`, so `g(2) = 2[ v₂(m) + v₂(b) + v₂(b−1) ] − 2 ≥ 0`, since for `b ≥ 2` one of
`b,b−1` is even (`v₂(b)+v₂(b−1) ≥ 1`); for `b ≤ 1`, `M₂ = 0` and `j=2` is excluded.

Thus the minimum of `val` is `val(0)`, attained only at `j ∈ {0,2}`. The tie `g(2)=0` happens iff
`v₂(m)+v₂(b)+v₂(b−1) = 1`, i.e. `m` odd and the even member of `{b,b−1}` is `≡ 2 (mod 4)`. ∎

This is a complete, unconditional proof for an infinite family. It also *re‑derives* non‑vanishing on
hooks two independent ways: directly via `|J*| odd ⇒ G ≠ 0` when `J* = {0}`, and via the Corollary of
Lemma 1 (closed form for `G_λ(i)`).

---

## 5. Verification

- `M_j = C(2m−1−j,a−1)` and the GF identity of Lemma 1: brute‑force symmetric‑function check, all
  hooks `m ≤ 7`; GF identity symbolic in `x` confirmed `m ≤ 7`. ✓
- Proposition 2 identity and `g(j) ≥ 0` (`>0` off `{0,2}`): all hooks `m ≤ 200`, 0 violations. ✓
- Lemma K (`v₂C(m−1,t) − v₂C(m,2t) ≤ t`, equality `⟺ t=1`) and the bound `c₃−c₁−c₂ ≤ ⌊j/2⌋` with
  equality at even `j` forcing `j=2`: all `m ≤ 400` (Lemma K to `m ≤ 3000`), 0 violations. ✓
- `|J*| ≤ 2` for all hooks `m < 60`; `G_hook ≠ 0` for all hooks `m < 12`. ✓
- (General λ, for context.) `e(y)` a power of `(1+y)` / `|J*|` a power of 2: all `λ ⊢ 2m`, `m ≤ 10`,
  267 ties at `m=10`, 0 exceptions. ✓

---

## 6. What remains open (general λ), stated precisely

The hook proof works because `M_j` has the closed form `C(2m−1−j,a−1)`, turning `val(j)−val(0)` into a
pure Kummer expression (Prop. 2) bounded by two clean binomial identities. For general `λ`:

1. **No closed form for `M_j`.** `M_j = Σ_{μ⊢2(m−j)} (#vertical‑2‑strip chains λ↓μ)·f^μ`, and
   `v₂(M_j)` has no global closed form (the 4‑core obstruction). By §1.4 the entire tie phenomenon
   lives in these valuations.
2. **The precise open statement** (sharper than the requested theorem, equivalent to it on the tie set):
   the residual polynomial `e(y)` of the slope‑`(−1/2)` edge of `P_λ` is a power of `(1+y)` in
   `𝔽₂[y]`; equivalently the on‑line set `T` is a 2‑adic box. By §1.3 this cannot follow from
   Newton‑polygon/ramification formalism; it must be a symmetry of the chain model. The hook case is
   the instance where that symmetry is visible as the Kummer identities of §3.
3. **Most promising general route.** Extend the falling‑factorial cancellation of Prop. 2: find, for
   each `λ`, a "comparison" generating function (a `λ`‑deformation of `(1+z)^{m−1}(1+x+z)^m`) whose
   2‑adic Newton edge reproduces `T`. Two‑row shapes `(a,b)` are the natural next family — they have a
   trinomial‑type generating function and should admit a Prop‑2 analogue.

---

### Files
`code:` `/tmp/jstar.py`, `/tmp/box.py`, `/tmp/hook.py`–`/tmp/hook5.py`, `/tmp/lemK.py`, `/tmp/final.py`
(to be migrated into `projects/code/`).
