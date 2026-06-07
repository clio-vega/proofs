# Two-row d=4 law: a 2-adic proof for all b ≡ 1 (mod 4)

**Date:** 2026-06-07 (prove session)
**Target (from PROVE.md):** close `(♦)` — for `b ≥ 5`, the polynomial `Q_b(m)` has no
integer root `m ≥ b` — which is the sole remaining gap in the two-row case of the d=4
fiber-vanishing law `G_{(2m−b,b)}(i) = 0 ⟺ (2,2)`.

**Result.** A clean **2-adic valuation identity** proves the law outright for **all
`b ≡ 1 (mod 4)`** (an infinite family; prior cycles had only finite verification `b ≤ 24`).
The same identity is shown to hold for `b ≡ 0 (mod 4)` numerically (verified `b ≤ 40`,
random `m ≤ 3·10⁴`) and is reduced to a single clean sub-statement; and the reason the
method **must** be modified for `b ≡ 2,3 (mod 4)` is pinned exactly (those `q_b` acquire a
genuine root mod 2). This reframes the previous cycle's "unbounded `v₂` kills the 2-adic
route" verdict: the unbounded valuation lives **entirely** in the forced product
`∏(m−r)`; the quotient `Q_b` is a 2-adic *unit* of constant valuation.

---

## 0. Setup (proved last cycle; not re-derived)

For `λ = (2m−b, b)`, `0 ≤ b ≤ m`, with `s = 1+i`, `P(u) = (1+su+u²)^m`,
```
   I_b(m) := Im G_{(2m−b,b)} = Im [u^b]((1−u) P(u)) ∈ ℤ.
```
Last cycle (`2026-06-06-tworow-d4-no-rational-root`) proved: `I_b(m)` is a polynomial in `m`;
its only integer roots are the **forced** ones `m ∈ {0,1,…,R}` with `R := ⌊(b−1)/2⌋`; and
writing
```
   I_b(m) = (m)_{R+1} · Q_b(m),     (m)_{R+1} := ∏_{r=0}^{R} (m−r)   (falling factorial),
```
the two-row law is exactly:

> **(★)  `Q_b(m) ≠ 0` for every integer `m ≥ b`.**

Since `m ≥ b > R`, the forced product `(m)_{R+1} ≠ 0`; thus `(★) ⟺ I_b(m) ≠ 0` for `m ≥ b`.

---

## 1. The exact binomial expansion of `I_b(m)`

Split the trinomial as `1 + su + u² = (1+u²) + s·u` (real part `1+u²`, the entire `i`
sitting in `s = 1+i`). Then
```
   P(u) = ((1+u²) + su)^m = Σ_{j≥0} C(m,j) s^j u^j (1+u²)^{m−j},
```
so, taking imaginary parts (only `s^j = (1+i)^j` is complex),
```
   Im P(u) = Σ_{j≥1} C(m,j) · Im((1+i)^j) · u^j (1+u²)^{m−j}.
```
Applying `Im [u^b]((1−u)·—)` and writing `T(l) := [u^l](1+u²)^{m−j} = C(m−j, l/2)` for `l`
even (and `0` for `l` odd), we obtain the **exact** finite sum

> **(E)**  `I_b(m) = Σ_{j=1}^{b} τ_j(m)`,  where
> ```
>    τ_j(m) = C(m,j) · Im((1+i)^j) · ( T(b−j) − T(b−1−j) ).
> ```

Of the two extractions `T(b−j), T(b−1−j)` exactly one is nonzero (their arguments
`b−j, b−1−j` are consecutive, so exactly one is even); hence
`T(b−j) − T(b−1−j) = ± C(m−j, β_j)` with `β_j := ⌊(b−j)/2⌋`.

**Valuation of the complex factor.** `Im((1+i)^j) = 2^{j/2} sin(jπ/4)`, an integer, with
```
   v₂(Im((1+i)^j)) = ⌊j/2⌋   (j ≢ 0 mod 4),    Im((1+i)^j) = 0   (j ≡ 0 mod 4).
```
In particular **`Im((1+i)^j)` is odd iff `j = 1`** (the sequence is `0,1,2,2,0,−4,−8,−8,0,…`).

### 1a. Corollary (the mod-2 law), rigorous for every `b`

Reducing (E) mod 2, every `j ≠ 1` term is even, so the `j=1` term survives:
`τ_1(m) = C(m,1)·1·(T(b−1)−T(b−2))`. For `b` odd this is `m·C(m−1,R)`; for `b` even it is
`−m·C(m−1,R)`; in both cases `|τ_1(m)| = m·C(m−1,R) = (m)_{R+1}/R! = (R+1)·C(m,R+1)`. Hence

> **(M2)**  `I_b(m) ≡ (R+1)·C(m, R+1)  (mod 2),    R = ⌊(b−1)/2⌋.`

This already shows `I_b(m)` even whenever `(R+1)C(m,R+1)` is even, recovering the prior
"`Im G` even for `m` even" and the `b ≡ 0,3 (mod 4)` evenness, with a one-line proof.

---

## 2. Main Theorem: the law for all `b ≡ 1 (mod 4)`

> **Theorem.** Let `b ≡ 1 (mod 4)`, `b ≥ 5`, and `R = (b−1)/2`. Then for **every** integer `m`,
> ```
>      v₂(I_b(m)) = v₂((m)_{R+1}) − v₂(R!).                         (V)
> ```
> Consequently `v₂(Q_b(m)) = −v₂(R!)` is constant; `Q_b(m)` is a 2-adic unit times
> `2^{−v₂(R!)}`, in particular **`Q_b(m) ≠ 0` for all integer `m`**, so `(★)` and the two-row
> d=4 law hold for all `b ≡ 1 (mod 4)`.

The proof is: the `j=1` term of (E) has valuation exactly the RHS of (V), and **every other
term has strictly larger `v₂`**, so it cannot disturb the bottom 2-adic digit.

### 2.1 The `j=1` term is exact

For `b ≡ 1 (mod 4)`, `b` is odd, `b−1 = 2R` is even and `b−2` is odd, so `T(b−2)=0` and
```
   τ_1(m) = m·C(m−1,R) = (m)_{R+1}/R!,    v₂(τ_1(m)) = v₂((m)_{R+1}) − v₂(R!).   (= RHS of V)
```

### 2.2 Reducing each `j ≥ 2` term to an `m`-independent inequality

Fix `2 ≤ j ≤ b`, `j ≢ 0 (mod 4)` (else `τ_j ≡ 0`). Put `h := ⌊j/2⌋`, `d_j := ⌊(j−1)/2⌋`,
`β_j = ⌊(b−j)/2⌋`, and `s_j := j + β_j`. For `b` odd one computes `s_j = ⌊(b+j)/2⌋` and
`d_j = s_j − (R+1)`. The Vandermonde/multinomial identity
```
   C(m,j)·C(m−j, β_j) = C(m, s_j)·C(s_j, j)
```
turns the term valuation into
```
   v₂(τ_j(m)) = v₂(C(m, s_j)) + v₂(C(s_j, j)) + ⌊j/2⌋.
```
Telescoping `C(m, s_j)/C(m, R+1) = ∏_{t=1}^{d_j} (m−R−t)/(R+1+t)` gives
```
   v₂(C(m,s_j)) = v₂(C(m,R+1)) + Σ_{t=1}^{d_j} v₂(m−R−t) − Σ_{t=1}^{d_j} v₂(R+1+t).
```
Since for `b ≡ 1 (mod 4)`, `R+1 = (b+1)/2` is **odd** (so `v₂(τ_1) = v₂(C(m,R+1))`), subtracting
`v₂(τ_1(m)) = v₂(C(m,R+1))` yields the clean splitting

> **(Δ)**  `v₂(τ_j(m)) − v₂(τ_1(m)) = D_j + Σ_{t=1}^{d_j} v₂(m−R−t)`,
> where the **`m`-independent** part is
> ```
>    D_j := v₂(C(s_j,j)) + ⌊j/2⌋ − Σ_{t=1}^{d_j} v₂(R+1+t).
> ```

### 2.3 A closed form for `D_j`, and `D_j ≥ 0`

Using `s_j − j = R − h` (valid for `b` odd), `v₂((2h)!) = v₂((2h+1)!) = h + v₂(h!)`, and
`Σ_{t=1}^{d_j} v₂(R+1+t) = v₂((R+1)!) − v₂((R−h)!)` (note `d_j = h` if `j=2h+1`, and the same
final identity holds for `j=2h` because `s_j−j = R−h`):
```
   D_j = ⌊j/2⌋ − v₂(j!) + v₂((R+1)!) − v₂((s_j−j)!)
       = h − (h + v₂(h!)) + v₂( ∏_{i=R−h+1}^{R+1} i )
       = v₂( (h+1)! · C(R+1, h+1) ) − v₂(h!)
```
(`∏_{i=R−h+1}^{R+1} i` is a product of `h+1` consecutive integers `= (h+1)!·C(R+1,h+1)`), so

> **(D)**  `D_j = v₂(h+1) + v₂( C(R+1, h+1) ) ≥ 0,    h = ⌊j/2⌋.`

### 2.4 Strict domination, and the conclusion

Combine (Δ) and (D). For each `j ≥ 2` (with `j ≢ 0 mod 4`):

- **If `D_j ≥ 1`:** then `v₂(τ_j) − v₂(τ_1) ≥ D_j ≥ 1`.
- **If `D_j = 0`:** by (D) this forces `h+1` odd (i.e. `h` even) and `C(R+1,h+1)` odd. With
  `h` even and `j ≥ 2`:
   * `j` even `= 2h` with `h` even would give `j ≡ 0 (mod 4)` — excluded. So `j` is **odd**,
     `j = 2h+1` with `h` even, i.e. **`j ≡ 1 (mod 4)`**, and `j ≥ 5` (since `h ≥ 2`; `h=0`
     is `j=1`). Then `d_j = h ≥ 2`.
   * Among the `d_j ≥ 2` consecutive integers `m−R−1, …, m−R−d_j` at least one is even, so
     `Σ_{t=1}^{d_j} v₂(m−R−t) ≥ 1`. Hence `v₂(τ_j) − v₂(τ_1) ≥ 0 + 1 = 1`.

(The small cases are consistent: `j=2 ⇒ D_2 = v₂(C(R+1,2))+1 ≥ 1`; `j=3 ⇒ h=1`,
`D_3 = v₂(2)+v₂(C(R+1,2)) ≥ 1`.)

So **`v₂(τ_j(m)) ≥ v₂(τ_1(m)) + 1` for every `j ≥ 2` and every integer `m`.** Therefore in
`I_b(m) = τ_1(m) + Σ_{j≥2} τ_j(m)` the tail has `v₂ ≥ v₂(τ_1)+1`, and
`v₂(I_b(m)) = v₂(τ_1(m)) = v₂((m)_{R+1}) − v₂(R!)`, which is (V). ∎

**Verification.** The decomposition (Δ) with the closed form (D), and the strict inequality
`v₂(τ_j) > v₂(τ_1)`, were checked exactly for all `b ≡ 1 (mod 4)`, `5 ≤ b ≤ 41`, all valid
`j`, at 300 random `m ≤ 10⁵` per `b` (`verify_Dj.py`: "ALL … OK"). The identity (V) itself was
checked directly for `b ≤ 21` over `m ∈ [b, b+200]` plus 400 random `m ≤ 2·10⁴` (`bigm.py`).

---

## 3. The `b ≡ 0 (mod 4)` case: same identity, one remaining estimate

For `b ≡ 0 (mod 4)`, `b` even, the `j=1` term is again `τ_1(m) = −(m)_{R+1}/R!`,
`R = (b−2)/2`, but now `R+1 = b/2` is **even**, so `v₂(τ_1) = v₂(R+1) + v₂(C(m,R+1))` and the
`m`-independent part picks up a correction `−v₂(R+1) < 0`:
```
   D_j^{(0)} = D_j − v₂(R+1)            (b ≡ 0 mod 4).
```
This can be negative, so individual terms can **tie** `τ_1` in valuation (verified: ties but
never a strictly smaller term). The numerical identity nevertheless holds:

> **Proposition (verified, not yet proved).** For `b ≡ 0 (mod 4)` and every integer `m`,
> `v₂(I_b(m)) = v₂((m)_{R+1}) − v₂(R!)`, hence the two-row law holds for `b ≡ 0 (mod 4)`.
> *(Checked `b ≤ 40`; 500 random `m ≤ 3·10⁴` for `b = 8,12,16,20,24`.)*

The whole gap is the single clean statement: **the sum of the valuation-minimal terms of (E)
is a 2-adic unit times `2^{v₂(τ_1)}`** (equivalently: dividing `I_b` by `2^{v₂(τ_1)}` and
reducing mod 2 never gives `0`). Concretely this is a mod-`4` refinement of (M2): only
`j = 1,2,3` contribute mod 4 (since `Im((1+i)^j) ≡ 0 mod 4` for `j ≥ 5` and `j ≡ 0`), giving
```
   I_b(m) ≡ −m C(m−1,R) + 2 C(m,2) C(m−2,R) − 2 C(m,3) C(m−3,R−1)   (mod 4),
```
and one must show its `v₂` equals `v₂((m)_{R+1}/R!)`. This is left as the precise open piece.

Equivalent face of both cases (verified `b ≤ 40`): as a polynomial over `𝔽₂`,
```
   q_b(m) ≡ (m² + m + 1)^{⌊(b−1)/4⌋}   (mod 2)        for b ≡ 0, 1 (mod 4),
```
where `q_b` is the monic integer "irreducible part" of `Q_b`; since `m²+m+1` is odd at every
integer, `q_b` has **no root mod 2**. (For `b ≡ 1` this is the Theorem; for `b ≡ 0` it is the
Proposition.)

---

## 4. Why `b ≡ 2, 3 (mod 4)` needs a different prime (obstruction pinned)

For `b ≡ 2, 3 (mod 4)` the empirical factorization is
```
   q_b(m) ≡ m · (m² + m + 1)^{⌊(b−1)/4⌋}   (mod 2),
```
i.e. `q_b` has a **genuine root `m ≡ 0 (mod 2)`**: the mod-2 certificate is *unavailable*,
not merely unproven. Correspondingly the `j=1` domination fails — the offset
`v₂(I_b(m)) − v₂((m)_{R+1})` is **non-constant** in `m` (e.g. `b=6`: it takes values
`{−1,0,1,2,3,6}` over a short `m`-window), because the higher `j` terms now genuinely
interfere. These `b` require an **odd** prime `p` with a root-free reduction / Newton-polygon
slope obstruction (the Filaseta `Lemma 17`, `k=1` route named in PROVE.md). That is the
natural next target. The discriminant of `q_b` is non-square for all `b ≤ 30` (so `q_b` is
non-square-ramified, consistent with — but not implying — irreducibility), and `q_b` is in
fact irreducible over `ℚ` for all `b ≤ 35` (sympy); a uniform irreducibility proof would close
all residue classes at once but is the hard arithmetic core.

---

## 5. What is proved, and the precise status

**Proved (rigorous, new):**
- (M2) `I_b(m) ≡ (R+1) C(m,R+1) (mod 2)` for every `b` (§1a).
- **Theorem (§2): the two-row d=4 law for all `b ≡ 1 (mod 4)`** via the exact 2-adic
  identity `v₂(I_b(m)) = v₂((m)_{R+1}) − v₂(R!)`. This is an *infinite* family; previously
  only `b ≤ 24` was known (by direct factorization).
- The structural reframing: the unbounded `v₂(I_b)` of the prior cycle is *entirely* the
  forced product `(m)_{R+1}`; `Q_b` is a 2-adic unit of constant valuation `−v₂(R!)`. The
  "no finite 2-adic truncation closes it" obstruction was an artifact of not factoring out
  the forced roots.

**Verified, reduced to one clean estimate:**
- `b ≡ 0 (mod 4)` (§3): same `v₂` identity; gap = "minimal terms sum to a 2-adic unit"
  (a mod-4 refinement of (M2)). Verified `b ≤ 40`, random `m ≤ 3·10⁴`.

**Open (obstruction identified):**
- `b ≡ 2, 3 (mod 4)` (§4): `q_b` has a real root mod 2; the 2-adic method is structurally
  inapplicable. Needs an odd-prime Newton-polygon / irreducibility input.

**Net:** `(♦)`/`(★)` is now proved for the full residue class `b ≡ 1 (mod 4)` and reduced to
one mod-4 lemma for `b ≡ 0 (mod 4)` — i.e. the two-row d=4 law holds unconditionally for half
of all `b` up to a single clean sub-statement, and rigorously for `b ≡ 1 (mod 4)` outright.

## 6. Scripts
`~/projects/scratch/2026-06-07-prove/`: `qb.py`, `explore2.py` (factor `Q_b`, discriminants),
`modp.py` (irreducibility / `mod p`), `eisen.py` (no Eisenstein prime), `smallp.py`
(`no-root-mod-p` by class), `qmod2.py` (the `(m²+m+1)^k` factorization), `v2test.py`/`v2exact.py`
(offset constancy), `terms.py`/`bigm.py`/`verify_Dj.py` (term domination + the (Δ),(D)
decomposition), `check_b0.py` (`b ≡ 0` aggregate at large `m`).
