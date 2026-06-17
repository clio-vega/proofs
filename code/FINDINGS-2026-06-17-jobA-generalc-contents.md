# FINDINGS — Job A (CODE 2026-06-17): general-`c` boundary content table `g[c][k]`

**Headline verdict.** The 2-adic content `g[c][k]` of the boundary deficits is **NOT
`c`-uniform** as an equality — at deep depth `k≥4` it grows with `c`. **But** there is a
clean `c`-uniform-up-to-parity **lower bound**

> **`g₀(k) = 2⌊k/2⌋ + [c odd ∧ k odd]`**

that is a **valid floor for all `c≥4`** (0 violations) and **exactly sharp for shallow
depth `k≤3`** — which is the only regime the boundary lemma actually consumes. The deep
surplus (`true − g₀`, growing with `c`) is unused slack. So the PROVE.md "Content Conjecture
`g(k)` uniform in `c`" is **FALSE as an equality** (confirming the 06-17 residual) but **TRUE
and provable as the lower bound `g₀`**. This is what the prover should carry.

Engine: `generalc_contents.py` (reuses `threerow-boundary/mn.py`, `generalc_master.py`,
`theta_scout.py`). Two cross-checked methods — see the **convention note** (§5), which is
itself a finding.

`λ=(a,b,c) ⊢ 2m`, `a=2P+b+c`, `m=P+b+c`. Box interior `b≥2c`. `val(j)=j+2(v₂C(m,j)+v₂M_j)`,
`Δ(j)=val(j)−val(0)`. Depth `k=c−i` for boundary index `j=b+i` (`k=0` is the top `N_c=1`).

---

## 1. The content table `g[c][k]` (the deliverable)

Min over `b`-parity, box interior `b≥2c`, `m≤84` (mins hit on 100s of shapes — true
content, not a grid edge):

```
        k=0  k=1  k=2  k=3  k=4  k=5  k=6  k=7
  c=2    0    0
  c=3    0    1    1
  c=4    0    0    2    2
  c=5    0    1    2    3    4
  c=6    0    0    2    2    5    5
  c=7    0    1    2    3    5    7    6
  c=8    0    0    2    2    6    6    8    8
```

`b`-even / `b`-odd split:

```
  b even:                          b odd:
        k0 k1 k2 k3 k4 k5 k6 k7           k0 k1 k2 k3 k4 k5 k6 k7
  c=4    0  1  2  4                  c=4   0  0  2  2
  c=5    0  1  2  3  4               c=5   0  2  2  4  4
  c=6    0  1  2  3  5  7            c=6   0  0  2  2  5  5
  c=7    0  1  2  4  5  7  6         c=7   0  1  2  3  5  7  7
  c=8    0  1  2  4  6  7  8 11      c=8   0  0  2  2  6  6  8  8
```

## 2. `c`-uniformity verdict — NOT uniform, but `g₀` is a sharp floor

Per-depth values across `c` (min over parity):

| `k` | values `(c, g)` | uniform? |
|---|---|---|
| 0 | all 0 | **UNIFORM** |
| 1 | `c` even→0, `c` odd→1 | `c mod 2` |
| 2 | 1 (c=3), else 2 | c-dependent (small-c) |
| 3 | even→2, odd→3 | `c`-parity |
| ≥4 | grows: e.g. `k=4`: 4,5,5,6 (c=5,6,7,8) | **`c`-DEPENDENT (grows)** |

The candidate floor `g₀(k)=2⌊k/2⌋+[c odd ∧ k odd]`:

- **valid lower bound** (`true ≥ g₀`) for **all `c≥4`, all depths** — 0 violations;
- **exactly sharp** (`true = g₀`) for **all shallow depth `k≤3`, `c≥4`**;
- deep deviations `true−g₀ ≥ 0` grow with `c` (e.g. `c=8`: `+2` at `k=4,5,6,7`).

(`c=3` is the lone small-`c` exception, `g[3][2]=1<g₀=2`; `c=3` is independently complete.)

Reading of `g₀`: `2⌊k/2⌋` is the **even-`c` depth content** (the `b(b+1)`-even / integral
2-content peeling, `c`-uniform at shallow depth). The `+[c odd ∧ k odd]` is exactly the
**`(a−b+1)|N_i` divisibility** for odd `c` (see §5): `a−b+1=2(P+(c+1)/2)` is always even when
`c` is odd, and it divides `N_i` for `i≥2` (`k≤c−2`), donating `+1`. At `k` even / `c` even
there is no smooth `(a−b+1)` deficit, so no `+1`. This matches the c=3/c=4 closed proofs.

## 3. `θ(c)` interior offset (prover carries this symbolically)

```
  c :  a even    a odd
  2 : (θ,j0)=(0,0)  (0,0)
  3 : (0,0)         (3,3)
  4 : (0,0)         (0,0)
  5 : (0,0)         (3,3)
  6 : (0,0)         (0,0)
  7 : (0,0)         (3,3)
  8 : (0,0)         (0,0)
```

**`θ∈{0,3}`, uniform in `c`: `θ=3, j₀=3` iff `c` odd ∧ `a` odd; else `θ=0, j₀=0`.** NOT
`τ(τ+1)/2`. So the boundary requirement (`Δ(b+i) > −θ`) never gets harder as `c` grows.

## 4. Certify `Δ(b+i) > −θ` — boundary always loses, `c=4..8`

True `Δ` (not the hand bound), box interior, `m≤100`: **0 violations**, min true slack
`Δ+θ` per `(c, b-parity)`:

| `c` | `b` even (slack @ shape) | `b` odd (slack @ shape) |
|---|---|---|
| 4 | 8 @ (12,8,·,i=2) | 8 @ (13,9,·,i=3) |
| 5 | **4** @ (19,12,·,i=1) | 8 @ (24,11,·,i=1) |
| 6 | 8 @ (20,14,·,i=2) | 8 @ (25,13,·,i=1) |
| 7 | **4** @ (21,14,·,i=1) | 14 @ (22,15,·,i=1) |
| 8 | 18 @ (24,16,·,i=2) | 18 @ (27,17,·,i=3) |

Smallest slack `4` (odd `c`, `b`-even, deepest index `i=1`) — comfortable; we need only
`> 0`. The proof's **hand bound** (`generalc_certify.py`, master + Lemma P + compensation)
is also valid and sufficient `c=4..8`, `m≤110`, **0 fail**, min hand margin `2`.

## 5. Convention note on `N_i` (cross-check resolution — a finding for the prover/Lean)

The two engines disagreed at exactly **one** entry (`c=5, k=3, i=2`): symbolic `2/3` vs
direct `3/4`. **Not an error — a convention difference**, and the prover must use the right one:

- **DIRECT / MASTER convention** (matches the `c=4` memory table *and* the master `Δ`
  formula): for `i≥2`, `N_i = M_{b+i}·c!·k! / [(b−c+1)∏_{t=2}^i(b+t)]`, **keeping**
  `(a−b+1)` and `(a−c+2)` inside `N_i`; for `i=1` additionally strip `(a−b+1)`. In the
  master `Δ`, `−2v₂(a−c+2)` and `−2v₂(a−b+1)` are **separate** terms (compensated by
  `+2v₂(Π_i)` via Lemma P, and by `(a−b+1)|N_i`). **This is the content the prover needs**,
  and it is the §1 table.
- **SYMBOLIC full-strip** also removes `(a−b+1)` whenever it is a *polynomial* factor.

They agree everywhere except where `(a−b+1)|N_i` polynomially — which happens for **odd `c`,
`i≥2`** (`(a−b+1)|N_2^{(5)}`, memory-confirmed), lowering the full-strip content by `1`. **The
`+[c odd ∧ k odd]` in `g₀` is precisely this `(a−b+1)`-divisibility.** Lean must fix the
`N_i` definition to the direct/master convention or the contents shift by 1 on odd-`c` indices.

---

## Verdict for the prover / Lean

1. **Carry `g₀(k)=2⌊k/2⌋+[c odd ∧ k odd]` as the content lower bound** (master convention).
   It is `c`-uniform up to parity, sharp at the shallow depths the boundary uses.
2. **Do not aim to prove the content `=` a depth-only `g(k)`** — it is genuinely
   `c`-dependent at deep `k` (grows). You don't need that; `g₀` + Lemma P (absorbs `a−c+2`) +
   `(a−b+1)|N_i` (odd `c`) close `Δ > −θ` with slack `≥4`.
3. **`θ(c)∈{0,3}`**, `=3` iff `c` odd ∧ `a` odd; carry symbolically.
4. **Two proof obligations remain to fully generalise** (the named residual, unchanged):
   (a) the even-`c` depth content `2⌊k/2⌋` from an integral-2-content peel (`b(b+1)`-even
   type) — needs a `c`-uniform statement; (b) `(a−b+1)|N_i` for all odd `c≥7` (verified
   `c=5,6` at `i=2`; candidate rep-theoretic vanishing on `a=b−1`).

## Files
`code/generalc_contents.py`, `code/results_jobA_generalc.txt`. Reuses
`threerow-boundary/{mn,generalc_master,generalc_certify,theta_scout,generalc_content}.py`.
Updates `2026-06-17-generalc-master-and-c5` residual; supports
`proofs/2026-06-17-generalc-boundary-master-and-c5.md`.
