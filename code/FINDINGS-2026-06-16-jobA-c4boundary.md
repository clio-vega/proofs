# FINDINGS — Job A (2026-06-16 code session): de-risk the `c=4` three-row boundary

**One-line verdict.** The `c=4` three-row boundary is **comfortable, not tight**: over
**75 855 box-interior shapes to `m≤400`**, every boundary index loses by **`Δ(b+i) ≥ 8`**
(we need only strict `Δ>0`, `θ=0`). The "wall" — the irreducible a-quadratic `M_{b+2}` — is
**genuinely irreducible** (confirmed symbolically) but **illusory as an obstruction**: the
crude 2-adic **content bounds** `v₂(N₂)≥2`, `v₂(N₁)≥2` (b odd) / `≥3` (b even), `v₂(N₃)≥1`
(b even) already over-deliver against a `≥8` actual deficit, so **no factoring and no `F_3`
are needed**. Every closed form cross-checked against Murnaghan–Nakayama: **0 mismatches**
(3154 forms `m≤46`; 4900 boundary `Δ` indices `m≤60`). This **confirms the closed proof**
`proofs/2026-06-16-c4-boundary-complete.md`.

Object: `λ=(a,b,4)`, `a+b=2m−4` (**a,b same parity**, c=4 even), `P=m−b−4`, `a=2P+b+4`.
Box interior `b≥2c=8`; `val(j)=j+2(v₂C(m,j)+v₂M_j)`, `Δ(j)=val(j)−val(0)`. Need `Δ(b+i)>0`
(`θ=0`, `j₀=0`, box `{0,2,4,6}`).

---

## 1. Closed forms (SymPy, MN-validated)

| `j` | `M_j` | deficit factor `N_i` |
|---|---|---|
| `b+4` (top) | `(b−3)(b+2)(b+3)(b+4)/24` | — (Lemma T) |
| `b+3` | `(b−3)(b+2)(b+3)·N₃/24` | `N₃ = ab+4a−b²−3b−8` (**a-linear**) |
| `b+2` | `(b−3)(b+2)·N₂/48` | `N₂ = a²(b+3)(b+4) − a(2b³+9b²+21b+20) + (b⁴+2b³+13b²+16b−8)` (**a-quadratic, irreducible**) |
| `b+1` | `(b−3)(a−b+1)·N₁/144` | `N₁ = a²(b+2)(b+3)(b+4) − a(2b⁴+7b³+13b²+14b) + (b⁵−2b⁴+17b³−4b²−84b−96)` (**a-quadratic, irreducible**) |

The second subtop `N₂` and the third subtop `N₁` are **irreducible quadratics in `a` over
`ℤ[b]`** (SymPy `factor_list`), with leading coefficients `(b+3)(b+4)` and `(b+2)(b+3)(b+4)`.
This reproduces the scout's claim: the factor-in-product machinery (which needs an `a`-**linear**
factor, as at c=3) does **not** apply to `M_{b+2}` for `c≥4`. The right tool is the **content**
of `N_i` after `a=2P+b+4`, not its factorisation.

## 2. The content bounds (the engine) — 0 violations, `m≤400`

Substituting `a=2P+b+4` and peeling powers of 2 (divide by `2^v`, reduce mod 2):

| `N_i` | after `a=2P+b+4`, `/2^v` then mod 2 | proved bound | **empirical min `v₂`** |
|---|---|---|---|
| `N₃` | `≡ b` | `v₂(N₃)≥1` iff `b` even | **1** (b even); `0` (b odd) — sharp |
| `N₂` | `/2 ≡ Pb(b+1) ≡ 0` | `v₂(N₂)≥2` (all) | **2** — sharp |
| `N₁` | `/2 ≡ (P+1)b²(b+1) ≡ 0` | `v₂(N₁)≥2` (all); `≥3` (b even) | **2** (b odd, sharp); **4** (b even) |

`b(b+1)` and `(P+1)b²(b+1)` are even unconditionally, which is *why* the content lifts past the
bare gcd (which is only `2`). **The b-even `N₁` bound is even better than the memory's `≥3`:
the true minimum is `4`.** All bounds: **0 violations over 75 855 shapes** (`N₂≥2`, `N₁≥2`,
`N₁≥3` b-even, `N₃≥1` b-even).

## 3. `Δ(b+i)` slack — the decisive de-risking number

> **`min Δ(b+i) = 8` for every `i∈{1,2,3,4}`**, over all 75 855 shapes to `m≤400`.
> (We need only `Δ>0`.) Extremal shapes are all small, at the box corner:

| `i` | min `Δ(b+i)` | extremal `(a,b,m,P)` |
|---|---|---|
| 1 | 8 | `(23, 9, 18, 5)` |
| 2 | 8 | `(12, 8, 12, 0)` |
| 3 | 8 | `(13, 9, 13, 0)` |
| 4 | 8 | `(12, 8, 12, 0)` |

Contrast **c=3**, where the boundary was *tight* (`Δ(b+1)` slack `2`, F2 surplus `0` at
`(17,10,15)`). At **c=4 there is no tight shape**: the content bound certifies `Δ≥1` (per the
proof) while the truth is `Δ≥8` — **7 units of unused slack everywhere**. This is the structural
reason the `c≥4` wall is illusory: the deficit `−2v₂(a−2)`-type terms enter as **positive
surplus** `+2v₂(N_i)`, and the surplus dwarfs what is needed.

## 4. `M_{b+2}` residue law (the "wall" is rich but irrelevant)

`v₂(M_{b+2})` is **unbounded** (observed range `0…17`), so there is *no* residue bound on
`M_{b+2}` itself — confirming the wall is a real arithmetic object. Per `b mod 8`:

- `b≡3, 6 (mod 8)` **force `v₂(M_{b+2}) ≥ 1`** (M always even);
- `b≡0,1,2,4,5,7 (mod 8)` allow `v₂=0` (M can be odd).

Of the 64 classes `(a%8, b%8, P%8)`, **32 force `M_{b+2}` even, 32 allow it odd**. The point:
the proof must **not** try to bound `M_{b+2}` by residues — it should keep the deficit attached
to its consecutive product and use the **content of `N₂`** (`≥2`), which is all the `Δ≥8` budget
ever requires.

## 5. `F_k` premise threshold — `F_3` needs `Q≥8` (confirmed, but unneeded)

Among `Q` consecutive integers, drop the `k` highest-`v₂` members; worst-case remaining `v₂`:

| `Q` | F1 (k=1) | F2 (k=2) | F3 (k=3) |
|---|---|---|---|
| 4 | **1** | 0 | 0 |
| 6 | 2 | **1** | 0 |
| 8 | 4 | 2 | **1** |

So **F1 sharp at `Q≥4`, F2 sharp at `Q≥6`, F3 sharp at `Q≥8`** — the scout's `F_3` guess is
**confirmed**. But the `c=4` boundary **does not need `F_3`**: the content route (§2) plus the
`Δ≥8` slack (§3) closes it outright.

---

## Verdict for the prover / Lean

1. **Engine = content, not factoring.** Use `v₂(N₂)≥2`, `v₂(N₁)≥2`/`≥3`, `v₂(N₃)≥1` (b-even),
   each from the one-line peel `a=2P+b+4 ⟹ /2 ⟹ Pb(b+1)≡0` etc. No irreducible-quadratic
   factoring, no `F_3`.
2. **There is no tight shape at c=4.** Min `Δ=8` against threshold `0`; the proof's certified
   `Δ≥1` is loose by `≥7`. Any of the four content bounds carries the day with margin.
3. **Do not bound `M_{b+2}` by residues** — its `v₂` is unbounded (0…17); the structure lives in
   `N₂`'s content, not `M_{b+2}`'s class.
4. **Base shapes to feature:** `(12,8,12)` (mins `Δ(b+2)`,`Δ(b+4)`), `(13,9,13)` (`Δ(b+3)`),
   `(23,9,18)` (`Δ(b+1)`) — all small, all at the `P`-small box corner.
5. **General `c`:** the content+peel recipe is `c`-uniform; the only `c`-dependent inputs are the
   subtop polynomials `N_i^{(c)}` and `θ(c)`. `c=5` (odd, `θ` returns) needs `θ(5)` + the `N_i^{(5)}`
   contents — same machinery.

## Files
`code/c4_boundary_sweep.py` (forms + MN validation + content/Δ/residue sweep, `m≤400`),
`code/c4_boundary_residue.py` (sharp per-class mins, `F_k` thresholds). Reuses
`threerow-boundary/mn.py`, `c4_forms.py`. Proof: `proofs/2026-06-16-c4-boundary-complete.md`.
