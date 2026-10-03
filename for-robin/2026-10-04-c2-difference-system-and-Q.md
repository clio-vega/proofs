# Condition (A) is a lattice count on an alcoved polytope — and the nest's two ends are where it goes affine

**PROVE 2026-10-04 cycle 2.** Paper: `proofs/2026-10-04-c2-difference-system.tex` (compiles).
Code: `proofs/code-1004c2-nest/`. Registry: 6 new nodes + 6 children under
`conj-A-logconcave` in `proofs/registry/cylindric-lorentzian.json`, trustcheck exit 0.

## The one thing to look at

Translating the three horizontal-strip conditions of the chain
`μ ⊆ ν ⊆ κ ⊆ λ` into bead coordinates `y_i = ν_i − (λ_{i−1}+1)`, `r_i = κ_i − (λ_{i−1}+1)`
gives, with `g_i = λ_i − λ_{i−1} − 1` and `σ = Σ y_i`:

```
k(a,b) = # { (y,r) ∈ Z^m × Z^m :  y ∈ B,  Σ y_i = σ,
                                  0 ≤ r_i ≤ g_i,
                                  y_i ≤ r_i,  r_i − y_{i+1} ≤ g_i,
                                  Σ r_i = a + σ }
```

Apart from the two coordinate-sum equalities, **every constraint bounds a single
coordinate or a difference of two coordinates**. So the ambient set is the lattice-point
set of an alcoved polytope, closed under componentwise max and min, and condition (A) is

> the fibre counts of an alcoved lattice set, along one of two coordinate-sum
> functionals, are log-concave.

Two things I had proved the hard way become one-liners here: **concentricity** is
`Σ_i ((y_i)₊ − (−y_i)₊) = σ` and the **ℓ¹ half-width formula** is
`Σ_i ((y_i)₊ + (−y_i)₊) = ‖y‖₁`. Both are the scalar identity `(t)₊ − (−t)₊ = t`.
Also: the support **shift** of each `Trap_ν` is exactly the horizontal-strip **defect**
`Σ_i(−y_i)₊`, which is why the summands are concentric at all.

Calibrated before anything else: exact on **155 402** slices against the box-slice
formula and **2491** against an independent chain enumeration using neither `L_i,R_i`
nor `y,r`; two refusal controls bite (1231/2493, 1956/3065).

## What it buys: the nest's two ENDS

Shift is constant on a slice **iff** there is a single half-width. At the bottom of the
nest (`defect ≡ 0`, i.e. `y_i ≥ 0` throughout) the width map is **affine**,
`w_i = g_i + 1 − y_i`, so the width-vector set is a **box slice** — M-convex, single
half-width, exactly the regime where the briefed hypothesis (H1) is *satisfiable*.
Dually at the top (`y_i ≤ 0`). So on those slices condition (A) is precisely

> **(Q)**  `F(z) = Σ_{w ∈ Box ∩ {Σ w_i = T}} Π_i [w_i]_z ∈ PF₂`,  `[w]_z = 1+z+…+z^{w−1}`

equivalently: `#{(x,e) ∈ Z_{≥0}^{2m} : l_i ≤ x_i+e_i ≤ h_i, Σx = a, Σe = D−a}` is
log-concave in `a`. Measured reach: one-sided on **12 665** of 25 768 slices; and `W` is
a box slice on **exactly** the 20 322 single-half-width slices, all of which (Q) covers.

## What is proved about (Q)

- **Interval support is free**: `supp F = [0,D]` exactly, since every summand has that support.
- **`Ω = {(x,e) ≥ 0 : l_i ≤ x_i+e_i ≤ h_i, Σ = D}` is M-convex**, elementary, two cases.
  The only hypothesis consumed is that the constraint family (singletons, the pairs
  `{x_i,e_i}`, the ground set) is **laminar** — and crossing (non-laminar) bounds refuse
  517/897, 1817/1987, 2128/2899 against 0/3227 for the laminar family.
- **(Q) at m=2, by CONCAVITY** — strictly stronger than log-concavity:
  `F(a) = Σ_j n_j (min(M(a),j)+1)` with `M(a)=min(a,D−a)` concave. 1456/1456.
  The mechanism **stops at m=3 for an exact reason**: at `L=U=(3,3,3)`, `T=9`, `W` is a
  *single point* and `F = [3]_z³ = (1,3,6,7,6,3,1)` with `1+6 > 2·3`. A single `Trap`
  stops being concave at m≥3, so the argument is intrinsically two-bead.
- **The free-D version is PF₂ outright** from `pf2-convolution`:
  `Σ_{w ∈ Box} Π_i [w_i]_z = Π_i (Σ_{w=L_i}^{U_i} [w]_z)`, each factor a positive
  nonincreasing concave sequence. 271 440 cases. **So the entire difficulty of (Q) is
  the single constraint `Σ w_i = T`** — the second hyperplane.

(Q) itself: **20 471 exhaustive + 60 000 random** instances, `m ≤ 6`, `|W| ≤ 32 765`,
`deg F ≤ 59`, zero failures.

## The question for you

(Q) at `m ≥ 3` is open, and I have two routes with a named obstruction in each.

1. **Exchange injection.** For `z ∈ K_{a−1}`, `z' ∈ K_{a+1}`, pick `p` with `x'_p > x_p`
   and `q` with `e'_q < e_q` (both exist), taking `p=q` when possible. Then
   `z + e_{x_p} − e_{e_q}` and `z' − e_{x_p} + e_{e_q}` are **both** in `K_a` — proved,
   and re-checked case by case on 805 269 pairs. **Only injectivity is missing**: the
   transfer is determined by `(z,z')` but not recoverable from the image. The natural
   repair (a bipartite degree count) is **refuted**: 7556 of 10 427 configurations have
   min-degree < max-degree, smallest witness `m=2`, `(l,h)=((0,0),(1,3))`, `D=3`, `a=1`.
   *Is there a canonical choice of `(p,q)`, or a known normality/IDP argument that gives
   log-concavity of fibre counts of a laminar polymatroid base set?*

2. **Lorentzian, and this one is blocked by bookkeeping, not mathematics.** Put
   `P̃_i(X,E,W) = Σ_{x+e+w=h_i, l_i ≤ x+e ≤ h_i} X^x E^e W^w` — all coefficients 1, support
   a box ∩ hyperplane hence M-convex — so `k(a) = [X^a E^{D−a} W^{H−D}] Π_i P̃_i`, and by
   my own proved `raw-hessian-lemma` the `{X,E}` 2×2 Hessian minor at
   `β = (a−1, D−a−1, H−D)` **is literally** `k(a)² ≥ k(a−1)k(a+1)`. So (Q) follows if
   `N(Π_i P̃_i)` is Lorentzian, which follows from Brändén–Huh 1902.03719
   Cor `normalizedcoefficients` (log-coefficients M-concave ⇒ `N(f)` Lorentzian) plus
   Cor `CorollaryConvolution` (`N(f),N(g)` Lorentzian ⇒ `N(fg)` Lorentzian).

   **Why I did not claim it.** `memory/reading/sources.json` grades 1902.03719
   `deep-read`, but its `locators` field — the field that holds that read's output —
   lists eight labels and *none* of them is `normalizedcoefficients`,
   `CorollaryConvolution` or `flow`. Those three appear only in registry-node prose, and
   the 2026-09-30 WAKE adjudication inside that same `sources.json` entry says in terms
   that closure theorems are uncited until settled against the source. The local `.tex`
   is not on disk and this was a no-browsing session. **If you can confirm those two
   corollaries, (Q) is proved and with it condition (A) on every single-half-width slice.**

   I verified the chain independently of the citation, using `raw-hessian-lemma` with the
   lean-verified `l3-det-reduction`, exact integer arithmetic: `N(P̃_i)` Lorentzian on all
   66 pairs `0 ≤ l ≤ h ≤ 10`, `N(P̃₁P̃₂)` on all 441 cases, `N(P̃₁P̃₂P̃₃)` on all 1000
   cases; refusal control (delete one interior support point) refuses 5/5.

   One wrong version, recorded so nobody re-walks it: applying `normalizedcoefficients`
   *directly* to `Ω ⊂ Z^{2m}` and aggregating gives log-concavity of `Σ_{z ∈ K_a} 1/z!`,
   **not** of `|K_a|` — the `α!` normalisation weights the count. The per-bead
   homogenisation above is the version that avoids it.

## Two routes I closed, so you need not

- **Dual form.** Summing over `y` first gives `k(a,b) = Σ_{Σr = a+σ} N(r)` for one fixed
  function `N` on a box — the slice-sum transform. The obvious hypothesis "`N` is
  M♮-concave" is **false**: 8 of 1842 real slices fail, smallest at `m=4`, `n=7`,
  `μ=(0,1,3,5)`, `λ=(1,2,4,6)`, `b=1` — and its slice sums `(1,5,5,1)` are still PF₂.
- **The brief's own first move** was the third consecutive vacuous pass. Neither
  `thm:blind` witness is *realisable* in the class the hypothesis quantifies over: at
  `m=2`, `w₁ − w₂ = (g₁−g₂) + σ − 2y₁` is injective on a slice, so realisable width
  multisets have pairwise distinct differences, and witness A has `(0,0,−4,−4)`, B has
  `(0,0,0,0)`. 0 realisations among 43 008 exhaustively enumerated abstract slices.

## Honest notes on my own instruments

Three of my controls this session were not controls, and I would rather you heard it from me:

- C3–C5 on the abstract conjecture never refuse — because each perturbation **preserves
  concentricity**, so each is another *instance* of the statement, not a perturbation.
  That is also how I found the right level of generality.
- "Drop the pair structure" for Theorem `Ω`-M-convex produced a **box ∩ hyperplane**,
  which is M-convex: 1037/1037 green, and I briefly read it as confirmation.
- The two-hyperplane control returned `0/6000` that was **meaningless** —
  `centred_sum` was silently keeping one parity class of the half-integer grid and
  dropping the other. Guard added; the control now returns *undefined*.
