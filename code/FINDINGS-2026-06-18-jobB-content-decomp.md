# FINDINGS — Job B (2026-06-18 CODE): `c!·k!/den` floor arithmetic for the general-`c` Content Lemma

**Scripts:** `threerow-boundary/jobB_content_decomp.py`, `jobB_route1_qlift_stub.py`
**Output:** `threerow-boundary/results_jobB_content_decomp.txt`
**Target:** decompose the exact 2-adic content of the master deficit `N_i^{(c)}` at deep depth
`k=c−i ∈ {4,5,6}`, `c=4..12`, both `b`-parities, against the master `c!·k!/den` factor; isolate the
residual by parity of `i`; settle the floor for the next Claim-A PROVE; q-lift route-1 stub.

Convention (matches `jobA_deep.py`): `N_i^{(c)} = M_{b+i}·c!·k!/[(b−c+1)∏_{t=2}^i(b+t)]`, `M` the
alternant. Slices `a=2P+b+c`, `b=2B`(even)/`2B+1`(odd).

## Gate (all passed)
- **Optimized alternant vs MN:** `tested=409, bad=0` (new closed-form `E^j` coefficient
  `al=(A+B−C)/2,…`; same `Malt` value as the abacus χ engine).
- **Regression vs last cycle's published content table:** **0 mismatches** (grid `B≤64,P≤64`
  reproduces the stability-verified `g[c][k]` of `FINDINGS-2026-06-18-jobA-claimA-deep.md`).
- **Pointwise identity** `g = v₂(M) + v₂(c!k!) − v₂(den)` at the content minimizer: **match=True**
  on every sampled case (e.g. `c=8,k=6` even: `8 = 1+11−4`; `c=7,k=5` even: `7 = 4+7−4`).

## Verdict (one line)

The uniform floor **`v₂(N_i) ≥ 2⌊k/2⌋` (both slices), `≥ k` (even slice)** holds with **0 violations**
to `c≤12, k≤6, b≤200` — it is the right target for the next Claim-A PROVE. But it is **NOT tight at
deep `k`** (surplus 1–4), and the dream's "odd-`i` floor-tight at `k≤3`" **does NOT extend to
`k=5,6`** — at `k≥4` every slice has surplus `≥1`, odd `i` included. No clean `(k, i mod 2)` floor
sharper than `2⌊k/2⌋` is uniform in `c`. **Build on the floor; do not chase the exact content.**

## The decomposition table `g = v₂(c!k!) + μ`, `μ = min_slice[v₂(M_{b+i}) − v₂(den)]`

```
  c  k  i i%2 | g_even  mu_e |  g_odd  mu_o | v2(c!k!)
  5  4  1   1 |     6     0 |      6     0 |     6
  6  4  2   0 |     5    -2 |      5    -2 |     7
  6  5  1   1 |     7     0 |      5    -2 |     7
  7  4  3   1 |     5    -2 |      5    -2 |     7
  7  5  2   0 |     7     0 |      7     0 |     7
  7  6  1   1 |     7    -1 |      8     0 |     8
  8  4  4   0 |     6    -4 |      6    -4 |    10
  8  5  3   1 |     7    -3 |      6    -4 |    10
  8  6  2   0 |     8    -3 |      8    -3 |    11
  9  4  5   1 |     6    -4 |      6    -4 |    10
  9  5  4   0 |     7    -3 |      8    -2 |    10
  9  6  3   1 |     8    -3 |      8    -3 |    11
 10  4  6   0 |     5    -6 |      5    -6 |    11
 10  5  5   1 |     7    -4 |      5    -6 |    11
 10  6  4   0 |     8    -4 |      8    -4 |    12
 11  4  7   1 |     5    -6 |      5    -6 |    11
 11  5  6   0 |     7    -4 |      7    -4 |    11
 11  6  5   1 |     7    -5 |      8    -4 |    12
 12  4  8   0 |     6    -7 |      6    -7 |    13
 12  5  7   1 |     7    -6 |      6    -7 |    13
 12  6  6   0 |     8    -6 |      8    -6 |    14
```

**Reading the decomposition.**
- The explicit constant `v₂(c!k!)` *over*shoots the content at deep `k`: `μ = g − v₂(c!k!) < 0`
  (down to `−7`). So the `den` factors `(b−c+1)∏(b+t)` carry **substantial positive 2-content** that
  the master valuation must credit back — the content is *not* "`c!k!` minus a small block", it is
  `c!k!` heavily discounted by `den`. This is why a pure `c!k!`-counting bound is loose.
- The residual `μ` (what a Content Lemma must lower-bound) is itself **`c`-dependent and not a clean
  function of `(k, i mod 2)`** — e.g. odd slice `k=5`: `μ_o = 0,0,−4,−2,−6,−4,−7` for `c=6..12`.

## What survives as a clean statement
- **`k=5`, even slice: `g_even = 7` uniformly for all `c=6..12`** (surplus exactly `3` over `2⌊k/2⌋`).
  The *only* fully-uniform-in-`c` content in the deep range; the other slices oscillate between two
  adjacent values (even `k=4 ∈{5,6}`, `k=6 ∈{7,8}`).
- **Floor certificate** (`§3`, `0` violations, `c≤12, k≤6, b≤200`):
  `F1: v₂(N_i) ≥ 2⌊k/2⌋` both slices; `F2: v₂(N_i) ≥ k` even slice. These remain the largest
  *uniform* bounds — confirmed still valid at `k=5,6`.

## Does "odd-`i` floor-tight" extend to `k=5,6`?  **NO.**
Surplus `g − 2⌊k/2⌋` for odd-`i`, odd slice: `k=4`→`{1,2}`, `k=5`→`{1,2}` (e.g. `c=10,i=5`: surplus 1
i.e. `g=5`, not tight), `k=6`→`{2}`. The floor-tightness the dream observed at `k≤3` is a
**small-depth** phenomenon; by `k≥4` the alternating-sum cancellation lifts content off the floor on
every parity. So the next PROVE cannot assume odd-`i` deficits sit at `2⌊k/2⌋` for `k≥4`.

## Cleanest floor statement for the next Claim-A PROVE
> Prove `v₂(N_i^{(c)}) ≥ 2⌊k/2⌋` (both slices), `≥ k` (even slice), as an **arithmetic statement
> about `v₂(c!·k!) − v₂(den)` plus a *lower* bound `v₂(M_{b+i}) ≥ ⌈floor − (v₂(c!k!)−v₂(den))⌉`**.
> Do **not** target the exact content (provably c-dependent, no low-modulus closed form — last
> cycle's `c=5` vs `c=9`), and do **not** rely on a per-term Kummer bound (cancellation-born content,
> `c=7,k=5` even: `v₂(M)=4` lifts `+3`). The `μ<0` structure shows `den` is the heavy 2-content
> carrier at deep `k`; the Lemma is really "`den`'s 2-content ≤ `v₂(c!k!) + v₂(M) − floor`".

## Route-1 q-lift stub — status: **insufficient at `q=-1` alone, browse needed** (`jobB_route1_qlift_stub.py`)
A clean *negative-with-reason*:
- I first guessed `ord_{q=-1}[n,k]_q = v₂(C(n,k))` — **FALSE** (46 mismatches). Correct fact (0
  mismatch): `ord_{q=-1}[n,k]_q = ⌊n/2⌋−⌊k/2⌋−⌊(n−k)/2⌋ ∈ {0,1}`. A single Gaussian binomial's
  `Φ₂`-multiplicity is **≤ 1**.
- `ord_{q=-1}[m]_q! = ⌊m/2⌋` (verified) — only the **leading** binary layer `⌊m/2^1⌋` of
  `v₂(m!)=Σ_r⌊m/2^r⌋`. So `q=-1` evaluation **cannot** reproduce contents of size 5–8; it sees the
  `r=1` layer (≈ the floor `2⌊k/2⌋`) only.
- **What the route needs (browse `2502.06032`):** the Gatzweiler–Krattenthaler **q-determinant**
  (LGV lift of the 3×3 Jacobi–Trudi minor, so no alternating cancellation hides content) **and** the
  full cyclotomic tower `Φ_{2^r}` at `2^r`-th roots of unity (not just `Φ₂`). The surplus
  `content − 2⌊k/2⌋` the table shows is precisely the `r≥2` layers. Machinery (`qbinom`, `phi2_mult`,
  `qfact`) is in place for the post-browse plug-in.
