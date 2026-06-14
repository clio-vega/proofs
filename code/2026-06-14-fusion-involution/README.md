# 2026-06-14 code session — level-2 fusion probe (Job A) + top-generator involution (Job B)

Two jobs from `state/CODE.md`. Both reach a plain verdict; every script carries a cross-check.

## Job A — is `ord_{q=−1}G_λ = ⌊|2-core(λ)|/2⌋` a level-2 fusion multiplicity? (Dobner cylindric lead)

**Verdict: DECISIVE MISMATCH — prune the "level-2 fusion multiplicity" form.**

The order law takes unbounded values `0,1,3,5,7,…`, occurring exactly on wide shapes `λ_1 ≥ 3`. But a
level-2 weight of `sl_n` has `λ_1 ≤ 2`, and on that entire domain `⌊|2-core|/2⌋ ∈ {0,1}`. So the
function is not in the image of any level-2 fusion bracket. First gap witnesses: `(3,2,1)` (value 3),
`(4,3,2,1)` (value 5). The `d→level-d` test (Part 2) confirms `ord_{ζ_d}` is `d=2`-only (rich at d=2,
collapses to `{0,1}` for d≥3) and does not track `⌊|d-core|/d⌋` for d≥3 either. Exact carrier stays the
**2-core** (rank-2 affine ŝl₂ / Littlewood `t=2`), not level-2 Verlinde fusion.

- `fusion_lib.py` — LR coefficients + level-k Kac-Walton fusion for `sl_n`. **Validated** against the
  Ising (sl₂ level 2) rules, the sl₂ level-k closed Verlinde formula (139/139), and sl₃ level-1 ℤ/3.
  Run `python3 fusion_lib.py`.
- `2026-06-14-jobA-fusion.py` — order-law data (formula vs `ord_{q=−1}Σ_T q^{s(T)}`, 0 mismatches),
  the d→level-d test, and the fusion search. See `FINDINGS-2026-06-14-jobA-fusion.md`.

## Job B — even-`|J*|` fpf involution = the TOP-generator toggle, paired by residue mod π²

**Verdict: CONFIRMED on every tie (1043 shapes, m ≤ 16). Sawin adjacency FALSIFIED.**

The fpf involution forcing `|J*|` even is XOR-with-the-top-generator of the J\* box (for the `|J*|=4`
family `{3,5,7,9}` this is `j↔j+4`). It pairs **residue-equal** leading units mod π² (`v_π(sum)≥2`),
while Sawin's adjacent `j↔j+2` pairs residue-opposite (`v_π(sum)=1`). Underlying law (0 failures):
`residue(u_j) = residue(u_{j'}) ⟺ j ≡ j' (mod 4)`.

- `2026-06-14-jobB-residue-involution.py` — residue tables, the involution check, Sawin comparison.
  See `FINDINGS-2026-06-14-jobB-residue-involution.md`.

## Bundled support libraries (copied for standalone reproducibility)

`job_jstar_engine.py` (M_j chains, J\*), `syt_lib.py` (SYT, `s(T)`, `G_λ(q)`, 2-core/quotient),
`dcore_lib.py` (d-core/quotient via abacus, `ord_at_root`). Scripts prepend their own directory to
`sys.path`, so they run from a fresh checkout: `python3 2026-06-14-jobA-fusion.py`.
