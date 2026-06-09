# SUMMARY — two-row d=4 law, b ≡ 2,3 (mod 4): code session 2026-06-09

Four jobs, all run with exact arithmetic to the **b ≤ 119 frontier** (60 values
b ≡ 2,3 mod 4). Reusable engine `_qbcore.py`; polynomials cached in
`results/Qb_cache.pkl`; per-job data in `results/*.csv`.

## The one-line headline

**Q_b is irreducible over ℚ for every b ≡ 2,3 (mod 4), 6 ≤ b ≤ 119** (Job 4,
certified). Since deg Q_b ≥ 3 there, irreducibility ⟹ no rational root ⟹ the
**two-row d=4 law holds for all those b** — no analytic estimate required. The
remaining task is a *uniform proof of irreducibility*, and Jobs 1–3 say which
route can deliver it.

## Per-job results

- **Job 1 (discriminant / Newton — the go/no-go).**
  - Coefficient Newton polygon: **no pure Eisenstein/Dumas prime for any b ≥ 6**
    → the classical recurrence-Newton **irreducibility** method is DEAD here
    (parameter-vs-degree mismatch confirmed).
  - Discriminant: **Hajir's lever (prime > deg, odd multiplicity) is PRESENT**
    for all b ≥ 6 (disc non-square ⇒ Galois ⊄ A_d ⇒ contains a transposition).
  - Clean 2-adic fact: **v₂(Q_b(0)) = 1** for all 58 values b ≥ 6.
  - Cross-resultants: gcd = 1 for all consecutive pairs (structured-coprime).
  - **Verdict: switch to the Galois-group route (Dedekind cycle types ⟹ S_d).**

- **Job 2 (central-trinomial identity).** Proved the exact map, symbolic in
  ℤ[m] for b ≤ 50:  `1+su+u² = (1+u+u²)+iu` gives
  **I_b(m) = Σ_{k odd} (−1)^{(k−1)/2} C(m,k)[τ(m−k,b−k) − τ(m−k,b−1−k)]**,
  τ = trinomial triangle A027907, with **A002426 = τ(n,n)** entering on the
  diagonal m = b. This is the cleanest closed form for I_b(m).

- **Job 3 (finite window).** **No integer root of Q_b in [b, ⌊0.33b²⌋]** for any
  b ≤ 119 (modular sieve, 0 survivors). Frontier extended past b ≤ 70/150.

- **Job 4 (irreducibility census).** All 60 Q_b irreducible over ℚ (mod-p
  certificates); disc square only for trivial b=2,3. Prop 1 verified
  symbolically: **g_b(m) = C_b^{(−m)}(−(1+i)/2)** (Gegenbauer in the parameter).

## Recommended next PROVE step

Prove irreducibility via a **Frobenius d-cycle**: for every b ≡ 2,3, exhibit a
prime p with **Q_b irreducible mod p** — then the Frobenius at p is a d-cycle,
Gal(Q_b) is transitive, Q_b is irreducible over ℚ, and (deg ≥ 3) the law follows.
(Holds for 9/60 b with a single small prime; by Chebotarev a d-cycle prime has
density 1/d when Gal = S_d, so one always exists — the task is a uniform
Artin–Schreier/reduction argument producing it for every b.) The non-square
discriminant + odd-multiplicity prime (Job 1) then upgrade Gal(Q_b) to the full
**S_{⌊b/2⌋}** via Jordan, but that transposition is a refinement, not needed for
the law. Irreducibility is *stronger* than (★) and closes the last open residue
class b ≡ 2,3 outright (b ≡ 0,1 already proved 2-adic).
```
files: _qbcore.py, build_cache.py, job1_discriminant.py, job1b_resultants.py,
       job1c_newton.py, job2_trinomial.py, job3_window_roots.py,
       job4_irreducibility.py, job4b_gegenbauer.py
data:  results/Qb_cache.pkl, results/job1_disc.csv, results/job1_resultants.csv,
       results/job3_window.csv, results/job4_irred.csv, results/*_run.txt
```
