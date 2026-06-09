# Two-row d=4 law, b ≡ 2,3 (mod 4) — code session 2026-06-09

Arithmetic / family-structure attack on **(★)**: *Q_b(m) has no integer root
m ≥ b for b ≡ 2,3 (mod 4)*, equivalent to the two-row d=4 order law for the
last open residue classes (b ≡ 0,1 already proved 2-adically).

`Q_b(m)` is the degree-⌊b/2⌋ cofactor of `I_b(m) = Im[u^b]((1−u)(1+su+u²)^m)`,
s = 1+i, after dividing out the forced roots {0,…,⌊(b−1)/2⌋}. Monic over ℤ.
All results to the **b ≤ 119 frontier** (60 values), exact arithmetic.

## Headline

**Q_b is irreducible over ℚ for every b ≡ 2,3 (mod 4), 6 ≤ b ≤ 119.** Since
deg Q_b ≥ 3 there, irreducibility ⟹ no rational root ⟹ the law holds for all
those b. The open problem reduces to a *uniform proof of irreducibility*.

## Files

| script | job | what it does |
|--------|-----|--------------|
| `_qbcore.py` | core | builds primitive integer Q_b (engine + cache loader) |
| `build_cache.py` | core | pickles all Q_b, b ≤ 120 → `results/Qb_cache.pkl` |
| `job1_discriminant.py` | 1A | disc(Q_b), Hajir odd-multiplicity-prime lever |
| `job1b_resultants.py` | 1C | cross-resultants Res(Q_b,Q_{b'}), gcd structure |
| `job1c_newton.py` | 1B | coefficient Newton polygons (Dumas/Eisenstein scan) |
| `job2_trinomial.py` | 2 | exact central-trinomial identity (TRI), A002426 |
| `job3_window_roots.py` | 3 | no integer root in [b, ⌊0.33b²⌋] (modular sieve) |
| `job4_irreducibility.py` | 4 | irreducibility census (mod-p), disc-square test |
| `job4b_gegenbauer.py` | 4 | Prop 1: g_b(m)=C_b^{(−m)}(−(1+i)/2), symbolic |

FINDINGS-*.md hold the per-job write-ups; **FINDINGS-2026-06-09-SUMMARY.md** is
the overview. `results/` holds the CSV data, run logs, and the polynomial cache.

## Verdict for PROVE

- **Coefficient Newton-polygon irreducibility (Hajir/Filaseta) is DEAD** for
  b ≥ 6: no pure Eisenstein/Dumas prime exists (job1c). Parameter-vs-degree
  mismatch, as feared.
- **Discriminant carries the Hajir lever** (prime > deg, odd multiplicity) and is
  **non-square** for all b ≥ 6 → Gal(Q_b) contains a transposition.
- **Recommended route:** prove irreducibility via a **Frobenius d-cycle** — for
  each b exhibit a prime p with Q_b irreducible mod p (holds for 9/60 with a
  small prime; Chebotarev guarantees existence when Gal = S_d). The non-square
  disc then upgrades Gal to S_{⌊b/2⌋} (Jordan), a refinement beyond the law.
- **Job 2 bonus:** `1+su+u² = (1+u+u²)+iu` ⟹ exact identity
  `I_b(m) = Σ_{k odd}(−1)^{(k−1)/2}C(m,k)[τ(m−k,b−k)−τ(m−k,b−1−k)]`
  (τ = trinomial triangle A027907; A002426 = τ(n,n) on the diagonal m=b),
  verified symbolically in ℤ[m] for b ≤ 50.
