# threerow-boundary — Boundary Lemma (Gap 3) verification

Verifies the Boundary Lemma for the three-row d=4 family (proof:
`projects/proofs/2026-06-15-boundary-lemma-threerow.md`).

- `mn.py`          — Murnaghan–Nakayama M_j (shared, from threerow-c1)
- `boundary_forms.py` — fits boundary M_{b+i} closed forms (Lemmas T, H), c=1..5
- `general_top.py` — verifies hook formula (H), top M_{b+c} (T), clean Δ(b+c), c=1..5
- `c1_clean.py` / `c1_bdry.py` / `c1_fullproof.py` — c=1 Prop 1 + boundary lemma conclusion
- `c2_check.py` / `c2_final.py` / `c2_prod.py` — c=2 Δ(b+1),Δ(b+2) identities + full boundary check
- `c3_aeventop.py` / `c3_status.py` — c=3 a-even top identity; full c=3 boundary scan
- `lemma_only.py`  — Lemma F (factor-in-product) pure number-theory check
- `jobA_widesweep.py` — **(2026-06-15 evening)** gated val-gap sweep, c=1,2,3, **m≤200**
  (pure-Python closed forms + MN cross-check): boundary loses by min Δ=2, 0 in-gate ties;
  the only boundary ties are the finite list c=1:b∈{1,2,4}, c=2:b=2, c=3:b=6. Also the
  boundedness probe showing `v₂(M_{b+i})` is NOT bounded (the real engine is Lemma F).
- `jobA_c45.py` — c=4,5 confirmation via the exact 3-row Jacobi–Trudi determinant (min Δ≥8).

### general-`c` master valuation + c=5 (2026-06-17 prove; proof `../../proofs/2026-06-17-generalc-boundary-master-and-c5.md`)
- `theta_scout.py` — interior offset θ∈{0,3} uniform in c (corrects PROVE.md's τ(τ+1)/2)
- `generalc_master.py` — **the master valuation** `v₂(R_i)=v₂(N_i)−v₂(k!)−v₂(a−c+2)−[k≤c−2]v₂(a−b+1)+v₂(Π_i)−E`; verifies `const_i=c!·k!`, vs MN c=4,5,6 (0 mismatch)
- `generalc_content.py` — deficit polys `N_i^{(c)}` and their 2-adic contents `g[c][k]`, c=2..5
- `generalc_Ndiv.py` / `generalc_polydiv.py` — compensation: which deficits are Π-absorbed vs `(a−b+1)|N_i` polynomial divisors
- `generalc_certify.py` — uniform hand bound (master + Lemma P + compensation) ≤ true Δ and > −θ, **c=4..8, m≤110, 0 fail**
- `c5_content.py` — explicit c=5 slice factorisations (the `B(B+1)`-even trick)
- `c5_certify.py` — c=5 boundary hand bound, m≤140, 0 fail
- `c5_fulltheorem.py` — no boundary index is a minimizer, ALL c=5 shapes m<50 (0/4390), |J*|∈{1,2}
- `generalc_subtop_scout.py` — (2026-06-16) top + first-subtop closed forms, F2 c-independence

All checks: 0 mismatches / 0 violations in stated ranges.
See `../FINDINGS-2026-06-15-jobA-boundary.md` for the earlier structural write-up.

### g₀ content floor + alternant tool (2026-06-17 prove; proof `../../proofs/2026-06-17-generalc-g0-content-floor.md`)
- `fast_alt.py`     — fast multinomial extractor for the alternant `M_j=[x^{a+2}y^{b+1}z^c] V E^j H^{2m-2j}` (polynomial continuation, valid for a<b); cross-checks vs MN
- `alternant_check.py` — verifies Lemma 0 (alternant=Mj) and the equal-exponent vanishing
- `oblig_b_large.py` — **Obligation (b) CLOSED**: `M_{b+i}|_{a=b-1}=0` (⟹ (a−b+1)|N_i), c≤15 odd, i≤(c−1)/2, **m≤200, 0 nonzero**
- `oblig_a_proof.py` — **Obligation (a) k≤3**: exact uniform N_{c−k}(a,b,c) + 4-parity substitution, coeff 2-content ≥ 2⌊k/2⌋
- `claimA_verify.py` — Claim A (b-even content =k, b-odd ≥2⌊k/2⌋), c=4..11, k≤4
- `g0_fixeddiv.py`  — full fixed-divisor floor `d(N_i)≥g₀(k)`, master convention, c=4..8 (0 violation)
- `Ni_uniform.py` / `Ni_k23.py` / `Ni_symbolic.py` — uniform closed forms N_{c−k}(a,b,c) per depth

### deep-`k` Claim A exact content census (2026-06-18 code; `../FINDINGS-2026-06-18-jobA-claimA-deep.md`)
- `jobA_deep.py` — exact **fixed-divisor** content `g[c][k]` (min v₂ over slice, NOT
  coeff-gcd), c=4..12, k≤6, both parity slices. Finds true content **exceeds** the old
  `even=k` (that was only the coeff-gcd lower bound).
- `jobA_verify_content.py` — verification gate: Malt-vs-MN (0 bad), grid stability
  64→512, fit-free direct-Malt cross-check (ALL pass).
- `jobA_certify.py` — certified content tables over **b≤200**, `(c mod 4)` grouping;
  **0 violations** of the floor.
- `jobA_law.py` — content indexed by `(c,i)`; H1 floor (`≥2⌊k/2⌋`, even-slice `≥k`,
  0 viol) + H2 **parity-of-`i`** law (odd-i odd-slice = floor for k≤3).
- `jobA_factor.py` — slice factorization `N_i = 2^σ·(a−b+1)·(irreducible block)`;
  no even-spaced linear run, **no `2^k` sibling of Theorem B**.
- `jobA_kummer.py` — pointwise `v₂(N_i)=v₂(M)+v₂(c!k!)−v₂(den)`; content sometimes
  term-by-term, sometimes cancellation-born (+3 at c=7,k=5).
- **Verdict:** exact content c-dependent with no low-modulus closed form (c=5 vs c=9);
  use the certified FLOOR as the Content Lemma, proven via `c!k!/den` arithmetic +
  i-parity, not factorization or pure-Kummer.
