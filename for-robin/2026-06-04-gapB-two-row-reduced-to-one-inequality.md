# Two-row d=4 fibre-vanishing reduced to a single real inequality

The problem: prove `G_λ(i)=0 ⟺ λ=(2,2)`, where `G_λ(i)=⟨s_λ,(h₂+i e₂)^m⟩`, n=2m. Today I cracked the structure of the two-row case ("Gap B") and reduced it to ONE explicit inequality.

Reduction (all verified exactly in ℤ[i] to m=300):
- For two rows, `G_{(a,b)}=r_a−r_{a+1}` with `r_l=[u^l]f(u)^m`, `f(u)=u²+(1+i)u+1`.
- `f` is palindromic, so `r_l=r_{2m−l}`, hence `G_{(a,b)}=r_b−r_{b−1}` (b = smaller part).
- The surprise: set `N_l=|r_l|²`. The sequence `N_0<N_1<…<N_m` is strictly increasing, with the SINGLE exception `N_1=N_2` at m=2 — which is exactly the shape (2,2). Since `|r_l|≠|r_{l−1}|` forces `r_l≠r_{l−1}`, monotonicity of N_l closes the entire two-row case.

This matters because the d=4 fibre value is genuinely complex (the d=3 Schur-positivity proxy is dead at ζ₄), and I had concluded only a (1+i)-adic valuation argument could work. For two rows that was too pessimistic — the modulus N_l is a real, nonnegative, cancellation-free proxy.

What's proved: a 3-term recurrence `(k+1)r_{k+1}=(1+i)(m−k)r_k+(2m−k+1)r_{k−1}`; a Gram system and a 4th-order P-recursion for N_l; monotonicity unconditionally on k ≲ 0.586m, and all but the endpoint elsewhere.

The one remaining gap (the endpoint inequality `D_m=N_m−N_{m−1}>0`, m≥3): the obstruction is the phase of `r_{m−1}r̄_{m−2}`; Cauchy–Schwarz alone is provably insufficient. Details and the proof skeleton are in the two attached tex files. This is my next prove-session target.
