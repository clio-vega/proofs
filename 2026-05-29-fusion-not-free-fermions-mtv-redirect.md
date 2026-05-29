# The order law's home: not free fermions, but maybe the Wronski map (2026-05-29)

Robin —

A clean negative and a promising redirect today, on the question of *where* the graded order law
$\mathrm{ord}_{x=q^2}Z_\lambda = \tau(\tau+1)/2 = \binom{\tau+1}{2}$ lives.

**The setup, briefly.** $Z_\lambda(x) = \mathrm{tr}_{V^\lambda} M(x,\dots,x)$, where $M$ is the
uniform-rapidity monodromy built from the Baxterised Hecke R-matrix $\check R(x) = P_+ + r(x)P_-$,
$r(x)=(q^2-x)/(q(x-1))$. Proved for hooks (chip-firing), verified on 13 shapes.

**The negative (decisive).** Last night's dream had me chasing the Bump–Hardt–Scrimshaw "Factorial
Fock free fermions" paper — their row transfer matrix *is* a half vertex operator, and its matrix
elements are KP tau functions, so I hoped $\binom{\tau+1}{2}$ would be a fermion-mode count (the order
to which a Slater determinant vanishes when $\tau+1$ modes collide — a Vandermonde). It's a beautiful
picture and it's **wrong for my object**: I computed the free-fermion defect of my R-matrix,
$$a_1a_2 + b_1b_2 - c_1c_2 = \frac{(q-1)(q+x)}{q(x-1)} \neq 0,$$
and this defect is *gauge-invariant*. My $\check R$ is a genuinely interacting six-vertex R-matrix
($\Delta\neq0$), so it cannot embed in the BHS/Naprienko free-fermionic lattice model. The
Slater/Vandermonde collision lives in *free*-fermion Fock space; my world is interacting. So that
specific route is closed — with a precise reason, not a vague mismatch.

**What $x=q^2$ actually is.** $\check R(q^2)=P_+$ — it projects onto the symmetric channel. The
degenerate point is a **fusion** point; the uniform-rapidity staircase fuses all $n$ sites onto the
symmetric line, and $\tau+1 = 2\lambda'_1 - n$ is the "antisymmetric excess" that can't fit. That's a
clean integrable statement.

**The redirect (promising, but I'm holding it at arm's length).** The triangular number $\binom{k}{2}$
*is* an exact theorem in the **Mukhin–Tarasov–Varchenko Wronski-map** framework — the Wronskian
vanishes to order $|\lambda(\text{ramification})|$, and the minimal collision of $k$ exponents off the
generic staircase gives exactly $\binom{k}{2}$. It has the right XXZ/trigonometric avatar
(Tarasov–Varchenko quasi-polynomials, arXiv 1210.2315), it's indexed by a Schur–Weyl multiplicity of
exactly my type, and it predicts what I see numerically — only the *minimum* branch is the clean
triangular number; the full exponent spectrum is the specific ramification sequence, not a lattice.
And — this is the part I like — the Wronski map *is* Schubert calculus on the Grassmannian, which is the
puzzle path of our own seed. If this is the home, the order law comes home to where we started.

**The honest catch.** The staircase mechanism as stated computes the *total* Wronskian order
$=\mathrm{ord}\det M$ ($12, 30, 60,\dots$ on my shapes) — not my *minimum* branch $\binom{\tau+1}{2}$.
The reconciliation has to be a small $(\tau+1)\times(\tau+1)$ confluent discrete Wronskian of just the
colliding modes. Building that and checking its order $=\binom{\tau+1}{2}$ is the decisive test, and I
haven't run it yet. Until then MTV is a candidate, not a home — I've been burned by ten matched-but-distinct
"shadows" and I'm trying not to make it eleven by enthusiasm.

The main proof (internal per-subset operator vanishing) is unchanged and remains the primary route; MTV
is a parallel lead to try if the internal machinery stalls.

Scripts and the hypothesis note are in `scratch/2026-05-29-fermion-collision/`; branch-exponent data
reused from `scratch/2026-05-23-omega-monodromy-verify/`.

— Clio
