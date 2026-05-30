# Two prunings, 2026-05-31 — the Bethe home is closed, and the order-law proof is split-independent

Hi Robin — short update on the order law `ord_{x=q²} Z_λ = τ(τ+1)/2` for off-hook
λ=(2,2,1^m). Today closed two loose threads, one external (where does this live?) and
one internal (the remaining gap in the proof).

## 1. The Bethe-subalgebra home is a shadow (decisive)

Last week's most promising lead was that my Baxterised factor `Ř(x)` equals an
Isaev–Kirillov Baxterised Hecke generator (exact, per-factor). The hope was that the
assembled monodromy `M(x)` (and hence `Z_λ = tr_{V^λ} M`) would be a Bethe-subalgebra
element, making the order law a spectral invariant of a named commutative algebra.

**It is not.** Exact computation (q=2,3) shows the family `{M_λ(x)}` is **non-commutative**
— `[M(x),M(y)] ≠ 0` for every shape I tested. A non-commuting family generates no
commutative subalgebra, so it cannot be a Bethe/Jucys–Murphy family. The per-factor
identity is real, but the assembled object is the q-KZ / full-twist (w₀) monodromy, not a
commuting transfer matrix. A web sweep confirmed there is no known formula expressing a
single irreducible character χ^λ via Bethe spectral data, and the vanishing-order-at-fusion
invariant appears genuinely novel — there's no off-the-shelf machine to borrow.

The silver lining is sharp: the order law is fundamentally **spectral** — with {d_j} the
orders to which the eigenvalues of M(x) collide at the fusion point x=q²,
`ord Z_λ = min_j d_j = τ(τ+1)/2` and `ord det M = Σ_j d_j`. Plus a clean new fact:
**M(q²) = 0 on V^λ exactly when τ ≥ 1** (rank 1 when τ=0). The right neighbourhood is
q-difference / q-KZ monodromy exponents, not Bethe/Gaudin.

## 2. The full order law's residual is resolved (conceptually)

The proof had a gap I'd called "placement not count": for the intermediate split
1 ≤ |S_T| < τ, my per-split count bound came out *below* D = C(τ+1,2), yet computationally
no subset ever survives below D. The cause turned out to be a precise conflation:

> My proved Deficit Laws are about factor **deletion** (delete u factors ⇒ E⁻-depth drops
> by u). The order-law word **inserts** P_{-1}, and *insertion is not deletion* — a tail
> P_{-1} generically does **not** lower the tail depth.

So the effective depth to consume is the full τ regardless of split; a depth-τ staircase
costs exactly D insertions; a tail P_{-1} merely *pre-pays* one of those D. Net:
**min |S| = D, split-independent** (verified exact for m=3,4; survivors at |S|=D: 11 and 56).
This makes the next prove session a clean target: a single depth-τ staircase-consumption
lower bound that should close the full off-hook order law (upgrading the proved `ord ≥ τ`).

Nothing here needs action from you — just flagging where the program stands. Scripts are in
`scratch/2026-05-31-{bethe-spectral,intermediate-split}/`. (Gmail's still locked on my end —
~14 sessions — so this note is the channel; a `/mcp` re-auth would reopen email when you get
a chance.)

— Clio
