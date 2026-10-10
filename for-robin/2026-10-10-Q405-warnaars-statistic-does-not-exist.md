# Warnaar's "additional statistic on cylindric tableaux" does not exist

**2026-10-10, PROVE session, Q405.** Robin — this one came out negative and I
think the negative is the better result, so I want you to look at it.

**Paper:** `proofs/2026-10-10-Q405-twisted-trace-cylindric-statistic.tex`
(9 pages, compiles clean: `errors=0 undefined_cs=0` from a counter I planted a
bad macro against to confirm it moves).
**Registry:** `proofs/registry/cylindric-statistic.json` — validates OK, and the
bogus-path control moved the count 1 → 2 before I fixed the one real violation.

## What I set out to do and what actually happened

The brief had me compute a twisted trace $Z = \mathrm{Tr}(D\,T^N)$ of a cylindric
transfer matrix with a diagonal twist, and ask whether the grading it induces is
the statistic Warnaar asks for in `2511.17034` §6. I expected the answer to be
"no, because a diagonal twist is *linear* in the winding number while Warnaar's
exponent $t^{M y_i^2 - j y_i}$ is *quadratic*." That expectation was right about
the arithmetic and wrong about where the obstruction lives, and I never had to
build the vertex model.

## The result

Warnaar asks (verbatim, last sentence of his HKKO paragraph): *"Is it possible to
extend this result to (Eq_generalise) by introducing an additional statistic on
cylindric tableaux?"* — "this result" being HKKO `2301.13117` Thm 3.3, which
evaluates the $t=1$ case as a signed count of cylindric tableaux.

**Theorem.** For every $\ell \ge 0$, put $M = \ell+2$ and take $k=1$,
$w = 2M-2$, $n = N = 2M$, content $\alpha = (1^n)$. Then
$$ c_\alpha(t) = C_M - (2M-1)\,t^{M-1} + t^{M+1}, \qquad A_\alpha = C_M, \qquad B_\alpha = 2M-2, $$
where $C_M$ is Catalan, $A_\alpha$ counts the cylindric tableaux of that content
carrying HKKO's $+1$ and $B_\alpha$ those carrying $-1$. A monomial weight
$t^{\mathrm{stat}(T)}$ can only *redistribute* tableaux among powers of $t$, so
the coefficient $-(2M-1)$ needs $2M-1$ tableaux on the minus side and there are
only $2M-2$. **It fails by exactly one tableau, at every $\ell$.**

At $\ell = 0$ it fails harder: there the coefficient is $2-3t+t^3$, whose
coefficients have absolute values summing to $6$, while the total number of
$(2,2)$-cylindric tableaux of content $(1^4)$ — over *all* shapes in
$\mathrm{Par}(2,2)$, not just the two HKKO uses — is $4$. So at $\ell=0$ you
cannot even rescue it by re-choosing the signs.

$A_\alpha = C_M$ because the tableaux are length-$2M$ Dyck paths and the ceiling
$w = 2M-2$ never bites; $B_\alpha = 2M-2$ because those paths are all-up-except-
one-down and the down-step can sit in exactly $2M-2$ of the $2M$ positions. The
whole proof fits on a page and every number in it is derived, not computed.

## Why this kills the twist ansatz, and not for the reason I expected

A twisted trace attaches **one monomial in $t$ per configuration**. So it is
subject to the same count, and it dies on *cardinality* — not on the twist being
diagonal, not on linear-versus-quadratic. The thing I went in to measure turned
out to be downstream of the obstruction. I kept the linear-versus-quadratic
observation in the paper as Prop. 5.2 but flagged it as conditional (it needs $y$
to be the winding sector, which I did not prove) and noted that the corollary
does not need it.

## What I think the answer actually is — and this is the bit I want your eyes on

The deformation must attach a **polynomial** $\pi_T(t)$ to each cylindric
tableau, not a power of $t$. That is the shape of the classical answer one level
down: Macdonald III (5.11′) gives $P_{\lambda/\mu} = \sum_T \psi_T(t)\,x^T$ with
$\psi_T$ a product of factors $(1-t^{m_j(\mu)})$. And at $\ell=0$ Warnaar's
determinant *is* $\sum_{\lambda\text{ even},\lambda_1\le 2k}P_\lambda(x;t)$ — his
Conjecture (Eq_PCn), which I verified numerically at $(n,k)=(3,1),(4,1),(3,2)$
against an independent Hall–Littlewood engine. So $t$ here is a Hall–Littlewood
**bulk** parameter, which deforms a tableau count by products of $1-t^m$; it is
not a boundary fugacity, which would deform it by a power. The quadratic Gaussian
is what such a bulk weight *sums to* after the affine folding, not what it *is*
on the tableau side.

I have opened that as **Q410** with the normalisation trap written down
(`memory/questions/Q410-polynomial-weight-on-cylindric-tableaux.md`).

## What I am NOT claiming

- For $\ell \ge 1$ the *re-signed* form is **not** excluded. The capacity test is
  necessary, not sufficient, and at $\ell \ge 1$ it is simply not violated — the
  slack comes entirely from the cylindric shapes with $c^-_{2k,w}(\lambda) = 0$,
  which HKKO's theorem never mentions. At $(1,1,6)$ the shape $(4,2)$ supplies 9
  of the 18 tableaux. An empty violation column is an absence of evidence.
- I have not constructed $\pi_T(t)$.
- Prop. 5.2 is conditional, as above.

## The two questions for you

1. **Am I reading the question too literally?** Everything turns on "statistic"
   meaning $t^{\mathrm{stat}(T)}$ rather than a polynomial weight. I think that
   *is* the standard reading — it is what "statistic" means in Warnaar's own
   $q,t$-Rogers–Ramanujan identities, and his $t^{s^2}$ in the $k=\ell=1$ example
   is literally a power — but you know this literature's idiom better than I do.
   If he means it loosely, my theorem is a precise answer to a question he did
   not ask, and should be framed as "the monomial form is impossible, so here is
   the constraint on the polynomial form."
2. **Is the $\ell=0$ case in scope?** Eq_generalise is stated for $\ell$ a
   nonnegative integer and at $\ell=0$ it is his own Conjecture (Eq_PCn), so I
   believe yes. But $\ell=0$ is where the strong version of the obstruction lives,
   so if he intends $\ell \ge 1$ the headline weakens from "no signed monomial
   weight at all" to "not with HKKO's signs."

## Provenance, stated plainly

The *idea* of looking for the statistic on the integrable side came from a
Mittag-Leffler talk I did not have a written version of. It is **not** an input
to anything: no construction is attributed to it, and nothing in the paper cites
it. Neither Warnaar nor HKKO contains the words "transfer matrix" or "vertex
model"; the integrable framing was mine, and in the end it was not needed. I say
this in §1 of the paper so no reader mistakes the provenance.

## Verification trail (all controls, including the one that was vacuous)

- $\dot e_{n-r} = (x_1\cdots x_n)^{-1} f_r$: verified, all $r$, $n=3$.
- HKKO Thm 3.3 as I read it: exact at $(h,w,n) = (1,1,4), (1,2,3), (1,2,4),
  (1,4,4), (1,4,5), (2,2,4)$. **Three wrong $c^-$ variants** ($c^+$, unsigned,
  wall at $w-1$) **all fail** at three parameter points — the gate discriminates.
- My transcription of Eq_generalise including the exponent: pinned by matching
  Eq_PCn against the independent Hall–Littlewood engine. **Five of six exponent
  perturbations fail.** The sixth is **vacuous** and I say so in the paper: the
  gate cannot tell $-j y_i$ (column) from $-i y_i$ (row); both pass at $k=1,2$.
  Nothing in the result depends on that distinction.
- Tableau counts computed **twice** by independent mechanisms (cylindric
  Jacobi–Trudi determinant; HKKO's lattice-path model) — 476 contents compared,
  0 mismatches.
- Witness formula verified symbolically at $\ell = 0,1,2,3$: $2-3t+t^3$,
  $5-5t^2+t^4$, $14-7t^3+t^5$, $42-9t^4+t^6$ with $(A,B) = (2,2), (5,4), (14,6),
  (42,8)$ — matching $(C_M, 2M-2)$ in all four.
- Exhaustive per-monomial sweep, 15 parameter points, ~59k monomials, with
  $c_\alpha(1) = A_\alpha - B_\alpha$ asserted on **every** monomial.
- Both sources upgraded `deep-read` → `verified-quote` with file+line+md5
  locators. My edit added **0** problems to the sources index (baseline 60,
  after 60, measured in separate calls).
- Orevkov–Orevkov $n \le 14$ trap: recorded **not applicable**, nothing here is a
  polytope computation.
