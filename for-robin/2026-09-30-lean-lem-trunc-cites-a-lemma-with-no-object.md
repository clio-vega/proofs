# The three lemmas under the m=2 theorem are Lean-verified — and one citation inside them points at a lemma that has no object on the case it is cited for

2026-09-30 c2, LEAN session. Sorry count **0**.

## What is verified

`TworowD4Kernel/DiscreteConcavity.lean` in `clio-vega/tworow-d4-kernel`, commit `a991128`
(local `==` `git ls-remote`, both printed), imported from the root `TworowD4Kernel.lean` so it is
inside the axiom-audit closure. `lake build` exit 0, 3176 jobs, module compile time 39 s.

https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/DiscreteConcavity.lean

Everything in `proofs/2026-09-30-c1-cylindric-kostka-logconcavity.tex` — which proves condition
(A) at `m = 2` for 79.9% of slices — stands on three small results, and all three were proved
there by prose case analysis with because-clauses. They are now type-checked:

* `PFtwo_of_concave_on_interval_support` — `lem:conc` (l.106),
* `PFtwo_posPart` — `lem:trunc` (l.119),
* `GR_PFtwo` — `prop:regII` (l.347).

All 22 audited declarations return `[propext, Classical.choice, Quot.sound]` or less.
`PF₂` is stated as the paper writes it — nonnegative, interval support, log-concave — as three
separate structure fields, not transported through a Mathlib predicate, so the name `PFtwo` is
never load-bearing.

## The finding worth your time

`lem:trunc`'s proof says:

> "at the two ends of the support the argument of Lemma `lem:conc` applies verbatim."

`lem:conc` is stated for support a **bounded** interval `J = [j₀,j₁]`. `max(c,0)` need not have
bounded support: `c ≡ 1` is concave and the support of its positive part is all of `ℤ`. So on that
case `lem:conc` has no "two ends" to offer — the phrase has no referent. The Lean witness is
`posPart_support_unbounded`.

The **conclusion** of `lem:trunc` is true, and `PFtwo_posPart` proves it. But it proves it
*without* invoking `lem:conc`, through a quasiconcavity lemma — `min(c r, c t) ≤ c s` for
`r ≤ s ≤ t`, by strong induction on the span `t − r`, on top of the fact that the increment
`s ↦ c(s+1) − c(s)` is nonincreasing. That is the clause the paper leaves ungraded
("`{c>0}` is an interval **because** `c` is concave"), and it is the honest route.

So the defect is in the **pointer**, not in the mathematics. It is the reason I keep pointing Lean
at because-clauses rather than at theorems I trust: every instrument I own grades a *proposition*,
and a citation inside a proof is a pointer to one. The sweep that validated `lem:trunc`'s
conclusion on thousands of slices could never have seen this, because the conclusion is true.

Two clean repairs, both one paragraph, and the choice between them is editorial rather than
mathematical — which is why I did not make it in a LEAN session:

1. restate `lem:conc` for an interval `J ⊆ ℤ` not necessarily bounded (then the citation is valid
   as written), or
2. delete the citation from `lem:trunc` and give the quasiconcavity argument inline — two lines,
   and it is now Lean-checked.

I lean toward (1), because `prop:regI` cites `lem:trunc` for a `w_I` that *is* boundedly
supported, so (1) keeps one lemma serving both consumers instead of splitting them.

## A second, smaller one

`prop:regII` closes with "a nondecreasing concave function composed with a concave function is
concave, so `G_R = ρ_R ∘ H` is concave on `[ℓ,h]`". True, and **not needed**. On `ℓ < s < h` all
three of `H(s−1), H(s), H(s+1)` are `≥ 1`, so the positive part is inactive on every summand, and
what is left per summand is `min(H(s−1),k) + min(H(s+1),k) ≤ 2 min(H(s),k)` — linear arithmetic on
`min`s of affine functions, one `omega` each under `Finset.sum_le_sum`. Five lines instead of a
general composition lemma. The factorisation the paper names is still recorded, as
`GR_eq_rho_comp_H`, so the exposition loses nothing; I would just replace the closing sentence
with the certificate.

## What is deliberately not verified

`thm:regions` is **not** Lean-verified and I did not promote its registry node. `prop:regI`
(regions I and III) routes through closure of `PF₂` under convolution, and convolution on `ℤ` is
not defined in the file — `G_R` is *defined* by the trapezoid min-formula rather than derived from
it. Two of the four regions are verified, not four. The `∓∞`-extended version of `lem:trunc`
("the same proof applies") is a separate statement and is not proved.

## For the Lemma T attempt

The sharpest thing in the file for the parallel PROVE session is `Gsharp_not_logConcave`: the
`rem:sharp` witness as a kernel computation. `Φ = (−1,0,1)` is concave and 1-Lipschitz,
`Ψ = (−3,0,4)` is convex and **not** 1-Lipschitz (its step from 1 to 2 is 4), and
`G(s) = Σ_{i+j=s}(Φ(i) − Ψ(j) + 1)₊` is not log-concave — `G = (3,4,6,2,0)`, computed by `norm_num`
over the Finset sum rather than copied from the paper, with `4·4 = 16 < 18 = 3·6` at `s = 1`. The
1-Lipschitz hypothesis is load-bearing and the failure is attributable to `Ψ` alone. Any attempt
at Lemma T that does not spend that hypothesis somewhere is proving something false.

And `bump_PFtwo` together with `bump_not_concave` records that `lem:conc` is strictly one-way:
`(1,3,6,7,6,3,1)` is `PF₂` and is not concave. `prop:regII` proves `PF₂` by proving *concavity on
the support*, which is more than it needs — and that is precisely why the route cannot be expected
to extend past `m = 2`.
