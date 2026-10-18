# Pre-registered predictions, PROVE 2026-10-18 c1 — written BEFORE running `sc_census.py`

Theorem (derived this session, from rho^A = 1 + the colength count):

> For any chi such that chi(w_m) is a UNIT of O_ij (regular and non-vanishing at the
> generic point of {w_i = w_j}), every structure constant of the renormalised basis
> Psi_chi * Atilde is regular on every A-type divisor {w_i = w_j}.

Every chi in the eight-arm family is such a unit:
`id`=1, `psi`, `h1`, `sq`, `unrm`, `xpsi`, `x2psi`, `x` — each is a polynomial in w
(or a monomial) whose vanishing locus is a q-deformed or monomial divisor, never
{w_i = w_j}.

## P1 — the hard prediction
At N=3 type C, **ZERO** of the non-Laurent structure constants in **any** of the eight
arms has an A-type factor in its denominator.  Every non-Laurent denominator is
`pair`, `bdy`, or `qdef`.

This is falsifiable and the falsifier is cheap: `x`, `h1`, `sq`, `unrm` each produce
57/115 non-Laurent constants.  If even ONE of those 57 carries a (w_a - w_b) factor,
the theorem above is FALSE.

## P2 — the baseline counts must reproduce
`id` -> 13/115 non-Laurent, all pair.  `psi` -> 0/115.  (These reproduce the
2026-10-17 c2 arm logs; a disagreement means the engine or the solve has drifted.)

## P3 — what this does NOT predict
It says nothing about which arms are clean at pair/boundary divisors — that is
Theorem C and is already settled.  The A-type direction is INSENSITIVE to the whole
family, which is the exact complement of Theorem C: (COB) and (SR) are the only
conditions, and A-type adds none.  So P1 passing is NOT evidence for psi over x.

## P4 — the control that must FIRE
A plant with an A-type factor must be detected by the classifier+census pipeline, or
P1's green is decoration.  Arm `atp`: chi(v) = 1 is replaced by a basis lesion that
divides A_{M,I} by (w_1 - w_2) for |M| odd.  That is NOT a point-product
renormalisation — it is a deliberate violation of the unit hypothesis — and it must
produce A-type poles.
