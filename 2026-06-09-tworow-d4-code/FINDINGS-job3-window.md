# FINDINGS — Job 3: finite-window no-integer-root check

**Date:** 2026-06-09 (code session)
**Script:** `job3_window_roots.py`
**Data:** `results/job3_window.csv`

## Result

For every b ≡ 2,3 (mod 4), 2 ≤ b ≤ 119:

> **Q_b has NO integer root in the window [b, ⌊0.33 b²⌋].   ALL CLEAR (60/60).**

The largest window checked is b=119: [119, 4673]. This is the finite check the
no-integer-root proof leans on, now extended to the **b ≤ 120 frontier**
(Theorem B's previous frontier was b ≤ 70/150 depending on metric).

## Method

Q_b is monic, so any rational root is an integer dividing Q_b(0). We do **not**
factor the (hundreds-to-thousands-of-digit) constant term. Instead we
modular-sieve the window against four primes p ≈ 10⁶: a genuine integer root m
in [b, ⌊0.33b²⌋] must satisfy Q_b(m) ≡ 0 (mod p) for every p, so we intersect
the per-prime vanishing residues over the window. The intersection was **empty
for every b** (no survivor even reached the exact-check stage).

## Scope note (honest)

The window upper bound ⌊0.33 b²⌋ is imported from the analytic Prop 3 of the
PROVE track (all real roots of I_b lie in [b, 0.33 b²]). The integer-root check
is therefore complete **conditional on Prop 3's root-localisation**. Independently
of Prop 3, Job 4 certifies Q_b irreducible (deg ≥ 3 for b ≥ 6), which already
forbids *any* rational root — so the no-root conclusion holds unconditionally in
range, and this window scan is a cross-check that found zero violations.
