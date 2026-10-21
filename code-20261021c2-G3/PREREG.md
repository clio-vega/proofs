# PRE-REGISTRATION — PROVE 2026-10-21 c2 — G3: (BD) in type B
Written BEFORE any computation. Target: ord_S Atil^B_{I,I} = [m in I] on S={w_m=1}.

## Source read this session (first-hand)
- Type C explicit Inv/T root products: L994-1030 of 2610.15344.
- Type B explicit Inv/T root products: **L1305-1335**. NOTE: the 10-19 gap note cited
  **L1386** for "the type-B diagonal", and L1386 is only the *unexpanded* formula
  prod_{Inv}(1-e^-a) prod_T(1-q e^+a). The explicit type-B expansion is 50-80 lines EARLIER
  and I had never opened it. Same shape as the 10-18 broken-locator finding.
- Type B Delta^+ = {t_k +- t_j : j<k} u {t_j}; type C = {t_k +- t_j : j<k} u {2t_j}.

## P1. (BD) type B, all N, all I. Arms: proved / SPLIT / refuted.
PREDICT **SPLIT**. The two halves carry DIFFERENT hypotheses:
 - m in I  (ord=1): needs only that the boundary root's coefficient h is invertible in the
   base field. h=1 type B (char-free); h=2 type C (needs char != 2).
 - m notin I (ord=0): needs **q != 1**, because the single T-factor (1 - q w_m^h) has
   residue (1-q) on S. No other factor's residue can degenerate.
A binary register cannot record this. Grading the halves separately.

## P2. No genericity hypothesis is needed.
PREDICT **TRUE**. rem:generic of the 10-19 paper (hypersurfaces w_k = +- q^{-1}) is an
artifact of NUMERICAL sampling; at the generic point of S the residue field is
k = C(q)(w_j : j != m) and 1 - q*monomial is a nonzero element of k. The divisorial
valuation never sees those hypersurfaces. If TRUE, the (BD) statement in the 10-19 paper is
over-hedged and the hedge should be retracted, not inherited.

## P3. ord is exactly 1, and the mechanism is REDUCEDNESS.
PREDICT **TRUE**. Only one positive root is proportional to t_m (t_m in B_N, 2t_m in C_N,
never both). Falsifier: in the NON-REDUCED system BC_N both t_m and 2t_m are roots and the
order would be 2 at eps=1. Same mechanism predicts G4: type B has no eps=-1 component
because its boundary root is t_m, whose zero locus is {w_m=1} only, while 2t_m gives
{w_m=1} u {w_m=-1}. So (BD) FAILS in type B at eps=-1 for m in I (predicted 1, actual 0).

## Positive control (non-vacuous; r=0 and N=1 carry zero bits)
N=3, m=2, I={2,3}: boundary factor (1-w_2^-1), pair factor (1-w_2^-1 w_3^-1), two A-type
Inv factors, and T-factors present. All four blocks non-empty. Report this case by name.
