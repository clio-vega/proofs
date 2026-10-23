# Q559 PRE-REGISTRATION — written 2026-10-23 before testing the predicates

## P1: "the extremal-monomial argument closes the general case"
PREDICTED ARM: **SPLIT**.
Reasons, written before computing:
(a) The brief's strategy-1 as literally stated is **ill-posed**: the summands
    B(T)*det/V are NOT polynomials in w (B(T) carries the denominator
    prod_{t notin T} prod_{x in T}(1 - w_t/x)), so "a monomial occurring in
    exactly one T-summand" has no meaning until a common denominator is cleared,
    and clearing it destroys the T-indexing.
(b) I expect a leading-term argument to work in a DIFFERENT grading: the order of
    vanishing at w_n -> infinity in the single variable w_n, over the field
    K = Q(w_1,...,w_{n-1}).  I predict the T-containing-n block has order
    w_n^{max(m) - r + 1} and the T-avoiding-n block order exactly w_n^{-r}, so the
    first block strictly dominates iff max(m) >= 0.
(c) I therefore predict a residual case max(m) <= -1 that the leading-term
    argument does NOT close, to be covered by the reflection w -> 1/w.

## P2: "B(T) genuinely takes both signs for some r=2 instance in the test range"
PREDICTED: **TRUE** (B takes both signs).  B contains (1-x^{-2}), negative for
0<x<1 and positive for x>1, so a mixed sample should already do it.

## P3: "nonvanishing is generic, not universal"
PREDICTED: **CONFIRMED**, and I already hold a witness found while diagnosing a
58/60 negative control: w = (6, 5/3, 3/5) has the reciprocal pair w_2 w_3 = 1,
which makes EVERY B(T) vanish when n = r+1, hence Psi = 0 for EVERY theta.
Structural prediction to test: if M contains a,b with w_a w_b = 1 and |M| <= r+1,
then Psi_{(1^r),theta}(M) = 0 for all theta; if |M| >= r+2 it need not vanish.

## New prediction P4 (mine, not the brief's)
LEADING-COEFFICIENT LEMMA.  With m pairwise distinct and j* the unique argmax of m,
in K((1/w_n)):
    Psi_{(1^r),theta}(M) = w_n^{m_{j*}-r+1} [ (-1)^{r+j*} Psi_{(1^{r-1}),hat-theta}(M') + O(1/w_n) ]
                           + O(w_n^{-r}),
where M' = M minus w_n, hat-m = m with the j*-th entry deleted (order preserved),
hat-theta_l = hat-m_l - l + 1.  Signs in 1-indexed j*.
ARMS: CONFIRMED / wrong sign / wrong exponent / REFUTED.

## P5 (mine) REFLECTION IDENTITY
  Psi_{(1^r),theta}(1/w) = (-1)^r e^{2r} Psi_{(1^r),theta^vee}(w),  e = prod_{t in M} w_t,
  m^vee_j = -m_j, i.e. theta^vee_j = -theta_j - 2j + 2.
ARMS: CONFIRMED / REFUTED.

## Hypothesis the brief omits
The conjecture needs |M| >= r.  If |M| < r there are no distinct r-tuples and
Psi = 0 identically with m pairwise distinct.  Record as a boundary case, not a
refutation.
