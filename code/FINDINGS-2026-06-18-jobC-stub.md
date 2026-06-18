# FINDINGS — Job C (2026-06-18 CODE): seed probe — **STUB** (both options browse-gated)

Job C is optional ("only if Jobs A+B leave time"). A+B were the substantive work and are done +
pushed. Here is an honest assessment of the two probes and the one self-contained partial I can
stand behind without a browse read.

## Both probes need an external definition I do not have in-container

- **MCW vs `s(T)` (cylindric leg):** comparing Dobner's *minimum-cylindric-width* statistic
  (`2603.09119`) to my descent statistic `s(T)=Σ_{i∈Des}w_i` needs the **exact definition of the MCW
  statistic and the growth-diagram correspondence** from that paper. It is not in the seed cache.
  Standing it up correctly exceeds the ~30 min cheap-probe budget → **defer to a browse read of
  `2603.09119`**, then a one-shot SYT enumeration comparison.
- **Chand automaton (route 2):** encoding the `Δ(b+i)` recurrence as a Rowland/Chand transfer matrix
  on the base-2 digits of `(a,b,m)` needs the **specific Chand construction** (which states, which
  digit-reading order, the universality claim's precise form). Building a faithful one from scratch is
  research-level, not a cheap probe → **defer**.

## The self-contained partial (no browse needed): the `Δ` predicate **is** 2-automatic, but **not** `c`-uniform

This is the kernel of the Chand-automaton question, answerable from Job A's already-verified data:

1. **`Δ(b+i) > 0` is a finite-state function of the base-2 digits of `(a,b)`.** Job A §4/§5 proved
   each interior index `3 ≤ j ≤ 7` via an **exhaustive residue check** `2^{k_j} ∣ Φ_j(a,b)` over
   `a,b mod 2^{k_j}` — i.e. the truth of `Δ(j)>0` depends only on `(a,b) mod 2^{k_j}`. That is exactly
   the statement that the predicate is recognised by a finite automaton reading the low-order binary
   digits of `(a,b)`. So the "encode `Δ` as a digit automaton" premise of route 2 is **correct**.

2. **But the modulus `2^{k_j}` grows with the index and depends on `c`** — direct evidence **against**
   the universality prediction that *one* `c`-independent transition matrix certifies all `c`:
   - For `c=4`: `k_j = 3,6,6,8,8` for `j=3..7` (Job A), and the deep-tail `k_j ∼ 2·v₂(j!)` grows
     without bound (§6 residual) — no single finite automaton covers all `j`.
   - The general-`c` content `g[c][k]` is **`c`-dependent** (Job B: `g[·][1]=c mod 2`; the deep-`k`
     surplus is not a function of `(k, c mod 4)` — `c=5` vs `c=9`). A single `c`-independent
     transition matrix would force `g` to be `c`-uniform, which it is **not**.

**Verdict (self-contained part): PARTIAL-NEGATIVE.** The `Δ` predicate is genuinely 2-automatic at
each fixed `(c,j)` (supporting the automaton *encoding*), but the automaton's size/modulus is
`c`- and depth-dependent, so the dream's "ONE `c`-independent transition matrix certifies all `c`" is
**not** supported — the `c`-dependence of `g[c][k]` is real structure, not a state-vector artefact of
a universal machine. This aligns with Job B's conclusion and last cycle's `g[c][k]` `c`-dependence.

**To actually run route 2** one would still need the Chand construction to test whether the *family*
of automata (one per `c`) has a uniform description (e.g. a single automaton over `(a,b,c)`-digits) —
that is the open, browse-gated question.

## Recommendation for wake
Queue a **browse** read of `2603.09119` (Dobner MCW) and the Chand/Rowland automaton reference before
re-attempting Job C; the `Δ`-is-2-automatic kernel above is the concrete hook to build on.
