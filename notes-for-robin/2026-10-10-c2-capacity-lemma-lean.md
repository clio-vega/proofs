# The capacity lemma is machine-checked — and the stretch goal I was given is false

**Readable here:**
https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/CapacityL1.lean
Notes: https://github.com/clio-vega/tworow-d4-kernel/blob/main/NOTES-2026-10-10-c2-capacity-l1.md
Commit `8f8fd82`, verified reachable from `origin/main`.

`lem:cap` from yesterday's Q405 note is now 7 Lean declarations, 0 sorries, standard three
axioms on all 7. The pleasing part is how little of it is about cylindric tableaux: the
whole lemma collapses to one fibre count — *the fibres of a statistic over distinct values
are disjoint subsets, so their sizes sum to at most the size of the set* — and both parts
of the lemma plus both instantiations fall out as corollaries. A monomial weight can only
redistribute objects among powers of `t`. It cannot create them. That sentence is now a
thing I can `apply`.

**The part worth your attention.** The brief's stretch goal asked me to show the signed gap
is "exactly 2 for every `M ≥ 2`". It isn't. The arithmetic in the brief is correct —
`‖c_α‖₁ = C_M + 2M` and `C_M + 2M − 2` differ by 2 — but `C_M + 2M − 2` is `A_α + B_α`,
**not** `TOT_α`. The two coincide only at `M = 2`, because that is the only place where no
shape has `c⁻ = 0`. From `M = 3` the `c⁻ = 0` shapes — the ones HKKO's theorem assigns
weight zero, and therefore cannot see — start contributing to `TOT_α`, and at `M = 3` it is
`18` against `‖c_α‖₁ = 11`. The inequality reverses. The signed obstruction bites at
`ℓ = 0` and nowhere else, which is exactly what my own remark four lines below the theorem
says, under the heading *"why `ℓ = 0` is the only place the stronger form bites"*.

So I did not need a computation to catch this; I needed to read four lines past the theorem
I was quoting. The uniform-in-`ℓ` statement does exist — it is part (1), where the deficit
is exactly one tableau at every `ℓ` — and that is formalised as `no_strict_statistic`. The
brief had the right shape and attached it to the wrong half of the lemma.

I did not patch the false target. Per the session rules, I formalised the true neighbour and
recorded the refutation in the commit message, the module docstring, and here.

One small instrument finding you may care about: `registry_validate.py` has the same
double-prefix fault as `trustcheck.py` — `--proofs-dir` defaults to the registry's parent,
which is already `proofs/`, so every node reports "file not found". 27 false problems, 0
with `--proofs-dir .`. Third tool in this family with the fault.
