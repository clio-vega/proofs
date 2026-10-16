# PROVE 2026-10-16 c2 — the twisted collapse on a long root

Computational evidence for `proofs/2026-10-16-c2-longroot-multiplier.tex`.
Pure Python 3 + SymPy; SageMath is not installed in this container.

`q484_engine.py` is carried over unchanged from `code-1016-q492/`. It evaluates
`F_{N,n,r}(u|w|I)` by brute-force summation over `S_n` directly from Motegi's
symmetrisation formula (2610.15344, L474–L487). It **never uses the branching
formula** `eq:typeC-one-branch-explicit`, on which every proof in the paper
rests — so these tables are independent corroboration, not a second reading of
one engine.

| script | what it establishes | log |
|---|---|---|
| `gate_q0.py N` | Task-0 gate: generic pole census vs the `q→0` limit, plus `ord_q` of each residue | — |
| `gate2.py N` | the decisive columns: `ord_q(c)`, `ord_q(res)`, multiplicity. **25/25 all equal 1 at N=3** | `gate2_N3.log` |
| `collapse_test.py N X` | **Theorem A**, with 2 plants. 48/48 type C and 48/48 type B at N=3 over a 12-item non-degenerate domain; both plants fail 12/12. **N=4 C claim arm: 384/384 over 94 non-degenerate**; its plant arms had not returned at session end | `collapse_N3_C.log`, `collapse_N3_B.log`, `collapse_N4_C.log` |
| `psi_ratio_test.py N` | **Theorem A's *proof***: the factor-by-factor ratio table, with an **off-H control** that must fail. 1/1 at N=2, control fires | — |
| `nonvanH_test.py N` | **Gap nonvanH**, with a discriminating positive control. 12/12 + 52/52; control 41/41 over a domain holding both values | `nonvanH_N3_C.log` |
| `residue_test.py N` | **Theorem D**'s closed-form residue, with 3 plants in 2 blocks. **19/19 at N=3**, `ord_H = -1` in all | `residue_N3_C.log` |
| `coverage.py` | combinatorial coverage of Theorem D over jump-≥2 triples, N=2..6 | — |

## The two numbers that matter

`ord_q(c) = ord_q(residue) = 1` in **25/25** pole records. That is the whole
reconciliation with Rick's M2: the non-Laurent constants are exactly the ones
that *vanish* at `q=0`, so "Laurent at `q=0`" is satisfied by vanishing, not by
the pole cancelling. The gate script therefore divides by `q^{ord_q}` **before**
setting `q=0` — the naive substitution scores a vanishing constant as Laurent.

`rho - 1 = q (w_a - w_a^-1)^2 / ((1-q w_a^2)(1-q w_a^-2))` has exactly one
factor of `q`. That was derived from the branching formula and only then
compared with the 25 records: a prediction, not a fit.

## Domain hygiene

Every table reports its **non-degenerate domain size** beside the pass count.
At `N=2` the collapse test has 4 pairs but only **one** where either side is
nonzero on `H`, so the `N=2` run alone cannot distinguish any multiplier — it is
decoration, and `N=3` (12 items) is the first real test. Same reason the
`nonvanH` control prints both the zero count and the positive count: a control
returning one value over its whole domain is not evidence that it can detect a
violation.

## Still running at session end

`collapse_test.py 4 C` (claim arm done: **384/384**, 94 non-degenerate; plant
arms pending), `nonvanH_test.py 4 C` and `psi_ratio_test.py` at N≥3 had not
finished (SymPy `cancel` over 4 variables is the bottleneck; the `N=4` build
alone is 73 s and the `ord_H` loop is ~1500 rational-function reductions).
**No N=4 claim is made in the paper.** Logs `collapse_N4_C.log`,
`nonvanH_N4_C.log` will hold the answer if the processes complete.

## Timing

Bash tool timeout is 120 s and beats `timeout`. The N=3 runs were launched with
`nohup`, unbuffered (`python3 -u`), writing to the logs above, and polled.
`residue_test.py 3 C` takes 462 s; `gate2.py 3 C` a few minutes.
