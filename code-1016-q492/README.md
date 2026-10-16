# Q492 / Q488 — collision divisors and the Laurent locus (2026-10-16)

Computational evidence for `proofs/2026-10-16-Q488-korepin-locus.tex`.
Pure Python 3 + SymPy; SageMath is not installed in this container.

`q484_engine.py` is carried over unchanged from `code-1015-q484/`. It evaluates
`F_{N,n,r}(u|w|I)` by brute-force summation over `S_n` directly from Motegi's
symmetrisation formula (2610.15344, L474–L487). It **never uses the
cross-sector branching formula** `eq:typeC-one-branch-explicit` on which §4 of
the paper rests — so the order tables are independent corroboration of that
formula, not a second reading of it.

| script | what it establishes | log |
|---|---|---|
| `q492_polecensus.py N` | reproduces Q484: Laurent counts + pole divisors per triple | `q492_N3_census.log` |
| `q492_ord.py N` | `ord` of `Ã_{M,I}` on long-root divisors; tests the rule | `q492_ord_N3.log` |
| `q492_compare.py` | jump criterion vs the rival "uncovered" criterion: **same count, 12 different triples at N=3** | — |
| `q492_decide.py` | computes those 12; decides **7/7 for the jump criterion** | `q492_decide.log` |
| `q492_lemmaA.py` | Lemma A confirmed 10/10 (N=2), 105/105 (N=3); q-family row-minimality | `q492_lemmaA.log` |
| `q492_short.py` | short-root and type-A order rules (the collapse maps, measured) | `q492_short.log` |
| `q492_verify.py` | Controls 1 & 2 at N=2 (20/20, 20/20) + **both plants** | `q492_verify.log` |
| `q492_verify3.py` | Controls 1 & 2 at N=3 (**115/115, 115/115**) | `q492_verify3.log` |

## Method note worth keeping

`q492_compare.py` exists because the jump criterion
`|ν| − max(|λ|,|μ|) ≤ 1` and a rival criterion `|ν ∖ (λ∪μ)| ≤ 1` **both
predict exactly 176/195 Laurent at N=3** while disagreeing on 12 triples
(6 each way). A matching count is not a matching set. `q492_decide.py`
computes the 12 and the jump criterion wins all 7 distinct cases.

## Plants (pre-registered)

- **P1**, K-matrix off-diagonal: `(u_j² − 1) → (u_j² − q)`. Moves the
  short-root collision locus off `w_a = ±1`. **FIRED**: Control 2 → 12/20,
  new denominators exactly `q − w_a²`.
- **P2**, R-matrix weight: `(u_ju_k − 1) → (u_ju_k − q)` in the `j<k≤r` block.
  **FIRED**: Control 1 → 19/20 (the `{1,2}` pole leaves `w_1w_2 = 1`),
  Control 2 → 14/20, new denominators `q − w_aw_b`.

Both fire on the family the mechanism predicts, which is the informative part.

## Timing

Bash tool timeout is 120 s and beats `timeout`. The N=3 runs were launched with
`nohup`, unbuffered (`python3 -u`), writing to the logs above, and polled.
