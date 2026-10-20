# `code-1020-LS` — computations for `2026-10-20-hypothesis-LS.tex`

Pure Python + SymPy, exact rational arithmetic. Engine `q484_engine.py` is the
canonical copy (`md5 f53c0258…`), my own transcription of Motegi's
Motherfunction (`2610.15344` L473–490) and normalisers.

| script | what it measures | result |
|---|---|---|
| `cells.py` | gate 0b: split `hyp:L^S` into cells **before** measuring it | reproduces the 10-19 paper's `699 = 3×233` and plant `1179 = 3×393` from combinatorics alone |
| `probeR.py` | route (iv) as a *rank reduction* (the thing I did **not** need) | 135/135 nonzero; reduction is clean iff `m = N`, B-dependent for `m < N` |
| `degen.py` | the proof's mechanism: `Ã → [A≥B]` at `q=0` under one spectator + hierarchy | `--maxN 5`: **9558/9558**; `--plant-ord`: plant **fires** (12 deviations, errors grow) |
| `terms.py` | term-by-term: `ord_T(pref·t_σ) = Pen(σ)`, dominant coeff `+1`, multiplicity `(n−r)!` | **1119/1119**, 612/612 predicted zeros, 819 dominant |
| `spect.py` | the general multi-spectator form, every `Σ ⊆ [N]`, adversarial spectator values | **3108/3108** |
| `denom.py` | Lemma 2.1(2): `den(Ã)|_{q=0}` nonzero and supported on `(w_a)_{a∈A}` | **96/96** (all 96 in fact monomials — stronger than the lemma claims) |
| `fullq_N5.py` | **full-`q`** direct `ord_t` census at `N=5`, gated on (BD) at the sample point | **2622/2622**, plant **5058/5058** |

Logs: `degen_N4.log`, `degen_N5.log`, `plant_ord_N3.log`, `fullq_N5.log`.

Two instrument faults, both mine, both recorded in §4 of the paper:
`degen.py`'s first run reported 11 false failures from a fixed `1e-3` threshold
against a true `O(T)` decay (now graded on the decay *exponent*); `terms.py`'s
first run reported 243 false failures from my demanding leading coefficient `+1`
for *sub*dominant terms, which is not claimed and is false.
