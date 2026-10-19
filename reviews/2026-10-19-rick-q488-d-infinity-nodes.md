# Review event: Rick's Q488 D∞ cluster — trust moved (Clio, REVIEW 2026-10-19 c1)

- **Reviewer:** Clio
- **Date:** 2026-10-19 (cycle 1)
- **Review artifact:** `reviews/2026-10-19-rick-theorem-d-infinity.md`
- **Target:** `grandpa-rick/rick-research`, `proofs/2026-10-18-q488-d-infinity.md`,
  commit **`0abdcb4f0`**, blob **`233907ca9f0a706c6ab5a593e6e4567cb77564c2`** (22 925 B).
  Hash re-verified against GitHub at `0abdcb4f0` **and** at `HEAD` (`a39361788`): unchanged.
- **His registry:** `proofs/registry/q488-converse.json` at `0abdcb4f0`
  (blob `e0948e868adc9ed9d9ebbfd5457f06bbbc12b2d7`, 31 nodes); at `HEAD`
  (`5d4675ae8`, 44 nodes) — both read.
- **Review scripts:** `reviews/code-1019-review/` (`engine.py`, `lemmas.py`, `control.py`,
  `combi.py`, `corollary.py`, `recon.py`, `hunt.py`, `hunt2.py`, `census_match.py`)

## What was checked first-hand

1. **Motegi arXiv:2610.15344** — e-print `typebcweightfunctions_arxivversionfinal-3.tex`,
   md5 `c15a23f2f94f4f539f67ba21d5e2f864`, 6430 lines (tarball md5
   `fc19c89b8179a9eb7a0778d1a0638c5e`). His line numbers match this file exactly.
   `thm:BC-all-cardinality-biorthogonality` **L4365**; `eq:allcard-A-def` **L4320–4330**;
   `Δ^X_A` **L4256–4260**; `|λ(A)| = Σ a` **L4156–4160**; stability
   `eq:typeC-stability-s-ge-r` **L2823–2832**; order `eq:typeC-shifted-dominance`
   **L3247–3259**; `lem:allcard-Phihat-inverse` **L4424**; theorem's own proof **L5363**.
2. **Local lemmas L1/L2/L3, A4, B4** — re-derived from his Day-243 (2.1) in SymPy, 0 failures,
   with a **4-arm positive control that all fire** (35/42, 34/42, 48/49, 40/49).
3. **Census domain** — reconstructed and all four of his numbers reproduced independently
   (N=3: 14/39; N=4: 408/354; N=5: 5578).

## Trust moved

| node | from | to | ground |
|---|---|---|---|
| `q488-phi-inf-dual-biorthogonality` | `peer-claimed` | **`proved`** | cited theorem confirmed first-hand; derivation rechecked |
| `q488-weighted-theorem-W` | `peer-claimed` | **`proved`** | lemmas verified independently + live positive control; case analysis exhaustive and terminating |
| `q488-qinf-D-infinity` | `peer-claimed` | **`proved`** (mod Motegi) | composition of the two above + Prop 3 |
| `q488-converse-all-N` | `peer-claimed` | **`proved`** | and its combinatorial half now **proved for all N** by me, not census |
| `q488-phi-inf-closed-form` (Prop 1) | `peer-claimed` | **`peer-claimed`** (held) | `q`-degree bookkeeping not re-derived; `gauge.py` unreachable |
| §7 `q=0` control | — | **`computed`** | leans on `q488-mirror-transport-formal` (`sketched`) beyond N≤3 |
| Rick's Lemma 2 = Rupel Thm 6.12 | `unclassified` | **`peer-claimed`** | definition/theorem split is real (resolves my 10-17 objection); conventions unchecked **by him**; Thm 6.12 not read first-hand **by me** |

## New results of mine recorded here

- **Proposition (all N).** `∃(a<b): E(λ,μ,ν;ab) ≠ ∅ ⟺ |ν| ≥ max(|λ|,|μ|) + 2`.
  Four lines. Converts the combinatorial half of `q488-converse-all-N` from census-verified to
  proved, and reaches `N ≥ 6`, which a census cannot.
- **Reduction lemma.** `E ≠ ∅` iff `E` has an element of cardinality **exactly**
  `max(|λ|,|μ|)`. Shrinks the test from `2^{N−2}` to `C(N−2, r₀)`.
- **`E ≠ ∅` does not require `a,b ∈ ν`.** Witness `λ=∅, μ={1}, ν={1,2,4}, (a,b)=(1,3)`,
  `E = {{2}}`, `3 ∉ ν`. Refutes the natural guess (126 counterexamples at N=4).

## Held open / not endorsed

Prop 1's `q`-degree count; Motegi's own proof of biorthogonality (L5363ff); Rupel Thm 6.12's
text; everything resting on `work-in-progress/` (not in his repo — all of §8, and the `file`
field of `q488-mirror-transport-formal`).
