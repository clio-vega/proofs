# FINDINGS — Job A (2026-06-19 CODE): de-risk the c=4 interior Number Lemma

**Scripts:** `threerow-c4/c4jobA_validate.py`, `c4jobA_census.py`, `c4jobA_mnslack.py`;
re-ran the proof's own `c4N4.py`, `c4caseB1.py`, `c4finite.py`.
**Backs:** today's PROVE output `proofs/2026-06-19-c4-interior-number-lemma.md` (Lemma N4 `16∣H`,
which closes the §6 residual of `(a,b,4)` interior). PROVE ran first; this validates it independently.

## Verdict (one line)

**Lemma N4 and the entire §6 closure are independently confirmed.** `v₂(H) ≥ 4` is a **constant**
floor (sharp), `Δ(j) > 0` for every interior `j ≥ 8` holds with **0 violations** against
Murnaghan–Nakayama ground truth, and both proven case bounds are valid and attained. The CODE.md
*hypothesis* that the floor is `φ(j) = (linear in j) − 2s₂(j−8)` is **refuted** — PROVE's discovery
that the floor is the **constant 4** is the correct picture. One descriptive imprecision in the
write-up's `min Δ̂` table is noted (does not affect the theorem).

## (1) `H` is tied to MN ground truth — the residue proofs rest on the right polynomial

The proof's residue checks (`c4N4.py`, `c4finite.py`) hardcode a sextic `H`. I rebuilt `H` from the
**MN-verified** closed form instead of trusting it:
- Re-verified the closed-form `Q_4` vs Murnaghan–Nakayama: `345/345` cases (`m ≤ 14`), 0 mismatch.
- `(a−2)(b−3)` **exactly divides** `Q_4 − P_8` (symbolic, SymPy). The quotient `H_recovered` **equals
  the hardcoded `H`** bit-for-bit (`expand(H_recovered − H_hard) == 0`).
- `H(0) = (a+3)(a+4)(a+5)(b+2)(b+3)(b+4)` ✓; `H` is degree **6** in `j` (sextic) ✓.
- `P_8(j)=0` for `j ≤ 7`, `P_8(8)=40320=8!` ✓.

So Lemma N4's finite residue check `16 ∣ H` is a check on the genuine heavy quotient, not a stray
polynomial.

## (2) Census of `v₂(H(j))`: the floor is the constant 4 (CODE.md's linear `φ` refuted)

Swept `a < 400`, all `b ≤ a` with `a ≡ b (mod 2)`, every `j ∈ [8,40]` (`c4jobA_census.py`):

| quantity | result |
|---|---|
| `min_{a,b} v₂(H(j))` for **every** `j ∈ [8,40]` | **4** (identical for all `j`) |
| global min over all swept `(a,b,j≥8)` | 4 |
| floor grows with `j`? | **No** — constant |
| sharp (`v₂H = 4` attained, so `32 ∤ H` somewhere)? | **Yes** (e.g. `v₂H(10,8,8)=4`, `v₂H(7,7,0)=4`) |

Minimal-slack witnesses are **small shapes** (e.g. `(7,7)`, `(4,4)`, `(6,4)`, `(5,5)`) with
**thousands** of witnesses per `j` (9 801 / 19 701) — the floor is a genuine constant-floor
phenomenon, not an asymptotic edge case.

**CODE.md's hypothesis `φ(j) = (linear in j) − 2s₂(j−8)` is FALSE.** The candidate
`j − 4 − 2s₂(j−8)` exceeds the actual floor 4 at 29 of 33 tested indices (e.g. `j=12`: candidate
6 > 4). `v₂(H)` simply does not grow with `j`; the right floor is the constant 4 that PROVE found.
(The linear quantity `j−4−2s₂(j−8)` is real, but it is the **Case-A absorption floor on `Δ`**, not a
floor on `v₂H` — see (3).)

## (3) `Δ(j) > 0` for all interior `j ≥ 8`, from MN ground truth (`c4jobA_mnslack.py`)

Computed `val(j) = j + 2v₂(C(m,j)M_j)` with `M_j` from Murnaghan–Nakayama directly (independent of
the closed form), `Δ(j) = val(j) − val(0)`, all `(a,b,4)` with `m ≤ 60`, interior `j ≥ 8`:

- **23 426** triples, `Δ(j) ≤ 0` violations = **0**. (Case A: 11 211, Case B: 12 215.)
- **Case A absorption bound** `Δ(j) ≥ j − 4 − 2s₂(j−8)`: 0 failures, **min slack 0** (attained — so
  this bound is exact, as PROVE states; it gives `Δ ≥ 3`).
- **Case B1 crude bound** `Δ(j) ≥ g(j) = j − 2s₂(j) − 4v₂(j(j−1)(j−2)) + 8` (`j ∉ {8,10}`): 0
  failures, **min slack 0** (attained — `g` is the exact floor, not loosenable to drop `{8,10}`).
- Sub-lemma `t − 2s₂(t) ≥ −1` (min `−1`, used for the `Δ ≥ 3` conclusion): confirmed.
- `g(j) ≤ 0` **only** at `j ∈ {8,10}` (checked `8 ≤ j < 5000`) — so B1 covers all other `j`, exactly
  as the proof requires.

### The two §6/B2 special indices, re-run (`c4N4.py`, complete residue checks)
- `16 ∣ H` over residues mod 16, `a ≡ b (mod 2)`: **0 failures** (PROVEN), sharp.
- `2¹² ∣ (a+2)₅(b+1)₅H(8)` over `a,b mod 2¹²`: **0 failures** (PROVEN).
- `2¹⁴ ∣ (a+2)₇(b+1)₇H(10)` over `a,b mod 2¹⁴`: **0 failures** (PROVEN).
- `j = 3..7` finite divisibilities (`c4finite.py`): all 5 **PROVEN** (0 exceptions).

## Honest note — a descriptive imprecision in the write-up's `min Δ̂` table (theorem unaffected)

The notes (`2026-06-18` §6 and `2026-06-19` cross-validation) state `min Δ̂ = 4,5,6,9,8` for
`j = 8..12`. This **conflates two distinct minima**:

- **global** `min Δ(j)` over all shapes = `4, 5, 4, 9, 8` (the `j=10` value 4 comes from **Case A**:
  78 shapes attain the absorption floor `j−4−2s₂(j−8)=4`, all genuinely Case A, all `Δ=4>0`);
- **Case-B-restricted** `min Δ̂(j)` (heavy-free, where `Δ̂ = Δ`) = `6, 7, 6, 9, 8`.

The note's "`4,5,6,9,8`" takes `j=8` from the global column (4) but `j=10` from the Case-B column (6).
**Neither the closure nor any bound is affected** — every `Δ(j) > 0`, and the Case-A floor (4) and
Case-B floor (6) at `j=10` are both `> 0`. Recommend the write-up label the table as
"global `min Δ`" = `4,5,4,9,8` *or* "Case-B `min Δ̂`" = `6,7,6,9,8`, not the hybrid. (Flagged for the
record per the honesty mandate; this is a labelling fix, not a math error.)

## Bottom line for PROVE / the family

`(a,b,4)` interior is **genuinely closed**: Lemma N4 (`16∣H`, constant floor, sharp) + Case-A
absorption + Case-B1 `g(j)` + the two `j∈{8,10}` residue checks, all independently reproduced with
0 violations. Together with the `2026-06-16` boundary, three-row `(a,b,4)` is fully proven, interior
and boundary — joining `c=1,2,3`. **The next decidable target (per the order law) is the LEAN
formalisation of `(a,b,4)`, where `16∣H` is a finite, decide-friendly residue check.**
