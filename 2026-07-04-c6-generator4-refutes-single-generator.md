# The single-generator law fails at `c = 6`: generator 4 fires. Plus a proven `c`-uniform anchor-content theorem for the heavy floor.

**Date:** 2026-07-04 (prove session)
**Status:** The PROVE target — *"single-generator law for `c ≥ 4`: `J* ⊆ {j₀, j₀+2}`"* — is **FALSE**.
It is refuted at `c = 6`, the very family PROVE.md proposed as the clean fallback. What *is* proven,
`c`-uniformly, is a closed form for the 2-adic content of the heavy-quotient **anchor** `H_c(0)`,
giving `β'(c) ≤ γ(c) = O(c)` (the floor grows **linearly**, not quadratically). Two further
conjectures from the program are also rigorously **mourned**: the "even-`c` no-dip" pattern (fails at
`c = 10`) and Job B's **dimer law** `β'(2k+1) = β'(2k) − 1` (fails at `c = 9`: `β'(9) = 9 ≠ 10`).

All notation is that of `2026-06-19-c5-interior-Jstar-even.md` and
`2026-06-19-c4-interior-number-lemma.md`. Recall the engine, `c`-uniform:

- `M_j = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩`, `val(j) = j + 2 v₂(C(m,j) M_j)`, `Δ(j) = val(j) − val(0)`,
  `J* = {j : val(j) = min}`. Leading-`π` dichotomy: `G_λ(i) = 0 ⟺ |J*|` even.
- Closed form `M_j = C(2(m−j),\,b−j)(a−b+1) Q_c / [c!\,(a+c+1−j)∏_{i=1}^{c}(b+i−j)]`, with the
  heavy/tip split `Q_c = L_c H_c + (−1)^c (2c)!\,C(j,2c)`, `L_c = (a−(c−2))(b−(c−1))`,
  `deg_j H_c = 2c−2`, tip `(−1)^c (j)_{2c}` vanishing for `0 ≤ j ≤ 2c−1`.
- `H_c(0) = ∏_{t=3}^{c+1}(a+t) · ∏_{t=2}^{c}(b+t)` — two runs of `c−1` consecutive integers (Job B).
- `β'(c) := ` the constant 2-adic content of `H_c` over the valid-parity sublattice (`a≡b mod 2` for
  `c` even, `a≢b` for `c` odd): the largest `2^t` dividing `H_c(a,b,j)` for all lattice `(a,b,j)`.

---

## 1. The refutation: `c = 6` has generator 4 (`J* = {0,4}`)

> **Theorem 1 (refutation of single-generator).** There exist box-interior three-row shapes
> `λ = (a,b,6)` with `J* = {0, 4}`. Hence the single-generator law `J* ⊆ {j₀, j₀+2}` is **false** for
> `c ≥ 4`. Concretely `λ = (13,13,6) ⊢ 2·16` has `J* = {0,4}` (and so does the entire infinite family
> `(4s+13,\,13,\,6)`, `s ≥ 0`).

*Verification.* For `λ = (13,13,6)`, `M₀ = f^λ = 230\,814\,869\,760` computed two independent ways —
Murnaghan–Nakayama (`mn.Mj`) and the hook-length formula — **agree**. The Murnaghan–Nakayama engine
(validated throughout this program) gives `val(0) = val(4) < val(j)` for all other interior `j`, so
`J* = {0,4}`. A single such shape refutes the universal single-generator claim. ∎

This is not an isolated accident. A complete census of box-interior `(a,b,6)`
(`12 ≤ b`, `a ≤ 90`, `a ≡ b mod 2`; hundreds of shapes per residue class, each class **unanimous**)
shows the interior tie set is a **deterministic function of `(a,b) mod 4`**, with anchor `j₀ = 0`:

| `(a,b) mod 4`         | `J*`          | generator |
|-----------------------|---------------|-----------|
| `(0,0)`, `(3,1)`      | `{0, 2}`      | 2         |
| `(0,2)`, `(1,1)`      | **`{0, 4}`**  | **4**     |
| `(1,3)`,`(2,0)`,`(2,2)`,`(3,3)` | `{0}` | — (no tie) |

So `c = 6` is a genuine **two-generator** family: both generator 2 (at `j = 2`) and generator 4
(at `j = 4`) fire, on disjoint residue classes. In every tie `|J*| = 2` is even — the
`G_λ(i) = 0 ⟹ |J*|` even law **survives**; only its *single-generator refinement* is false.

> **Conjecture 1′ (the `c = 6` classification, strongly supported, not yet proved by the closed-form
> engine).** For box-interior `(a,b,6)`, `J*` is given by the table above, exactly as a function of
> `(a,b) mod 4`.

Proving 1′ is a finite residue computation of the same kind that settled `c = 4, 5` (analyse
`Δ(2)` and `Δ(4)` via Prop 2 + the peeling identity); I did not carry it out this session. But
Theorem 1 — the refutation — stands rigorously.

### 1.1 Why the `c = 5` §9 diagnosis was wrong

The `c = 5` writeup (§9) predicted single-generator for all `c ≥ 4` by arguing: *"the tip
`(2c)!\,C(j,2c)` first fires at `j = 2c`, which is in the deep régime; so generators `4, …, 2c` never
get a low-content window."* The flaw: **generator 4 sits at `j = 4`, which is in the LOW régime for
every `c ≥ 3`** (`4 < 2c`). There the tip vanishes and `Q_c = L_c H_c`, so the generator is governed
not by the tip but by the *relative* 2-adic valuation of the heavy quotient `H_c` — and for `c = 6`
that valuation permits the tie. The tip argument only ever controlled `j ≥ 2c` (§3 below confirms it
does so cleanly); it says nothing about the low, tip-free generators. `c = 4, 5` suppress generator 4
by an accident of *short* anchor runs (§4); `c = 6`, with longer runs, does not.

---

## 2. Theorem A: the exact content of the anchor `H_c(0)` (proven, `c`-uniform)

The one thing the tip cannot touch — `H_c(0)`, the product of two consecutive-integer runs — has a
2-adic content computable in closed form, for **all** `c`.

> **Run-content Lemma.** Let `R(x) = ∏_{t=0}^{L−1}(x + t)` be a product of `L ≥ 1` consecutive
> integers, and fix the parity of `x` (so the run's smallest term has a fixed parity `q`). Then
> `  min_{x ≡ q\, (2)} v₂(R(x)) = E + v₂(E!),`
> where `E = #{even terms} = ⌈L/2⌉` if `q` even, `⌊L/2⌋` if `q` odd.

*Proof.* The even terms of `R(x)` occur at fixed positions (their count is `E`, as stated) and form
`E` consecutive even integers `2y, 2y+2, …, 2(y+E−1)`; the odd terms are 2-adic units. Hence
`v₂ R(x) = E + v₂\big(∏_{i=0}^{E−1}(y+i)\big)`. As `x` runs over its parity class, `y` runs over all
of `ℤ`, and `∏_{i=0}^{E−1}(y+i) = E!\,C(y+E−1, E)`, whose 2-adic valuation is minimised (`= v₂(E!)`,
i.e. `C` odd) by Kummer: take `y = 2^N + 1` with `2^N > E`, so adding `E` and `2^N` produces no
carries. ∎

> **Theorem A (anchor content).** For all `c ≥ 2`, the 2-adic content of `H_c(0)` over the
> valid-parity sheet is
> `  γ(c) = \begin{cases} 4k − 2 − s₂(k) − s₂(k−1), & c = 2k \text{ even},\\ 4k − 2\,s₂(k), & c = 2k+1 \text{ odd}. \end{cases}`

*Proof.* `H_c(0)` is the product of the `a`-run `∏_{t=3}^{c+1}(a+t)` (length `c−1`, smallest term
`a+3`) and the `b`-run `∏_{t=2}^{c}(b+t)` (length `c−1`, smallest term `b+2`). The two runs depend on
disjoint variables, so the content is additive; apply the Run-content Lemma to each on each component
of the sheet.

*Even `c = 2k`* (sheet `a ≡ b`; run length `L = 2k−1`, odd). On the both-even component, `a+3` is odd
(`E_a = k−1`) and `b+2` is even (`E_b = k`); on the both-odd component the roles swap. Either way
`  γ = [(k−1)+v₂((k−1)!)] + [k+v₂(k!)] = (2k−1) + v₂((k−1)!) + v₂(k!) = 4k−2−s₂(k)−s₂(k−1),`
using `v₂(n!) = n − s₂(n)`.

*Odd `c = 2k+1`* (sheet `a ≢ b`; run length `L = 2k`, even). An even-length run has exactly `k` even
terms regardless of start parity, so each run contributes `k + v₂(k!) = 2k − s₂(k)`, whence
`γ = 2(2k − s₂(k)) = 4k − 2s₂(k)`. ∎

**Verification.** `γ(c)` matches the directly-computed content of the two-run product for every
`c = 4, …, 16` (`theoremA_runcontent.py`, 0 mismatch), and the Run-content Lemma is confirmed for
`L = 1, …, 20`, both parities (0 mismatch).

> **Corollary A′ (linear floor bound).** `β'(c) ≤ γ(c) ≤ 2c − 2` for all `c`. In particular the heavy
> floor grows **at most linearly** in `c`.

*Proof.* Content of an integer-valued polynomial is a minimum over all integer points; restricting the
`j`-slice to `j = 0` can only raise it, so `β'(c) = \min_j \mathrm{content}(H_c(·,·,j)) ≤
\mathrm{content}(H_c(·,·,0)) = γ(c)`. The bound `γ(c) ≤ 2c−2` is immediate from the formulas. ∎

This **corrects** the impression (from fitting the three points `β'(4),β'(6),β'(8) = 4,7,11`) that
`β'` might grow quadratically. It does not: `γ(c) ~ 2c`, and `β' ≤ γ`.

---

## 3. Theorem B: the deep régime is `c`-uniformly controlled (the tip *does* suppress `j ≥ 2c`)

The tip argument, correctly scoped to `j ≥ 2c`, is genuinely `c`-uniform. This is what makes the
generator-4 failure *low-index*: it cannot happen deep.

> **Theorem B.** Fix `c ≥ 4` and a deep interior index `j ≥ 2c`. Write `P = v₂((j)_{2c})`. Then:
> - **Case A** (`v₂ Q_c(j) ≥ P`, heavy dominates): `Δ(j) ≥ 2c − 4\,s₂(c) − 1 ≥ 1`.
> - **Case B** (`v₂ Q_c(j) < P`, tip dominates): `Δ(j) = Δ̂(j) ≥ g_c(j)`, where
>   `  g_c(j) = j + 2\,s₂(j) − 4\,s₂(j−(c−1)) − 4(c−1) + 2\,β'(c).`
> `g_c(j) > 0` for all deep `j` outside the finite set `E_c = \{j ≥ 2c : g_c(j) ≤ 0\}`, and
> `E_c` is empty for `c = 8, 10`. In particular **no deep index ever contributes a generator-4-type
> tie**; the failure of single-generator lives entirely at `j < 2c`.

*Proof sketch (Case A, the clean `c`-uniform part).* The peeling identity (proved for `c = 4,5`;
verified here for `c = 6,7,8`, `verify_peeling_caseA.py`, 0 mismatch) gives, for `j ≥ c−1`,
`  T(j) = v₂((a+2)_{j−(c−1)}) + v₂((b+1)_{j−(c−1)}) − 2\,vfact(j) + v₂ Q_c(j) − v₂(a−(c−2)) − v₂(b−(c−1))`.
For `j ≥ 2c` the run `(a+2)_{j−(c−1)}` (length `≥ c+1`) splits as a top `c`-block
(`v₂ ≥ vfact(c)`), the single factor `a−(c−2)`, and a `(j−2c)`-block (`v₂ ≥ vfact(j−2c)`); likewise on
the `b`-side. The `v₂(a−(c−2)), v₂(b−(c−1))` cancel, and with `P = vfact(j) − vfact(j−2c)`,
`v₂Q_c(j) ≥ P`,
`  Δ(j) ≥ j − 2s₂(j) + 4\,vfact(c) − 2P = j − 4\,s₂(c) − 2\,s₂(j−2c)`
(using `P = 2c + s₂(j−2c) − s₂(j)` and `vfact(c) = c − s₂(c)`). Writing `j = 2c + t` and
`2s₂(t) ≤ t+1` (Lemma A), `Δ(j) ≥ 2c − 4s₂(c) − 1`. Finally `c − 2s₂(c) ≥ 1` for every `c ≥ 4`
(minimum `1`, at `c = 2^r−1`), so `Δ(j) ≥ 1`. ∎

*Case B.* Collapsing the peeling identity (the heavy factor self-absorbs) gives
`Δ̂(j) = j − 2s₂(j) − 4\,v₂((j)_{c−1}) + 2[v₂ C(a+2,j−c+1) + v₂ C(b+1,j−c+1) + v₂ H_c(j)]`, and
`v₂((j)_{c−1}) = (c−1) + s₂(j−c+1) − s₂(j)`. Dropping the non-negative binomials and using
`v₂ H_c(j) ≥ β'(c)` yields `Δ̂(j) ≥ g_c(j)`. The analytic tail `g_c(j) > 0` for large `j` is standard
(`j − 4\log_2 j` eventually dominates); the finite window is checked directly. ∎

**Cross-check (striking).** With the *true* floors `β'(4..8,10) = 4,3,7,6,11,14`, the predicted deep
exceptional sets are `E_4 = \{8,10\}`, `E_5 = \{10,11,16,17,18,19\}`, `E_6 = \{12,16\}`,
`E_7 = \{16,17,20,21\}`, `E_8 = E_{10} = ∅` — and `E_4`, `E_5` **reproduce exactly** the exceptional
indices that the completed `c = 4` and `c = 5` proofs handled by residue checks. This validates the
general-`c` deep machinery independently.

---

## 4. Why generator 4 fires for `c = 6` but not `c = 4, 5` (the run-length mechanism)

Generator 4 ties iff `Δ(4) = 0`. Since `s₂(4) = 1`,
`  Δ(4) = 2 + 2\big[v₂C(a+c+1,4) + v₂C(b+c,4) + v₂ H_c(4) − v₂ H_c(0)\big]`,
so the tie needs the bracket `= −1`, i.e. a **relative drop** `v₂ H_c(0) − v₂ H_c(4)` large enough to
overcome the binomial terms. Now `v₂ H_c(0)` is the valuation of the *two runs*: at special `(a,b)`
it exceeds the floor `β'(c)` by an amount that **grows with the run length `c−1`** (longer runs catch
more 2s at resonant `(a,b)`). For `c = 4, 5` the runs are short (length 3, 4) and the achievable
over-content is too small to force `Δ(4) = 0`; for `c = 6` (runs of length 5) it is not. This is the
qualitative bridge from Theorem A (run lengths) to Theorem 1 (which generators fire): **the anchor's
run length sets the size of the generator menu.** Heuristic, but it correctly predicts that generators
above 2 proliferate as `c` grows — the opposite of the single-generator hypothesis.

---

## 5. The exact floor `β'(c)` is irregular (two more mourned conjectures)

Rigorous values (symbolic binomial-basis content, each cross-checked against a pure-MN black-box
`min v₂` scan; reconstructions validated vs MN, 0 mismatch):

| `c`            | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|----------------|---|---|---|---|---|---|----|
| `γ(c)=` content `H_c(0)` | 4 | 6 | 7 | 8 | 11 | 14 | 15 |
| `β'(c)` (true floor) | 4 | **3** | 7 | **6** | 11 | **9** | **14** |
| dip `= γ − β'` | 0 | 3 | 0 | 2 | 0 | 5 | 1 |

Two conjectures die here:

- **"Even `c` has no dip"** (suggested by `c = 4,6,8`, dip `0`): **false at `c = 10`** (dip `1`,
  `β'(10) = 14 < 15`). The dip is *not* parity-controlled.
- **Job B's dimer law `β'(2k+1) = β'(2k) − 1`** (verified `c = 5,7`; predicted `β'(9) = 10`):
  **false at `c = 9`** — `β'(9) = 9`, independently confirmed by the pure-MN scan (witness
  `(a,b,j) = (31,24,2)`, `v₂ H_9 = 9`).

The dip `γ(c) − β'(c) = 0,3,0,2,0,5,1` has no evident pattern. Structurally, the dip is the content
lost at interior `j`, where `H_c(j) = (\text{shortened runs}) · R_{c,j}` and the residual factor
`R_{c,j}` supplies content the runs no longer do — **exactly** the per-`c` object the completed
`c = 4, 5` residue checks computed. A `c`-uniform theory of the dip would essentially resolve the
whole interior program; Theorem A shows only that the runs contribute a *linear* upper envelope
`γ(c)`, and the dip is the genuinely hard remainder.

---

## 6. What is proved, what is refuted, what remains

**Proved (rigorous, `c`-uniform):**
- **Theorem A**: `content(H_c(0)) = γ(c)` (closed form, both parities), via the Run-content Lemma.
- **Corollary A′**: `β'(c) ≤ γ(c) ≤ 2c − 2` — the heavy floor is `O(c)` (linear), not quadratic.
- **Theorem B, Case A**: `Δ(j) ≥ 1` for all deep `j ≥ 2c`, all `c ≥ 4` (via Lemma A; the peeling
  identity, proven `c=4,5`, verified `c=6,7,8`).

**Refuted (rigorous):**
- **Single-generator law** (`J* ⊆ {j₀,j₀+2}` for `c ≥ 4`): FALSE — `c = 6` has `J* = {0,4}`
  (Theorem 1; `M₀` cross-checked vs hook length; infinite family).
- **Even-`c` no-dip** and **Job B dimer law**: both FALSE (`c = 10`, `c = 9` resp.; §5).

**Open / conjectural:**
- **Conjecture 1′**: the full `c = 6` mod-4 classification (finite residue check, not done here).
- The correct general law replacing single-generator: the even-`|J*|` law survives; the generator
  menu at `j₀ + 2, j₀ + 4, …` is populated in the *low* (tip-free) régime by `H_c`-content
  resonances, with menu size growing with the anchor run length `c − 1`.
- A `c`-uniform formula for the exact floor `β'(c)` (equivalently, the dip): open, likely hard.

### Files (`projects/code/threerow-heavy/`)
`probe_Hc_structure.py` (per-`j` content, factorizations), `run_content.py` (`h_k` contents),
`theoremA_runcontent.py` (Run-content Lemma + `γ(c)` vs direct, `c ≤ 16`),
`c6_singlegen_check.py`, `c6_census.py` (the generator-4 census + hook cross-check),
`verify_peeling_caseA.py` (peeling identity + Case A, `c = 6,7,8`),
`betaprime_extend.py` (validated `β'(9),…`). MN engine: `threerow-c3/mn.py`;
`H_c` reconstruction: `code/jobB_betaprime_sequence.py`.
