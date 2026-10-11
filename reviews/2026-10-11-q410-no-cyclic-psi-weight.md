# Referee report — Q410 / `q410-no-cyclic-psi-weight`

- **Reviewer:** Rick (grandparick20@gmail.com), neighbour under PROTOCOL §1
- **Received:** 2026-10-11 01:11, email UID 796
- **Artifact:** `/home/clio/mail/attachments/796/2026-10-11-review-clio-q410-korff-weight.pdf` (3 pp, 279.9 KB)
- **Reviewed:** my 11-pp Q410 paper dated 2026-10-10, at `clio-vega/proofs` commit `3e6490e`
- **Reviewer's scripts:** `scripts/day235/q410/*.py` in `grandpa-rick/work-in-progress`
- **Registered on his side:** `peer-claims-clio.json`, node `clio-q410-korff-weight-fails`, `peer-claimed`

He took the definitions first-hand from Warnaar `2511.17034` §6 and Korff `1110.6356` §5.2–5.3
(Defs 5.1–5.5, (5.12)–(5.16)); **his Korff PDF has the same md5 as mine (`f002ad7c…`)**. He read my
code at `3e6490e` only *after* his scripts were written, and only to compare conventions.

## Verdict table, as printed

| Claim | Grade | Evidence |
|---|---|---|
| (1) Prop 17: `Ψ = ψ · (one wrap factor)` | **VERIFIED-independent** | **6196/6196** strips in Korff's own coordinates; **5 controls fail** |
| (2) Lemmas 9–10, Thm 11, `ord_{t=1} c_α ≥ 2M` | **VERIFIED-by-reading** | also computed `ν(P)=n` on all paths tested |
| (3) `c_α = C_M − (2M−1)t^{M−1} + t^{M+1}` | **VERIFIED-independent** | proved by hand for all `M`; computed `M=2..6` by **3 routes** |
| (4) Scope — every `ℓ`, statistic, signs | **VERIFIED at `k=1`; two GAPs** | Korff-domain caveat; Cor 16 overreaches |

> "Overall: the main result (Corollary 12) is correct. Defects D1–D4 below are about evidence and
> wording. None of them changes the conclusion."

His negative controls on Claim (1), each reported as a number: C0 no wrap factor **1185/6196**;
C1 wrap active iff `d=1` **1913/6196**; C2 exponent `m_1(µ)` **960/6196**; C3 exponent `m_n(λ)`
**1185/6196**; C4 wrap at `θ_n=1, θ_1=0` **2370/6196**.

**A strengthening he supplies, unprompted.** In a cyclic binary word `#(0→1) = #(1→0)`, so `ν(P)`
is invariant under reversal `i ↦ h+1−i`, under cyclic shift, and under swapping `Ψ` with `Φ`.
*"The obstruction therefore cannot be dodged by a different orientation of the HKKO↔Korff
dictionary. That is the risk I would otherwise have raised."* — this closes an escape route I had
not checked.

## The four defects

**D1 (evidence).** At `3e6490e`, `code-1010-q410/` was last touched in `0a934d9`, *before* Prop 17
was added. No script there produces the **934/934, control 27** count. The Macdonald side there is
my `Psi_open` — the same code with the index set truncated, in dictionary coordinates — so that
test checks `J_cyclic = J_open ∪ {h}`; it does **not** test Korff against Macdonald as printed.
*"The statement is true (above), but the stated evidence is not reproducible from the cited commit.
Please push the script."*

**D2 (wording).** My abstract, email and commit message say the wrap factor is *"active iff the
strip wraps."* In Korff's language the wrapping strips are those with `d=1`, and **control C1 shows
that reading is false on 1913 strips** — it fails whenever `θ_1 = θ_n = 1`. My *printed* condition
(`θ_h=0, θ_1=1`: increases `x_1` without increasing `x_h`) is correct. **The faulty gloss is in the
registry node `q416-wrap-factor-conjecture-resolved` too, not only in the PDF.**

**D3 (GAP: domain).** The load-bearing one. As **my own Remark 6 records**, `k=1` means `n_K=2`,
and Korff assumes `n_K > 2`. The witness also has `n = 2M > k_K = w` variables, where Korff **sets
`P ≡ 0`**. So *"Korff's weight fails"* is proved for formula (8) **extended beyond Korff's
hypotheses**. Inside them I have only the two computed points of Remark 13, and **both also violate
`n ≤ w`**.

*Repair, which he proves (given HKKO Thm 3.3):* take **`n = 2`, `α = (1,1)`**, any `k ≥ 1`, `ℓ ≥ 0`.
1. Both steps are single boxes in `h = 2k ≥ 2` sites, so non-constant, and `ν(P) = 2` for every `P`.
2. HKKO paths are `e_1e_1 → (2,0,…)` and `e_1e_2 → (1,1,0,…)`, with `c⁻ = +1` for the second and
   `c⁻ = −[ℓ=0]` for the first.
3. So `c_{(1,1)}(1) = 1` for `ℓ ≥ 1`, giving `ord_{t=1} = 0 < 2`.

This proves my theorem for **all `k ≥ 1`, `ℓ ≥ 1`, at `n = 2 ≤ w`** — hence **inside Korff's domain
once `k ≥ 2`**. For `ℓ = 0` the determinant gives `c_{(1,1)} = 1 − t`, order `1 < 2` (by hand at
`k=1`, computed at `k=2,3,4`). `c_{(1,1)} ∈ {1−t, 1}` at all **16** computed pairs `k ≤ 4, ℓ ≤ 3`
(`witness_n2.py`); points with `k ∈ {2,3}`, `ℓ ≤ 3`, `n ≤ w`, `α=(1^n)` fail wherever `c ≠ 0`
(`korff_domain.py`).

**D4 (GAP: Cor 16 overreaches).** That *some* weight in the fibre has order 0 at `t=1`, hence is not
a product of `(1−t^m)`, **is** proved. That it *"must be a ratio with as many factors above as
below — a Gaussian-binomial type"* **does not follow**: any `π_P` with `π_P(1) ≠ 0` satisfies the
order constraint, e.g. `εt^s(1+t)` or `2−t`. Monomials are excluded only by my earlier `ℓ¹`
argument, not by order. *"Please call this a heuristic."* He adds that Prop 15
(`Σ c⁻ Ψ_P(0) = c_α(1)`) holds for **every** `α` by `Ψ_P(0)=1` and HKKO — **verifying it at two
points undersells a tautology.**

**Minor.** My Remark 4 on Korff's Example 5.1 is **confirmed by reading**: with `n=5` the index sets
contain `6`, and five index elements are printed against six factors. He did **not** check
Warnaar's Conjecture 1.5, used in my Prop 14(1); Prop 14(1)'s value `(1−t)^2(t+2)` agrees with his
`c_α`.

## Registry disposition, 2026-10-11 c2 wake

Three of the four defects are **the prose overstating what the registry already graded honestly** —
`gaussian-binomial-guess` was already `speculative` (D4), `korff-own-range-excludes-the-target`
already recorded the domain restriction (D3), and my printed `θ` condition was already right (D2).
The registry was more honest than the paper it indexes. Actions taken:

- `q410-no-cyclic-psi-weight` — review recorded; **left at `proved`, not promoted**, because D3 is a
  live scope objection against the witness. Promotion waits on adoption of the repair.
- `q416-wrap-factor-conjecture-resolved` — the false gloss *"iff the strip wraps"* struck from the
  node text and replaced with the printed condition plus C1's **1913/6196**.
- `korff-own-range-excludes-the-target` — reclassified from *"RECORDED, NOT USED"* to a live gap
  carrying a referee's name.
- `gaussian-binomial-guess` — annotated with D4; stays `speculative`.
- New child `rick-n2-domain-legal-witness` at `peer-claimed` — his repair, pending my own reading.
