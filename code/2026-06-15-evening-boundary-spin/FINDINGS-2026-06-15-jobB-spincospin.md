# FINDINGS — Job B (2026-06-15 evening): spin/cospin home — PROBE-FIRST, FAIL

**Verdict: FAIL the (3,2,1) probe. Prune cleanly. Do not scaffold.** The spin/cospin
generating function of domino tableaux **cannot** be a home for
`ord_{q=−1} G_λ = ⌊|2-core(λ)|/2⌋`, because spin reads the **2-quotient** (the dominoes
stacked above the core) while the order law reads the **2-core** (which by definition
carries *zero* dominoes). They read complementary halves of the abacus. The candidate
is **identically trivial on every 2-core** — i.e. bounded `≡ 0` on the exact shapes
(the staircases) where the law must be unbounded. This is the same failure mode that
killed the fusion route, now pinned to a one-line structural reason.

**The dream-cycle premise was backwards.** CODE.md called spin/cospin "the only object
that *reads the 2-core*." The probe shows the opposite: spin/cospin **does not read the
2-core at all** — it is a statistic on the 2-quotient, orthogonal to `|2-core|`.

Script: `jobB_spin_probe.py` (self-contained; standard domino tableaux built up from
the 2-core via `core_quotient.py`; spin = #vertical dominoes; machinery sanity-checked
— (2,2) → spin_gen `{0:1,2:1}` = two-horizontal vs two-vertical tilings).

---

## The candidate

For `λ` with 2-core `C`, a **standard domino tableau** is a chain
`C = μ₀ ⊂ μ₁ ⊂ … ⊂ μ_k = λ`, each `μ_i/μ_{i−1}` a single domino (2 adjacent cells);
`spin(D) = #vertical dominoes`, `S_λ(q) = Σ_D q^{spin(D)}`. (Carré–Leclerc /
Stanton–White / Schilling–Shimozono–White — the cospin/ribbon-LR object, the named hub.)
Three natural readings as `candidate(λ)`: top degree `deg = max spin`, order of vanishing
`ord_{q=−1} S_λ`, and value `S_λ(−1)`.

## The probe table (candidate vs. KNOWN order law)

| shape | order law `⌊\|2-core\|/2⌋` | cand.deg | cand.ord@−1 | cand.val@−1 | #dominoes | `S_λ(q)` |
|---|---|---|---|---|---|---|
| δ₃=(3,2,1) | **3** | 0 | 0 | 1 | 0 | `{0:1}` |
| δ₄=(4,3,2,1) | **5** | 0 | 0 | 1 | 0 | `{0:1}` |
| δ₅=(5,4,3,2,1) | **7** | 0 | 0 | 1 | 0 | `{0:1}` |
| δ₆ | **10** | 0 | 0 | 1 | 0 | `{0:1}` |
| (2,1) [a 2-core] | **1** | 0 | 0 | 1 | 0 | `{0:1}` |
| (2,2) [bounded] | 0 | 2 | 0 | 2 | 2 | `{0:1,2:1}` |
| (3,1) [bounded] | 0 | 1 | 0 | −1 | 2 | `{1:1}` |
| (4,4) [bounded] | 0 | 4 | 0 | 6 | 4 | `{0:2,2:3,4:1}` |

**Bar (a) — reads the 2-core? NO.** **Bar (b) — unbounded? NO** (≡0 on staircases).
Both bars failed. On the bounded shapes the candidate is *nonzero and growing* while the
order law is `0` — anti-correlated with the target.

## Why — the decisive discrimination

**(b) Different 2-core, same (empty) 2-quotient ⇒ order law differs, spin identical.**
Every 2-core has empty 2-quotient and `S_λ(q) ≡ 1`:

| λ (a 2-core) | 2-core | `⌊\|core\|/2⌋` | 2-quotient | `S_λ` |
|---|---|---|---|---|
| (1) | (1) | 0 | (∅,∅) | `{0:1}` |
| (2,1) | (2,1) | 1 | (∅,∅) | `{0:1}` |
| (3,2,1) | (3,2,1) | 3 | (∅,∅) | `{0:1}` |
| (4,3,2,1) | (4,3,2,1) | 5 | (∅,∅) | `{0:1}` |

The order law takes 0,1,3,5,… but the candidate is the *constant* 1. A statistic that
cannot separate these shapes cannot equal a quantity that does.

**(a) Same (empty) 2-core, different 2-quotient ⇒ order law constant (0), spin varies:**
`(4,)→deg 0`, `(2,2)→deg 2`, `(3,3)→deg 3`, `(4,4)→deg 4`. Spin varies exactly where the
law is constant. So `S_λ` is a faithful invariant of the 2-quotient and *blind* to the
2-core — the precise orthogonal complement of what the order law needs.

## Structural reason (the prune, in one line)

`#dominoes(λ) = (|λ| − |2-core|)/2 = Σ_r |λ^{(r)}|` is a 2-**quotient** quantity; `spin`
is a statistic on those dominoes, so `S_λ` factors entirely through the 2-quotient and
is `≡ 1` whenever the quotient is empty — i.e. on every 2-core. But `⌊|2-core|/2⌋` is a
2-**core** quantity, maximal (and unbounded) exactly on the 2-cores where `S_λ ≡ 1`.
No spin/cospin generating function — degree, value, or order at `q=−1` — can recover it.
**Spin/cospin is dead as a home for this order law.** Schilling–Shimozono–White is not
worth opening for this; the obstruction is upstream of any branching identity.

## Status of the survivor list

All four external routes were already dead (odd-content, CTT-domino, Sawin-adjacency,
Dobner fusion). Spin/cospin was the last dream-crowned survivor; this probe kills it.
**No external CSP/algebraic home for `ord_{q=−1}G_λ = ⌊|2-core|/2⌋` remains on the
board.** The law's identified home stays the **Littlewood `t=2`** instance
(Albion 2212.07343; factorisation, not a sieve) — i.e. the grade is genuinely *ours*,
not an off-the-shelf cyclic-sieving phenomenon. Recommend the dream cycle stop seeking a
domino/quotient CSP and instead lean on the Littlewood-factorisation framing.

---

## Secondary — MO#509068 backbone (reached; partial, honest scope)

The reusable, hard part is **verified**: for the hook `λ=(a,1^{2m−a})⊢2m`,
> `M_j := ⟨s_λ, e₂^j h₁^{2(m−j)}⟩ / C(m,j) = C(2m−1−j, a−1)`  (j=0..m),

confirmed against Murnaghan–Nakayama, **1366 cases, 0 mismatches**
(`jobB_hook_mo509068.py`). As a function of `j` this is a ratio of falling factorials
`(2m−1−j)!/[(a−1)!(2m−a−j)!]` — the Kummer/telescoping handle the answer will use.

**What I could NOT verify (and why I'm not fabricating it).** MO#509068's exact prose —
*which* signed sum vanishes, in which ring / at which specialisation, with "j=0 and j=m
survive, interior vanish" — is **not in this session's files**. I tested the natural
2-adic reading `val(j)=j+2(v₂C(m,j)+v₂M_j)`: it gives `J*={0}` or `{0,2}`, **not** the
endpoints-only pattern — so that is the *wrong* sum. Drafting a public MO answer around
a guessed statement would be irresponsible. **Recommendation:** finish MO#509068 in a
**write/prove session with the MO question open**; the closed form above is the solid
foundation it rests on, now machine-verified and ready.

## Files
`jobB_spin_probe.py` (probe + both discriminations), `core_quotient.py` (2-core/quotient),
`jobB_hook_mo509068.py` (verified hook closed form `M_j=C(2m−1−j,a−1)`).
