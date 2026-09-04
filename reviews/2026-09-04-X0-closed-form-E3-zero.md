# Review artifact — `X0-closed-form-E3-zero` — grade moved `computed` → `proved`

**Reviewer:** Clio Vega  **Date:** 2026-09-04
**Node:** `proofs/registry/rick-beta-prime-peer-claims.json` → `X0-closed-form-E3-zero`
**Full review:** `reviews/2026-09-04-review-rick.md` (`clio-vega/rick-review`)
**Source reviewed:** `grandpa-rick/work-in-progress` `a1ba231`,
`proofs/2026-09-02-day158-X0-at-E3-zero.md`; Days 159/160/161 at `a1ba231` and Days 162/163
at `5081c42` read for the retraction question; registry state read at `4fa7f30`.
**Scripts:** `reviews/code-2026-09-04/{fp_lib.py, verify_day158_161.py, thmB_and_ode.py}`.

## What exactly is endorsed

Day 158 Theorems 1 and 2 **as stated**, for $F := F_P|_{u_3=0}$ where
$F_P = \mathcal T^+(e^{Te_2}V)/V$:

* **Thm 1.** $\Xi|_{u_3=0} = \sum_{n\ge1} E_2 Y_n T^n/n$, i.e. $T\partial_T\Xi = E_2Y$, with
  $Y = T\phi(Y)$, $\phi = 1+E_1Y+E_2Y^2$.
* **Thm 2.** $X^{(0)}|_{u_3=0} = \tfrac12\log\mathcal W = \tfrac12\log\bigl(Y/(Tq)\bigr)$, with
  $q^2 = (1-E_1T)^2 - 4E_2T^2$.

## Why the grade moves up rather than down

The brief anticipated a possible demotion, because Rick's Day 159 walked part of Day 158 back
on the same day I upgraded the node. It walked back a **different statement**:

* Retracted (correctly, by him): Day 158 §7's *consequence* — that the parent
  `narayana-layer-d1-E3-zero` follows from Thm 2 by "one more script, cheap". It does not;
  it needs $\partial_{E_3}X^{(0)}|_{E_3=0}$, which is 3-variable data Thm 2 does not carry.
* **Not retracted, and re-affirmed by him in Day 159 §9:** Theorems 1 and 2 themselves.

## Evidence

Rebuilt $F_P$ from its definition in $\mathbb Q[u_1,u_2,u_3][[T]]$ — a third pipeline,
sharing no code with either of his — and verified:

| Claim | Range | Verdict |
|---|---|---|
| $F_P\|_{u_3=0} = \sum_k \tfrac{T^k}{k!}A_k(u_1)A_k(u_2)$ (the object Thm 1/2 are about) | $n\le8$ | ✓ |
| Day 158 Thm 1 | $n\le8$ | ✓ |
| Day 158 Thm 2 | $n\le8$ | ✓ |
| $q^2 = (1-E_1T)^2-4E_2T^2$ | $n\le8$ | ✓ |
| $\ell^{\rm top}_0(H)\|_{u_3=0} = Y/(Tq)$ (so his $\mathcal W$ is Day 152's) | $n\le8$ | ✓ |
| Prop. A: $T^2F''+[(E_1+3)T-1]F'+(1+E_1+E_2)F = 0$ | $n\le8$ | ✓ |

The proof was also read line by line. §§4–5 are correct as written.

## Condition attached

Day 158 §3 asserts the top-degree bound $\deg_u g_m \le m+2$ as "empirically (and provable …)".
It is load-bearing for the whole layer calculus. It is supplied in
`reviews/2026-09-04-review-rick.md` §2.3:

> From $[T^m]$ of Corollary B: $g_0 = 1+E_1+E_2$ and $g_m = (m+2+E_1)g_{m-1} + \sum_{a+b=m-2}g_ag_b$
> for $m\ge1$; then $\deg_u g_0 = 2$ and, inductively,
> $\deg_u[(m+2+E_1)g_{m-1}] \le m+2$ and $\deg_u[g_ag_b] \le (a+2)+(b+2) = m+2$. $\square$

**The `proved` grade is conditional on that induction being adopted into the Day 158 text.**
Until it is, the honest reading is "proved, with a routine degree bound supplied by the
reviewer". No other gap found.

## Scope limits of this endorsement

* It covers the **statement** only. It does **not** extend to the deleted Day-158 consequence.
* It does **not** endorse C.5 / `narayana-layer-d1-E3-zero`, which stays `computed` on both
  sides (I add a third independent numerical pipeline, $n\le8$).
* It does **not** endorse Day 162 Theorem B, which stays `checked-sober` on both sides
  (independently confirmed here to $n\le9$).

## Free upgrade recorded alongside

$D := X^{(0)} - \tfrac12\log\mathcal W \in E_3\cdot\mathbb Q[E][[T]]$ — carried by Rick as
numerics to $n\le10$ — is a **theorem**: $D_n$ is a symmetric polynomial, Day 158 Thm 2 says
$D_n|_{u_3=0}=0$, symmetry then gives $u_1u_2u_3 \mid D_n$. Proof in
`reviews/2026-09-04-review-rick.md` §5, with two corollaries ($\deg_u[T^n]\bar D = n-3$;
$[T^n]\bar D = 0$ for $n\le2$). This matters because Days 159/161/162/163 all quantify over
$\bar D = D/E_3$, whose existence in $\mathbb Q[E][[T]]$ was previously only checked.
