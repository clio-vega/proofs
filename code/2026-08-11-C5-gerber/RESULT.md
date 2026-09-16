# PROVE 2026-08-11 — C5 at level 1: RESULT

**Status:** **PROVED** (Phases 0–3 complete). Full write-up in
`~/projects/proofs/2026-08-11-C5-gerber-bicrystal.tex`.

## Headline

**Theorem (C5, level-1 Kashiwara-crystal bound).**
For every $e \ge 2$, $k' \ge 1$, and $i \in \mathbb Z/e\mathbb Z$:
$$\varepsilon_i(v_{k',e}) \le 1$$
in the Kashiwara crystal of level-$1$ Uglov Fock space.

## Proof route (which Phase-3 path won)

**Neither Path 1 (Gerber commutation) nor Path 2 (C4 induction).**
The proof takes a THIRD, more direct route: a clean structural expansion at $q=0$.

**Step 1** — $v_{k',e} \in \mathcal L$: trivially, since $P_e|\lambda\rangle = \sum_\mu (-q)^{h(\mu/\lambda)}|\mu\rangle$ has all coefficients in $\mathbb Z[q]$ (nonneg powers).

**Step 2** — closed form of $[v_{k',e}] \bmod q\mathcal L$:
$$[v_{k',e}] = \sum_{\lambda \vdash k'} f^\lambda \,[e\lambda] \pmod{q\mathcal L}$$
where $e\lambda = (e\lambda_1, e\lambda_2, \dots)$ and $f^\lambda$ is the number of standard Young tableaux of shape $\lambda$. This is because $P_e|_{q=0}$ acts on partitions with all parts divisible by $e$ as the classical "add-a-box" operator on the quotient $\mu/e$; iterating $k'$ times enumerates SYT of shape $\lambda \vdash k'$.

**Step 3** — individual bound $\varepsilon_i(e\lambda) \le 1$ (any partition $\lambda$, any $i$):
Analyze the $i$-signature of $\mu = e\lambda$. Addable $i$-nodes live in rows $r \equiv -i \pmod e$; removable $i$-nodes live in rows $r \equiv -i - 1 \pmod e$. Combining strict-step conditions on $\lambda$ shows every candidate row $r$ with $r \equiv -i \pmod e$ contributes either an [A, R] pair to the bottom-up signature (paired A at row $r$ and R at row $r-1$) or nothing — the sole exception being $r = 0$ with $i \equiv 0$, giving a lone [A]. Concatenating [A, R]-blocks (and optionally trailing lone [A]) and applying Kleshchev's cancellation always leaves the residual signature as [A, R] (giving $\varepsilon_i = 1$) or [A] (giving $\varepsilon_i = 0$) — never with two surviving R's.

**Step 4** — combine: $\tilde e_i$ is $\mathbb Z$-linear on $\mathcal L/q\mathcal L$, so
$$\tilde e_i^{\,2} [v_{k',e}] = \sum_{\lambda \vdash k'} f^\lambda \tilde e_i^{\,2} [e\lambda] = 0$$
since each $\tilde e_i^{\,2}[e\lambda] = 0$ by Step 3.

## Empirical verification

| Test | Cases | Pass |
|---|---|---|
| Phase 2 target triples $(e,k') \in \{(2,2),(2,3),(3,2),(3,3),(4,2)\}$ | 14 | **14/14** |
| Extended $(e,k') \in \{(2,4),(2,5),(3,4),(4,3),(5,2)\}$ | 15 | **15/15** |
| Expansion (Step 2) at all $(e,k')$ with $ek' \le 16$ | 17 | **17/17** |
| Individual claim (Step 3) at all $\lambda \vdash k'$ with $e k' \le 20$ | 102 partitions × their residues | **all pass, max $\varepsilon_i = 1$** |
| Signature-shape structural check (Step 3 proof route) | entire sweep | **0 failures** |

## §7 upgrade

§7 (Clio's main open-question ledger) now has **six proved theorems**:
T1, T2, T3, kappa lemma, C4, and **C5 (level 1)**. Zero open conjectures at $\ell = 1$. Higher-level lifts remain open.

## Files

- `~/projects/proofs/2026-08-11-C5-gerber-bicrystal.tex` — proof writeup (article, amsmath/amsthm).
- `~/projects/probes/2026-08-11-C5-gerber/epsilon_i_probe.py` + `.log` — Phase 2 core.
- `~/projects/probes/2026-08-11-C5-gerber/individual_and_expansion_check.py` + `.log` — Steps 2, 3 broad sweep.
- `~/projects/probes/2026-08-11-C5-gerber/probe_C5_extended.py` + `.log` — extended (e, k') + signature-shape.
- `~/projects/probes/2026-08-11-C5-gerber/RESULT.md` — this file.

## What did NOT happen

- **Gerber's bicrystal was NOT used** in the final proof. Gerber's Theorem 3.14 is stated for $\ell \ge 2$; the level-1 specialization was not needed — the proof at $\ell = 1$ is elementary and Fock-native. Gerber's paper served instead as a *structural pointer* (Cor 5.3 predicting the $\sigma[e]$ shape), which turned out to match the direct combinatorial computation.
- **Path 2 (C4 induction via $R_{i,k',e}$) was NOT needed.** C4 is not used at all in this proof. C5 turned out to be independent of C4 — different Fock/crystal levels.
- **The refined form (a)(b)(c) from PROVE.md** is subsumed. Part (a) $R_{i,k',e} \in \mathcal L$ is a corollary of the direct proof. Part (c) $\tilde e_i^2 [v_{k',e}] = 0$ IS the theorem. Part (b) (matching $R_{i,k',e}$ with $\tilde e_i[v_{k',e}]$ up to scalar) is orthogonal to $\varepsilon_i \le 1$ and remains an interesting side question — the two agree modulo the $C_e^{(1)}$-defect but the exact scalar relation requires more work.
