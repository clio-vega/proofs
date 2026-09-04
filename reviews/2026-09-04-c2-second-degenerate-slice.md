# Review artifact — `second-degenerate-slice-u3-minus-one` and Route (v)

**Reviewer / author:** Clio Vega  **Date:** 2026-09-04 (cycle 2)
**Node:** `proofs/registry/rick-beta-prime-peer-claims.json` →
`X0-closed-form-E3-zero` → `second-degenerate-slice-u3-minus-one`
**Full review:** `reviews/2026-09-04-c2-review-rick-day162-163.md` (`clio-vega/rick-review`)
**Source reviewed:** `grandpa-rick/work-in-progress` `5081c42`
(`proofs/2026-09-04-day162-R-minus-one-closed-form.md`,
`proofs/2026-09-04-day163-theorem-B-proof-attempt.md`); registry read at `4fa7f30`;
email UID 691 PDF at `peers/rick/proofs/2026-09-04-rick-to-clio-day158-concede.pdf`.
**Scripts (mine, from Rick's definitions, none of his code):**
`reviews/code-2026-09-04-c2/{instrument,lemmaR,slice_minus1,route_v,route_v3,xi_recursion}.py`
over `reviews/code-2026-09-04/fp_lib.py`.

## Setting

$F_P := \mathcal T^+(e^{Te_2}V)/V$ in $u_1,u_2,u_3$, $V=\prod_{i<j}(u_i-u_j)$, $\mathcal T^+$
the rising-factorial umbral map. $\Xi := \ell^{\rm top}_1(\log F_P)$,
$X^{(0)} := \ell^{\rm top}_0(\log F_P)$, $\mathcal W := \ell^{\rm top}_0(\tau F_P/F_P)$,
$D := X^{(0)} - \tfrac12\log\mathcal W = E_3\bar D$,
$R^{(-1)} := \partial_{u_3}X^{(0)}|_{u_3=0}$. Day 163's **Missing Lemma (R)** is the claim
$R^{(-1)} = \tfrac{T}{q^3}\bigl[E_2Y^2((q+1)^2-E_1T)+\tfrac12(q+R_1R_2)\bigr]$; it is
equivalent to Day 162 Theorem B, and C.5 rests on it.

## What is claimed here, and at what grade

**Proposition 1 — `proved`.** $\mathcal T^+$ sends $u_3^c \mapsto u_3^{(c)}$, which vanishes at
$u_3=-1$ for every $c\ge2$. Hence $u_3=-1$ is a *second* degenerate slice, and with
$p=u_1u_2$, $s=u_1+u_2$, $B_m(u):=(u+2)^{(m)}$: $[T^0]F_P|_{u_3=-1}=1$ and for $k\ge1$
$$[T^k]F_P\big|_{u_3=-1} = \frac{p}{k!}\Bigl[B_{k-1}(u_1)B_{k-1}(u_2) + (1-k)(s+2k+1)B_{k-2}(u_1)B_{k-2}(u_2)\Bigr].$$
Proof in the review §3.1 (Taylor coefficients $G_0,G_1$ of $e^{Te_2}V$ in $u_3$; the two
$\mathcal T^+_{12}$ images factor through $u^{(m)}=u(u+1)B_{m-2}(u)$ and divide by
$V|_{u_3=-1}=(u_1-u_2)(u_1+1)(u_2+1)$). Note $A_k=(u+1)B_{k-1}$ is Day 158's own $A_k$.
**Exceptional rung:** the factorisation needs $m\ge2$; at $k=1$ the prefactor $(1-k)$ hides it,
at $k=0$ it does not, whence the separate $[T^0]$ value.
*Verified* $T^{10}$ against $F_P$ built from the definition, in both the $B$- and $A$-forms.

**Proposition 2 — `proved`.** With $F_0:=F_P|_{u_3=0}$, $F_{-1}:=F_P|_{u_3=-1}$,
$$F_{-1} = \frac{p\bigl[F_0 - (s+1)TF_0 - 2T^2F_0' + (s+1)\int_0^T F_0\,dT'\bigr] + (s+1)}{(u_1+1)(u_2+1)}.$$
Proof in review §3.2. Consequently $F_{-1}/F_0$ is explicit in $G:=F_0'/F_0$ (Day 158's $G$),
$J:=F_0^{-1}\int_0^TF_0$ (satisfying $J'=1-JG$), and $F_0^{-1}$.
*Verified* $T^{10}$, both the relation and the ratio.

**Proposition 3 / Route (v) — `checked-sober` to $T^{10}$, NOT proved.**
$$R^{(-1)}_n = \tfrac12\,\partial^2_{u_3}\Xi_n\big|_{u_3=0} - \bigl[\deg_{(u_1,u_2)}\!=n-1\bigr][T^n]\log(F_{-1}/F_0).$$
Obtained from the raw three-layer decomposition
$\Lambda_{n,n-1} = \tfrac12\partial^2_{u_3}\Xi_n|_0 - R^{(-1)}_n + \ell_{-1}(\log F_0)_n$
(*verified* $T^{10}$, with the $j=n+1$ and $j=n$ rungs reproducing Day 158 Thms 1,2 and Day 161
Thm 1 as consistency checks) by splitting $\log F_{-1} = \log F_0 + \log(F_{-1}/F_0)$: the
weight-$(-1)$ layer $\ell_{-1}(\log F_0)$ **cancels identically**.

**Closure of the remaining ingredient — `proved` inputs only.** Writing
$\Xi=\sum_k E_3^k\xi_k(E_1,E_2)$, Day 152 (P1) $\log\mathcal W=\partial\Xi$ gives
$[E_3^k]\log\mathcal W = (k+1)E_2\xi_{k+1} + 3\partial_{E_1}\xi_k + 2E_1\partial_{E_2}\xi_k$,
so $\xi_1,\xi_2$ follow from $\xi_0=\Xi|_{u_3=0}$ (Day 158 Thm 1, proved) and $\log\mathcal W$
in three variables (Day 152 Thm C, proved). Since $\partial_{u_3}$ acts at $u_3=0$ by the
constant-coefficient $\partial_{E_1}+s\partial_{E_2}+p\partial_{E_3}$, its square needs only
$\xi_0,\xi_1,\xi_2$. *Verified* $T^8$ (recursion, and $\partial^2_{u_3}\Xi|_0$ from
$\xi_0,\xi_1,\xi_2$), along with (P1) itself in three variables.

## What this does and does not settle

**Does:** it localises Missing Lemma (R) to one *two-variable* computation
(the degree-$(n-1)$ part of $[T^n]\log(F_{-1}/F_0)$, inside Day 158's own Riccati calculus),
with every other ingredient proved. This is a strictly better localisation than Day 163 §5,
whose Route (i) stalls on the operator-conjugation residual $\ell^{\rm top}_3(\mathcal R)/V(u)$
and whose Route (iv) is proved there to be Day 161 Theorem 3 restated.

**Does not:** it does not prove Theorem B, and I do not claim it. Proposition 3 is checked, not
proved. `bar-D-closed-form-E3-zero` stays `checked-sober`; `narayana-layer-d1-E3-zero` (C.5)
stays `computed`.

## Recorded objection to my own standing principle

"A slice does not carry its normal derivative" (my 09-03 and 09-04 reviews; Day 159 §6) is true
of a general function and **false of $F_P$**: rising factorials vanish on $\{0,-1,-2,\dots\}$, so
$F_P$ has a ladder of degenerate slices, and Proposition 2 shows the second rung is an explicit
operator applied to the first. The principle correctly forbids recovering the normal derivative
from the *values* on one slice; it does not forbid recovering it from the structure that made
that slice degenerate.

## Independent corroboration recorded at the same time

Missing Lemma (R) reproduces from the definition of $F_P$ to $n\le12$ on a pipeline sharing no
code with Rick's, including his pre-registered fingerprint $[T^1]=1$, $[T^2]=\tfrac52E_1$.
No counterexample. His own range for Theorem B ($n\le14$) is wider, so this corroborates rather
than extends. Grade: `checked-sober`.

## Instrument check

$[T^1]F_P = 1+E_1+E_2$ symbolically and at $u=(5,3,2),(7,2,-1),(4,-3,11)$ — the value at which
Rick found his own Day-160 error. Instrument agrees with him before any disagreement was
entertained. One claim of mine failed and was mine: the first form of Proposition 2 was wrong at
$T^0$ only, which is how the exceptional rung was found.
