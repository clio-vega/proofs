# The Core Lemma reduces to an all-stay inequality: the min-$D$ Key Lemma

**Clio — 2026-05-22 (prove session, continuation)**

## Summary

This note records a structural reduction of the **Core Lemma** (the minimal-degree
lower bound for $\chi^\lambda(\Omega)$, equivalently $\mathsf{mincost}(\lambda)\ge n(\lambda)$).
The whole lower bound is reduced to a single, clean, **per-walk** statement — the
**Key Lemma** — which in turn bottoms out in the *already proven* all-stay inequality
$D(T)\ge n(\lambda)$.

Concretely:

- **(Proven)** A pairing lemma: in any closed walk the sorting and de-sorting swaps are
  equinumerous. Consequences: the two natural cost conventions (before-state / after-state)
  agree on every closed walk, and
  $$\mathrm{cost} \;=\; \#\{\text{cost-}1\text{ stays}\} + \tfrac12\,\#\{\text{swaps}\}
    \;=\; \#\{\text{descending encounters}\}.$$
- **(Proven, given the Key Lemma)** The reduction
  $$\mathrm{cost(walk)} \;\ge\; \min_t D(U_t)\;\ge\; n(\lambda),$$
  the second inequality being the all-stay lemma (Lemma 1 of
  `2026-05-22-core-lemma-allstay-lower-bound.tex`), which holds for **every** tableau
  regardless of feasibility.
- **(Conjectured, verified $n\le 7$; exhaustively $n\le 6$)** The **Key Lemma**:
  every feasible closed walk satisfies $\mathrm{cost}\ge\min_t D(U_t)$.
- **(Proven structurally + verified)** The Key Lemma splits at the argmin into two claims
  A and B (suffix/prefix), both verified exhaustively for $n\le 6$.

The Key Lemma is the sole remaining gap, and it is much cleaner than the original
Core Lemma: it is a statement about a *single* walk, not a minimum over walks, and it
reduces the swap-laden lower bound to a fact about a *fixed* tableau.

---

## 1. Setup (recap, with the cost convention pinned down)

Notation as in the companion papers. $\lambda\vdash n$, rows/cols $0$-indexed,
$\mathrm{content}(r,c)=c-r$, $\mathrm{cont}_T(v)$ the content of the cell of value $v$.
Axial distance $d_T(i)=\mathrm{cont}_T(i{+}1)-\mathrm{cont}_T(i)$; descent indicator
$\delta_T(i)=[\,d_T(i)<0\,]$. The staircase word is
$\mathbf w_0=(1)(2\,1)\cdots(n{-}1\cdots 1)$, length $L=\binom n2$; letter $i$ occurs $n-i$
times. A closed walk is $U_0\to U_1\to\cdots\to U_L=U_0$, one step per letter $i_k$, each
step a STAY ($U_k=U_{k-1}$) or a SWAP ($U_k=s_{i_k}U_{k-1}$, allowed only when
$|d_{U_{k-1}}(i_k)|\ge 2$). A walk is *feasible* iff no step has a before-state with
$d_{U_{k-1}}(i_k)=-1$ (such a state kills the contribution; equivalently the active pair
is never in the same column when its comparator fires).

**Cost convention (after-state).** The order-at-$q=0$ computation of the Hoefsmit matrix
$M=T_i+1$ shows the cost of a step is read off the **column** index, i.e. the *after*-state:
$$\mathrm{cost}_k=\delta_{U_k}(i_k)=[\,d_{U_k}(i_k)<0\,].$$
Per step: STAY at $d>0$ costs $0$; STAY at $d<0$ costs $1$; SWAP from $d\ge 2$
("de-sorting") costs $1$; SWAP from $d\le-2$ ("sorting") costs $0$.

*(Remark.* The library `walklib.py` charged cost by the *before*-state. Proposition 1
below shows the two agree on every closed walk, so all prior verifications stand; but the
after-state convention is the one dictated by the matrix orders, hence the one proved here.)*

The all-stay cost of a fixed tableau is
$$D(T)=\sum_{k=1}^L \delta_T(i_k)=\sum_{i=1}^{n-1}(n-i)\,\delta_T(i),$$
i.e. cost of the all-stay walk at $T$. **All-stay lemma (proven, cited):**
$D(T)\ge n(\lambda)$ for every $T\in\mathrm{SYT}(\lambda)$, with equality at
$T_{\mathrm{rs}}$.

---

## 2. The pairing lemma and cost reformulations (proven)

For a tableau $U$ let the **content-sequence** be $a_U=(\mathrm{cont}_U(1),\dots,\mathrm{cont}_U(n))$
and let $F(U)=\#\{\,1\le j<j'\le n:\ \mathrm{cont}_U(j)>\mathrm{cont}_U(j')\,\}$ be its number
of inversions.

**Lemma 1 (swap = adjacent transposition).** A SWAP at comparator $i$ exchanges exactly the
two entries $a_U(i),a_U(i{+}1)$ of the content-sequence, leaving all others fixed; hence
$F$ changes by exactly $\pm1$: $+1$ for a de-sorting swap (before-state ascending,
$a_U(i)<a_U(i{+}1)$), $-1$ for a sorting swap (before-state descending). A STAY leaves $F$
unchanged.

*Proof.* $s_i$ relabels values $i\leftrightarrow i{+}1$, i.e. swaps the cells holding $i$
and $i{+}1$; in the content-sequence this swaps positions $i,i{+}1$ and fixes the rest. An
adjacent transposition of a sequence changes its inversion count by $\pm 1$, the sign being
$+1$ iff the swapped pair was in increasing order. (Verified computationally for all SYT of
all shapes $|\lambda|\le 7$.) $\qquad\blacksquare$

**Proposition 1 (pairing).** In every closed walk, $\#\{\text{sorting swaps}\}
=\#\{\text{de-sorting swaps}\}$.

*Proof.* $F(U_0)=F(U_L)$ since $U_0=U_L$. By Lemma 1, $0=\sum_k \Delta F_k =
\#\text{de-sorting}-\#\text{sorting}$. $\qquad\blacksquare$

**Proposition 2 (cost reformulations).** For a closed walk,
$$\mathrm{cost}
=\underbrace{\#\{\text{cost-}1\text{ stays}\}+\#\{\text{de-sorting swaps}\}}_{\text{after-state}}
=\underbrace{\#\{\text{cost-}1\text{ stays}\}+\#\{\text{sorting swaps}\}}_{\text{before-state}}
=\#\{\text{cost-}1\text{ stays}\}+\tfrac12\#\{\text{swaps}\}.$$
Equivalently, writing $\mathrm{cost}=\sum_k\delta_{U_{k-1}}(i_k)$,
$$\boxed{\ \mathrm{cost}=\#\{\,k:\ \text{the active comparator }i_k\text{ is content-descending in the before-state }U_{k-1}\,\}.\ }$$

*Proof.* A cost-$1$ stay has before $=$ after descending. A de-sorting swap has after-state
descending (cost $1$); a sorting swap has before-state descending. So the after-state count
is (cost-$1$ stays)$+$(de-sorting swaps) and the before-state count is (cost-$1$ stays)$+$
(sorting swaps); these are equal by Proposition 1, and their average is
(cost-$1$ stays)$+\tfrac12\#$swaps. The boxed form is the before-state count, listing
exactly the steps whose before-state has $\delta_{U_{k-1}}(i_k)=1$. $\qquad\blacksquare$

---

## 3. The reduction: Key Lemma $\Rightarrow$ Core Lemma (proven)

**Key Lemma (conjectural; verified $n\le 7$).** For every feasible closed walk
$U_0,\dots,U_L=U_0$,
$$\mathrm{cost}\ \ge\ \min_{0\le t\le L} D(U_t).$$

**Theorem (Core Lemma, conditional on the Key Lemma).** For every $\lambda$,
$\mathsf{mincost}(\lambda)\ge n(\lambda)$.

*Proof.* Let $W$ be any feasible closed walk. By the Key Lemma,
$\mathrm{cost}(W)\ge\min_t D(U_t)$. By the all-stay lemma, $D(U_t)\ge n(\lambda)$ for **every**
visited tableau $U_t$ (the all-stay inequality holds for all SYT, with no feasibility
hypothesis). Hence $\min_t D(U_t)\ge n(\lambda)$ and $\mathrm{cost}(W)\ge n(\lambda)$.
Minimising over feasible closed walks gives $\mathsf{mincost}(\lambda)\ge n(\lambda)$.
$\qquad\blacksquare$

This is the entire content of the missing lower bound, now resting on a single per-walk
inequality whose right-hand side is the *already proven* all-stay quantity.

---

## 4. Toward the Key Lemma: the split at the argmin (verified)

Let $t^\*$ be a time minimising $D(U_t)$ and $W:=U_{t^\*}$, so $\min_t D(U_t)=D(W)$.
Split the cost at $t^\*$:
$$\mathrm{cost}=\underbrace{\sum_{k=1}^{t^\*}\delta_{U_k}(i_k)}_{\mathrm{cost}_P}
+\underbrace{\sum_{k=t^\*+1}^{L}\delta_{U_k}(i_k)}_{\mathrm{cost}_Q},
\qquad
D(W)=\underbrace{\sum_{j\le t^\*}\delta_W(i_j)}_{c_{t^\*}(W)}
+\underbrace{\sum_{j>t^\*}\delta_W(i_j)}_{R_{t^\*}(W)}.$$

**Claim B (prefix).** $\mathrm{cost}_P\ \ge\ c_{t^\*}(W)$.
**Claim A (suffix).** $\mathrm{cost}_Q\ \ge\ R_{t^\*}(W)$.

Either claim plus the trivial other half would already give the Key Lemma; together they
give it exactly:
$\mathrm{cost}=\mathrm{cost}_P+\mathrm{cost}_Q\ge c_{t^\*}(W)+R_{t^\*}(W)=D(W)=\min_t D$.

**Status.** Claims A and B are verified **exhaustively** for all feasible closed walks of all
shapes with $n\le 6$ (and the Key Lemma itself for $n\le 7$). Both hold with no exceptions.

**Why they are not yet proven.** The *generic* segment lemma — "if the start (resp. end)
of a walk segment along a contiguous sub-word is $D$-minimal over the segment, then the
segment costs $\ge$ the all-stay cost of that endpoint over the sub-word" — is **false**
(explicit counterexamples already at $(2,1),(2,2),(3,1)$). Therefore A and B are *not*
consequences of a standalone segment statement: their truth uses that $W$ is $D$-minimal over
the **entire** closed walk (globality), and almost certainly the **closure/feasibility**
constraint. The obstruction is sharp:

> At a single swap, $D$ can change by as much as $\pm 3(n-i)$ (it touches comparators
> $i{-}1,i,i{+}1$ with weights $n-i{+}1,n-i,n-i{-}1$), while the step's cost changes by only
> $0$ or $1$. So $D$ is **not** Lipschitz with respect to cost; a single sorting swap can
> drop $D$ far below the cost paid. This is exactly why no naive potential of the form
> $c_k(T)+\psi(T)$ or $n(\lambda)-R_k(T)$ works (see §5), and why the proof must exploit the
> global return constraint.

---

## 5. Negative results (ruling out the easy certificates)

The following candidate node-potentials $\phi_k(T)$ were each tested on all edges of the
layered graph for all feasible shapes $n\le 7$ and **fail**:

1. **$\mathrm{cost}\ge D(\text{start})$** — false; walks genuinely sort below the all-stay
   cost of their own starting tableau (e.g. $(3,2)$ start $124/35$: $D=4$ but a walk costs $3$).
2. **Natural potential $\phi_k(T)=n(\lambda)-R_k(T)$** ($R_k=$ remaining all-stay cost) —
   boundary correct, stay-edges tight, but violated on every $D$-decreasing (sorting) swap;
   the number of violated swap edges grows like $n-3$.
3. **Forward-propagated** potential from the natural boundary — collapses (the all-stay
   "cost-to-go" is wrong because the optimal continuation is *not* all-stay).
4. **$k$-independent** $\phi_k(T)=c_k(T)+\psi(T)$ — the difference-constraint system has
   negative cycles for every nontrivial shape.
5. **Symmetric inversion potential** $\phi_k(T)=c_k(T)+\tfrac12 F(T)$ — boundary and
   stay-edges fine, but swap-edge violations up to $\sim n/2$.
6. **Running-minimum invariant** $\Theta_k:=\mathrm{cost}_{[1,k]}+R_k(U_k)\ \ge\ \min_{t\le k}D(U_t)$
   — would prove the Key Lemma by induction, but is **false** mid-walk (e.g. $(4,2)$ at
   $k=12$). The endpoint statement ($\Theta_L=\mathrm{cost}\ge\min_t D$) is nonetheless true.

A valid **integral** node-potential certifying $\mathsf{mincost}=n(\lambda)$ does exist for
each shape (LP dual, computed), but it is genuinely $k$-dependent and has resisted a
closed-form description; the natural closed form is the value function of the two-parameter
DP on (current tableau, running-min-$D$), which proves the Key Lemma for each fixed $n$ but
not uniformly.

---

## 6. Verification

All checks use the validated Hoefsmit module
`scratch/2026-04-25-biswal-chebyshev/compute_p_lambda.py`.

- **Cost conventions agree** (before $=$ after) and both equal $n(\lambda)$: all feasible
  shapes $n\le 7$ (`2026-05-22-corewalk.py`).
- **Pairing mechanism** (swap $=$ adjacent transposition of content-sequence, $\Delta F=\pm1$):
  all SYT, $|\lambda|\le 7$ (`2026-05-22-corewalk.py`, inline check).
- **Key Lemma** $\mathrm{cost}\ge\min_t D(U_t)$: DP over (state, running-min-$D$) for all
  feasible shapes $n\le 7$; **direct exhaustive enumeration** of all feasible closed walks
  for $n\le 6$ (`2026-05-22-endpoint-AB.py`, `2026-05-22-minD-test.py`, direct enum). Zero
  violations.
- **Claims A and B**: exhaustive over all feasible closed walks, $n\le 6$. Zero violations.
- **Generic segment lemma**: false (counterexamples from $n=3$ up).
- **Negative potential results**: as listed in §5.

A striking empirical fact: feasible closed walks are **very rare** (e.g. $(2,2,2)$ admits a
single one; $(4,2)$ admits $24$). Feasibility (never letting a comparator fire on a
same-column pair) together with closure is extremely rigid. This rigidity is likely the
lever for a proof of the Key Lemma and deserves a direct combinatorial characterisation.

---

## 7. The precise remaining gap

> **Prove:** every feasible closed walk along the staircase word satisfies
> $\mathrm{cost}\ge\min_t D(U_t)$. Equivalently (via the argmin split), prove Claims A and B,
> using the globality of the $D$-minimal state $W=U_{t^\*}$ and the closure/feasibility
> constraint.

Promising angles, in order of appeal:
1. **Characterise feasible closed walks** directly (they are rare and rigid) and read off
   the bound. The forbidden-state ($d=-1$) avoidance is the binding structure and is
   currently unused in A/B.
2. **Two-parameter potential**: find a closed form for the DP value function
   $g_k(U,m)=\min\{\mathrm{cost}_{[1,k]} : U_k=U,\ \min_{t\le k}D(U_t)=m\}$ and show
   $g_L(U_0,m)\ge m$.
3. **Straightening / induction on swaps** using Proposition 1 (sorting/de-sorting pairing):
   excise a matched swap pair and track the change in both cost and $\min_t D$.
