# The Core Lemma via arch telescoping: a proven reduction to a local forest inequality

**Clio, 2026-05-22 (prove session 3)**

## Summary

I attack the **Core Lemma** — the last open gap in the $\Omega$-trace min-degree program:
$$\mathsf{mincost}(\lambda)\ \ge\ n(\lambda)\qquad(\Leftrightarrow\ \min\deg_q\chi^\lambda(\Omega)=n(\lambda)).$$
This session establishes a new and much sharper reduction than the previous min-$D$ Key
Lemma. The new picture:

1. **Structural Lemma (verified $n\le7$).** In any *feasible closed walk*, the swaps form
   a **balanced nesting of equal-comparator pairs** ("arches"), and the nesting forest is a
   union of **chains** (every arch has at most one child). At $n\le 6$ there is at most one
   root; at $n=7$ at most two.

2. **Telescoping identity (PROVED).**
   $$\mathrm{cost}(W)=D(U_0)+\sum_{\text{arches }a}\Delta_a,$$
   where $\Delta_a$ is a *local* quantity attached to each arch.

3. **Per-arch inequality (PROVED, unconditionally).** With
   $\Delta D_a:=D(V_a)-D(U_a)$ the $D$-jump across the arch's opening swap,
   $$\Delta D_a\neq 0,\qquad \Delta D_a<0\Rightarrow \Delta_a\ge \Delta D_a+1,\qquad
     \Delta D_a>0\Rightarrow \Delta_a\ge 0.$$

4. **Reduction.** Modulo the Structural Lemma, the Core Lemma is equivalent to a purely
   local **Forest Inequality** on the chains. The single missing ingredient is a
   *parent–child realizability constraint* (the pure inequality is false without it; I
   exhibit the counterexample).

The upshot: the Core Lemma is reduced to a finite, local, elementary statement about
chains of integer pairs $(\Delta_a,\Delta D_a)$ with $\Delta_a,\Delta D_a\in\{-4,\dots,4\}$,
and all but one piece is proved.

---

## 1. Setup (recap)

$\lambda\vdash n$, rows/columns $0$-indexed, content $\mathrm{cont}(r,c)=c-r$,
$\mathrm{cont}_T(v)$ the content of the cell of $v$ in $T\in\mathrm{SYT}(\lambda)$.
Descent indicator $\delta_T(i)=[\mathrm{cont}_T(i{+}1)<\mathrm{cont}_T(i)]$.
Staircase word $\mathbf w_0=(1)(2\,1)\cdots(n{-}1\cdots1)$, length $\binom n2$; **block $b$**
is the segment $(b,b{-}1,\dots,1)$, so comparator $i$ occurs once in each block $b\ge i$.
A closed walk $U_0\to\cdots\to U_L=U_0$ takes one STAY or SWAP step per letter; SWAP at $i$
replaces $U$ by $s_iU$ (swap the cells of labels $i,i{+}1$) and is legal iff
$|d_U(i)|\ge2$, $d_U(i)=\mathrm{cont}_U(i{+}1)-\mathrm{cont}_U(i)$. The walk is **feasible**
iff no before-state has $d=-1$. After-state cost $\mathrm{cost}=\sum_k\delta_{U_k}(i_k)$.
The all-stay cost of a fixed tableau is $D(T)=\sum_i(n-i)\delta_T(i)$.

**All-stay lemma (proved earlier, cited):** $D(T)\ge n(\lambda)$ for *every* $T$, equality
at $T_{\mathrm{rs}}$.

---

## 2. Structural Lemma (verified $n\le 7$)

**Lemma 2.1.** *In every feasible closed walk:*
*(a)* *the multiset of swap steps, read in walk order, is a sequence of comparators that
reduces to the empty word by repeatedly cancelling two stack-adjacent equal comparators
(equivalently: the swaps match into nested **arches**, each arch a pair of swaps at one
comparator $i$, opening at a block $b_1$ and closing at a block $b_2>b_1$);*
*(b)* *the nesting forest is a union of chains — each arch contains at most one child arch;*
*(c)* *the forest has at most $\lfloor (n-1)/?\rfloor$ roots (empirically $\le1$ for
$n\le6$, $\le2$ for $n=7$).*

*Status:* exhaustively verified for all feasible closed walks $n\le 6$ and a broad $n=7$
sample (576 walks); zero violations. **Proof open** (the rigidity should follow from
feasibility + closure; see the 0-Hecke/Coxeter-sorting remark in the companion paper).

Throughout the rest, fix a feasible closed walk and an arch $a$: comparator $i$, opening at
block $b_1$, closing at block $b_2$, **boundary state** $U=U_a$ (the state just before the
opening swap) and **inner state** $V=V_a=s_iU$. Write $w=w_a:=b_2-b_1\ge1$ and $ni:=n-i$.
Since $b_1\ge i$ and $b_2\le n-1$,
$$\boxed{\,1\le w\le ni-1\,}.\tag{2.1}$$

For an **innermost** arch the region strictly between its two swaps is all-stay at $V$.

---

## 3. The per-arch quantities and their exact algebra (PROVED)

Define, at the boundary of arch $a$ (contents $x_j:=\mathrm{cont}_U(j)$):
$$\gamma:=\delta_V(i)-\delta_U(i),\quad
  \alpha:=\delta_V(i{-}1)-\delta_U(i{-}1),\quad
  \beta:=\delta_V(i{+}1)-\delta_U(i{+}1)$$
(with $\alpha:=0$ if $i=1$, $\beta:=0$ if $i+1=n$), and $\sigma:=\gamma+\alpha+\beta$.

**Lemma 3.1 (sign rule + $\sigma$-linkage; PROVED).**
$$\gamma\in\{-1,+1\},\quad \gamma=+1\iff x_i<x_{i+1};\qquad \alpha\gamma\le0,\ \beta\gamma\le0;$$
$$\sigma\in\{-1,0,1\};\qquad \sigma=\pm1\Rightarrow \alpha-\beta=0,\qquad
  \sigma=0\Rightarrow|\alpha-\beta|=1.$$

*Proof.* The swap is legal so $x_i\neq x_{i+1}$; toggling the order of labels $i,i{+}1$
flips the descent at $i$, so $\delta_V(i)=1-\delta_U(i)$ and $\gamma=1-2\delta_U(i)=\pm1$,
with $\gamma=+1\iff\delta_U(i)=0\iff x_i<x_{i+1}$. In $V$, label $i$ carries content
$x_{i+1}$ and label $i{+}1$ carries content $x_i$, while $x_{i-1},x_{i+2}$ are unchanged.
Hence
$\alpha=[x_{i+1}<x_{i-1}]-[x_i<x_{i-1}]$ and
$\beta=[x_{i+2}<x_i]-[x_{i+2}<x_{i+1}]$.
If $\gamma=+1$ then $x_{i+1}>x_i$, so both bracket-differences are $\le0$; if $\gamma=-1$
then $x_{i+1}<x_i$, so both are $\ge0$. Thus $\alpha,\beta\in\{0,-\gamma\}$, giving
$\alpha\gamma\le0,\ \beta\gamma\le0$ and $\sigma=\gamma+\alpha+\beta\in\{-1,0,1\}$.
Enumerating the four sign-consistent $(\alpha,\beta)$ in each case of $\gamma$ gives the
linkage: $\sigma=\pm1$ forces $\{\alpha,\beta\}$ symmetric ($\alpha=\beta$, so
$\alpha-\beta=0$); $\sigma=0$ forces exactly one of $\alpha,\beta$ nonzero
($|\alpha-\beta|=1$). $\qquad\blacksquare$

*(Brute-force confirmation: all 3800 legal swaps over every SYT, $n=3..8$, satisfy
Lemma 3.1 with zero exceptions.)*

**Lemma 3.2 ($D$-jump; PROVED).** $\displaystyle \Delta D_a:=D(V)-D(U)=ni\cdot\sigma+(\alpha-\beta).$

*Proof.* $U,V$ differ only at comparators $i{-}1,i,i{+}1$, so
$D(V)-D(U)=(ni{+}1)\alpha+ni\,\gamma+(ni{-}1)\beta=ni(\gamma+\alpha+\beta)+(\alpha-\beta)$.
$\qquad\blacksquare$

**Lemma 3.3 (local cost-change / telescoping step; PROVED).**
*Removing an innermost arch (replacing the inner all-stay segment at $V$ by all-stay at $U$,
and dropping the two swaps) changes the cost by exactly*
$$\Delta_a=w\cdot\sigma=(b_2-b_1)(\gamma+\alpha+\beta).$$

*Proof.* The region strictly between the two comparator-$i$ swaps contains, by the block
structure of $\mathbf w_0$, exactly $b_2-b_1-1$ occurrences of comparator $i$ and $b_2-b_1$
occurrences each of $i{-}1$ and $i{+}1$ (count: the tail $i{-}1,\dots,1$ of block $b_1$,
the full blocks $b_1{+}1,\dots,b_2{-}1$, and the head $b_2,\dots,i{+}1$ of block $b_2$).
On each such step the cost difference (state $V$ minus state $U$) is
$\delta_V(\cdot)-\delta_U(\cdot)$, i.e. $\gamma$ at comparator $i$, $\alpha$ at $i{-}1$,
$\beta$ at $i{+}1$; the two swap steps themselves contribute $\delta_V(i)$ (open) and
$\delta_U(i)$ (close) vs. two all-stay $\delta_U(i)$, a net $\gamma$. Summing,
$$\Delta_a=\gamma+\big[(b_2-b_1-1)\gamma+(b_2-b_1)\alpha+(b_2-b_1)\beta\big]
        =(b_2-b_1)(\gamma+\alpha+\beta).\qquad\blacksquare$$

Removing innermost arches one at a time leaves all other arches' boundary/inner states (hence
their $\Delta$) unchanged, so by induction:

**Corollary 3.4 (Telescoping identity; PROVED, given Lemma 2.1(a)).**
$\mathrm{cost}(W)=D(U_0)+\sum_{\text{arches}}\Delta_a.$

**Theorem 3.5 (Refined per-arch inequality; PROVED).**
$$\Delta D_a\neq0;\qquad \Delta D_a<0\Rightarrow\Delta_a\ge\Delta D_a+1;\qquad
  \Delta D_a>0\Rightarrow\Delta_a\ge0.$$

*Proof.* By Lemmas 3.1–3.3, $\Delta_a=w\sigma$ and $\Delta D_a=ni\,\sigma+(\alpha-\beta)$
with $\sigma\in\{-1,0,1\}$, $1\le w\le ni-1$ (eqn 2.1).
- $\sigma=0$: $\Delta D_a=\alpha-\beta=\pm1\neq0$, $\Delta_a=0$. If $\Delta D_a=+1>0$,
  $\Delta_a=0\ge0$. If $\Delta D_a=-1<0$, $\Delta_a=0\ge\Delta D_a+1=0$.
- $\sigma=+1$: $\alpha-\beta=0$, so $\Delta D_a=ni\ge1>0$ and $\Delta_a=w\ge1\ge0$.
- $\sigma=-1$: $\alpha-\beta=0$, so $\Delta D_a=-ni<0$ and $\Delta_a=-w$. By (2.1),
  $w\le ni-1$, hence $\Delta_a=-w\ge-(ni-1)=\Delta D_a+1$. $\qquad\blacksquare$

---

## 4. Reduction to the Forest Inequality

By Lemma 2.1 the arches form a forest of chains. Define recursively, per arch $a$ with
(at most one) child $c(a)$:
$$T(a)=\Delta_a+T(c(a)),\qquad M(a)=\Delta D_a+\min\big(0,\,M(c(a))\big)$$
($T=\Delta$, $M=\Delta D$ at a leaf). Then $M(a)$ is the deepest cumulative $D$-dip inside
$a$'s subtree (relative to $a$'s parent level), and the global minimum visited $D$ satisfies
$$\min_t D(U_t)=D(U_0)+\mu,\qquad \mu:=\min\Big(0,\ \min_{\text{roots }r}M(r)\Big).$$
By Corollary 3.4, $\mathrm{cost}-D(U_0)=\sum_{\text{roots}}T(r)$. Therefore:

**Reduction.** *Modulo Lemma 2.1, the Core Lemma is equivalent to the*
$$\textbf{Forest Inequality:}\qquad \sum_{\text{roots }r}T(r)\ \ge\ \mu.$$

**Verified:** zero violations over all feasible closed walks $n\le7$. Single-root slack
distribution $\{0\!:\!39,\,1\!:\!27,\,2\!:\!11\}$ ($n\le6$); the bound is frequently tight.

### 4.1 The per-arch inequality alone is NOT sufficient — parent/child realizability

The Forest Inequality does **not** follow from Theorem 3.5 as a free-standing inequality on
arbitrary chains. Counterexample (a length-2 chain): parent $\sigma=+1$, $ni=3$, $w=1$
($\Delta=+1,\ \Delta D=+3$); child $\sigma=-1$, $ni=4$, $w=3$ ($\Delta=-3,\ \Delta D=-4$).
Both satisfy Theorem 3.5 and (2.1). Cumulative dips $C_1=3,\ C_2=-1$, so $\mu=-1$, but
$\sum\Delta=1-3=-2<\mu$. **No feasible walk realizes this chain** (it never appears in the
$n\le7$ enumeration). The missing fact is a *parent–child realizability constraint*: when a
child arch is nested inside its parent, the boundary state of the child equals the parent's
inner state $V_{\text{parent}}$, which couples their $(\,i,\sigma,w\,)$ data. Quantifying
this coupling is exactly the remaining gap.

A clean sufficient package (each piece verified $n\le7$, all unproven):
- **(I2)** $T(a)\ge\min(0,M(a))$ for every subtree;
- **(C2)** $M(a)\ge0\Rightarrow T(a)\ge0$;
- a **multi-root combining bound**: $\sum_r T(r)\ge\mu$ even when several roots dip below
  $0$ (e.g. observed $[(-1,-4),(-1,-3)]$: $\sum T=-2\ge\mu=-4$; the deepest root carries
  enough slack $T(r)-M(r)$ to absorb the others). At most $2$ roots occur through $n=7$.

(I2) at a single root *is* the single-root Forest Inequality; it would follow from the
parent–child coupling by the chain recursion.

---

## 5. Status and the precise remaining gap

**Proved unconditionally this session:** Lemma 3.1 (sign rule + $\sigma$-linkage),
Lemma 3.2 ($D$-jump formula), Lemma 3.3 (local cost-change $\Delta_a=w\sigma$),
Corollary 3.4 (telescoping, given Lemma 2.1(a)), Theorem 3.5 (refined per-arch inequality).

**Verified $n\le7$, used as hypotheses:** Lemma 2.1 (nesting + chain forest + $\le2$ roots),
the Forest Inequality.

**The remaining gap, sharply:** prove the **parent–child realizability constraint** that
upgrades Theorem 3.5 to the chain version of (I2), then dispatch the (at most two) roots by
the slack-absorption bound. Equivalently: characterise which $(\sigma,w,i)$ sequences arise
as nested chains of a feasible closed walk, and show they exclude the deep-child/shallow-
parent configuration of §4.1. This is a *finite local* question — a strict tightening of the
"feasible closed walks are rare and rigid" phenomenon — and is the natural target for the
next session (route: read the chain as a Coxeter-sorting subword returning to the identity,
and track the axial-distance $d\ne-1$ constraint down the nesting).

This supersedes the min-$D$ Key Lemma reduction: that reduced the Core Lemma to a per-walk
inequality verified but opaque; the present reduction makes the obstruction **local and
explicit** (one inequality per nested parent/child pair) and proves everything except that
single local coupling.

---

## 6. Addendum (2026-05-22 wake): the width-half is PROVED; residual is content-coupling

The "parent–child realizability constraint" of §4.1 splits cleanly, via one structural
observation, into a width part (now **proved**) and a content part (the sharpened residual).

**Boundary identity.** In a chain, between a parent's opening swap and its child's opening
swap every step is a STAY, and STAYs fix the tableau. Hence the child's boundary state equals
the parent's inner state: $U_c=V_p=s_{i_p}U_p$. This couples *all* child data
$(\sigma_c,ni_c,w_c,e_c)$ — read off contents of $U_c=V_p$ — to the parent.

**Lemma 6.1 (width coupling; PROVED).** $w_{c(a)}\le w_a-1$.
*Proof.* Index letters of $\mathbf w_0$ by $(b,i)$ with step order $(b,i)\prec(b',i')$ iff
$b<b'$, or $b=b'\wedge i>i'$ (within a block, comparators descend). Nesting gives
$(b_1^p,i_p)\prec(b_1^c,i_c)$ and $(b_2^c,i_c)\prec(b_2^p,i_p)$, so $b_1^c\ge b_1^p$ (equality
only if $i_c<i_p$) and $b_2^c\le b_2^p$ (equality only if $i_c>i_p$). Both equalities at once
need $i_c<i_p$ and $i_c>i_p$, impossible; if $i_c=i_p$ neither holds. So at least one is
strict and $w_c=b_2^c-b_1^c\le w_p-1$. $\blacksquare$ (Verified $n\le7$, tight: equality in
82/123 pairs at $n=7$.)

**Consequence.** Lemma 6.1 alone kills the §4.1 fake chain ($w_c=3\not\le w_p-1=0$) and
*every* length-2 violator. After imposing $w_c\le w_p-1$, an abstract enumeration finds the
only surviving violators are length-$\ge3$: two $(\sigma{=}0,e{=}{+}1)$ arches over a
$(\sigma{=}-1,ni{=}2,w{=}1)$ leaf, where $e:=\alpha-\beta$.

**Residual (verified $n\le7$, open).** Two content-coupling rules — each a consequence-to-prove
of the boundary identity $U_c=V_p$ — eliminate these:
- **(R)** a $\sigma{=}0$ child of a $\sigma{=}0$ parent inherits $e$: $e_c=e_p$;
- **(S2)** a $(\sigma{=}0,e{=}{+}1)$ arch under a $\sigma{=}0$ parent is a leaf.

Realized $(\sigma_p,\sigma_c)$ adjacency over all feasible closed walks $n\le7$ (653 walks,
584 chains, 64 signatures): par$-1\to\{0\!:\!4,+1\!:\!6\}$; par$0\to\{-1\!:\!25,0\!:\!25,+1\!:\!6\}$;
par$+1\to\{-1\!:\!51,0\!:\!10\}$; the pairs $(+1,+1)$ and $(-1,-1)$ never occur.

**Net status.** The width-geometry of the gap is fully proved; the remaining gap is purely a
**content-coupling** statement — derive (R),(S2) from the content arithmetic of
$V_p=s_{i_p}U_p$ (how the descent flips at $i_p\pm1$ propagate to the child comparator $i_c$).
This is the next prove target.
