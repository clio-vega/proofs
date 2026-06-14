# The general-$c$ Number Lemma (NL$_c$): a two-line falling-factorial 2-adic bound

**Date:** 2026-06-14 · **Status:** PROVED (all $c\ge 1$, all even $F\ge2$, all $j$) · **Corollary:**
c=3 Gap 1 (L3″-A) fully closed.

For a positive integer $n$, $v_2(n)$ is the 2-adic valuation; $v_2(0)=+\infty$. The falling factorial
is $j^{(k)} = j(j-1)\cdots(j-k+1) = k!\binom{j}{k}$, and $s_2(n)$ is the binary digit sum, so by
Legendre $v_2(n!) = n - s_2(n)$. We write $C(n,k)=\binom{n}{k}$, with $C(n,k)=0$ when $k<0$ or $k>n$.

---

## 1. Statement

> **Theorem (NL$_c$).** Let $c\ge 1$, let $F\ge 2$ be **even**, and let $2c \le j \le F+2c-1$. Then
> $$ v_2\binom{F+2c-1}{j} \;+\; v_2\!\big(j^{(2c)}\big) \;\ge\; v_2(F) + \beta(c), $$
> where
> $$ \boxed{\ \beta(c) \;=\; (c-1) + v_2\big((c-1)!\big) \;=\; 2(c-1) - s_2(c-1).\ } $$
> The constant $\beta(c)$ is **sharp**: equality is attained.

First values: $\beta(1)=0,\ \beta(2)=1,\ \beta(3)=3,\ \beta(4)=4,\ \beta(5)=7,\ \beta(6)=8,\dots$
(matching $v_2((c-1)!)=0,0,1,1,3,3,4,\dots$). The $c=1$ value $\beta(1)=0$ recovers the hook/two-row
bound $R\ge 0$; the $c=2$ value $\beta(2)=1$ recovers the proved-and-Lean-checked Number Lemma
`v₂C(F+3,j)+v₂(j(j−1)(j−2)(j−3)) ≥ v₂(F)+1`.

The quantity $j^{(2c)} = (2c)!\,C(j,2c)$ is exactly the **inhomogeneous tip** $-(2c)!\,C(j,2c)$ of the
sextic/$2c$-tic $Q_c$ in the three-row even-$|J^*|$ analysis. NL$_c$ is the single estimate that, in
every $c$ at once, controls that tip.

---

## 2. Proof

The proof is the $c=2$ proof lifted verbatim, with the only genuinely new ingredient — the **anchor
bound (move 3)** — replaced by an *exact identity* that makes $\beta(c)$ fall out and the sharpness
transparent.

### Move 1 — the subset-of-a-subset identity (load-bearing)

For $2c \le j \le F+2c-1$ (so $0 \le j-2c \le F-1$),
$$ \binom{F+2c-1}{j}\binom{j}{2c} \;=\; \binom{F+2c-1}{2c}\binom{F-1}{\,j-2c\,}. \tag{$\star$}$$
This is "choose $j$ of $F+2c-1$, then $2c$ of those $j$" $=$ "choose $2c$ of $F+2c-1$, then the
remaining $j-2c$ from the other $F-1$": both equal $\frac{(F+2c-1)!}{(2c)!\,(j-2c)!\,(F-1-(j-2c))!}$.
Taking $v_2$ and using $v_2\binom{F-1}{j-2c}\ge 0$,
$$ v_2\binom{F+2c-1}{j} \;=\; v_2\binom{F+2c-1}{2c} + v_2\binom{F-1}{j-2c} - v_2\binom{j}{2c}
\;\ge\; v_2\binom{F+2c-1}{2c} - v_2\binom{j}{2c}. \tag{1}$$

### Move 2 — falling factorial to central binomial

$j^{(2c)} = (2c)!\,\binom{j}{2c}$, so
$$ v_2\!\big(j^{(2c)}\big) \;=\; v_2\big((2c)!\big) + v_2\binom{j}{2c}. \tag{2}$$

Adding (1) and (2), the term $v_2\binom{j}{2c}$ — the only $j$-dependence in the lower bound — cancels:
$$ v_2\binom{F+2c-1}{j} + v_2\!\big(j^{(2c)}\big) \;\ge\; v_2\binom{F+2c-1}{2c} + v_2\big((2c)!\big). \tag{3}$$
The right side is **independent of $j$**. Everything now rests on it.

### Move 3 — the anchor identity, and $\beta(c)$

The numerator of $\binom{F+2c-1}{2c}$ is a product of $2c$ **consecutive integers starting at $F$**:
$$ \binom{F+2c-1}{2c} \;=\; \frac{(F+2c-1)!}{(2c)!\,(F-1)!} \;=\; \frac{F(F+1)\cdots(F+2c-1)}{(2c)!}. $$
Hence $\binom{F+2c-1}{2c}\cdot(2c)! = \prod_{i=0}^{2c-1}(F+i)$, and the right side of (3) is an **exact
sum of valuations**:
$$ v_2\binom{F+2c-1}{2c} + v_2\big((2c)!\big) \;=\; \sum_{i=0}^{2c-1} v_2(F+i)
\;=\; v_2(F) + \sum_{i=1}^{2c-1} v_2(F+i). \tag{4}$$

It remains to bound $\displaystyle S(F):=\sum_{i=1}^{2c-1} v_2(F+i)$ below by $\beta(c)$, for every even
$F$. Since $F$ is even, the odd offsets $i\in\{1,3,\dots,2c-1\}$ contribute $v_2(F+\text{odd})=0$; only
the $c-1$ even offsets $i=2,4,\dots,2c-2$ survive. Writing $F=2G$ ($G\ge 1$):
$$ S(F) \;=\; \sum_{k=1}^{c-1} v_2(F+2k) \;=\; \sum_{k=1}^{c-1}\big(1 + v_2(G+k)\big)
\;=\; (c-1) + v_2\!\Big(\prod_{k=1}^{c-1}(G+k)\Big). $$
The product $\prod_{k=1}^{c-1}(G+k) = \frac{(G+c-1)!}{G!} = \binom{G+c-1}{c-1}\,(c-1)!$, so
$$ S(F) \;=\; (c-1) + v_2\binom{G+c-1}{c-1} + v_2\big((c-1)!\big) \;=\; \beta(c) + v_2\binom{G+c-1}{c-1}
\;\ge\; \beta(c), \tag{5}$$
because $v_2\binom{G+c-1}{c-1}\ge 0$. This holds for **every** even $F\ge 2$. Combining (3), (4), (5):
$$ v_2\binom{F+2c-1}{j} + v_2\!\big(j^{(2c)}\big) \;\ge\; v_2(F) + \beta(c). \qquad\blacksquare$$

### Sharpness

Tracking the two discarded non-negative terms — $v_2\binom{F-1}{j-2c}$ in (1) and
$v_2\binom{G+c-1}{c-1}$ in (5) — gives the **exact slack**
$$ \Big[v_2\tbinom{F+2c-1}{j} + v_2(j^{(2c)})\Big] - \big[v_2(F)+\beta(c)\big]
\;=\; v_2\binom{F-1}{j-2c} + v_2\binom{F/2+c-1}{c-1}. $$
Both Kummer terms vanish simultaneously: take $G=F/2=2^t$ with $2^t \ge c$ (then $\binom{G+c-1}{c-1}$
is odd, since adding $c-1<2^t$ to $2^t$ carries nowhere), and $j=2c$ (then $\binom{F-1}{0}=1$).
So equality is attained and $\beta(c)$ cannot be raised. $\square$

**Remark (the broken assumption).** The PROVE brief warned the likely error was *"$v_2\binom{F+2c-1}{2c}$
behaves like the $c=2$ case uniformly in $c$."* It does not behave like a fixed shift — but the *clean*
object is not $\binom{F+2c-1}{2c}$ in isolation. It is the **combination** $v_2\binom{F+2c-1}{2c} +
v_2((2c)!) = \sum_{i=0}^{2c-1}v_2(F+i)$, which is exact and $c$-uniform. Pinning $\beta(c)$ is then a
one-line minimisation over $c-1$ consecutive integers $G+1,\dots,G+c-1$, whose excess valuation over
$v_2((c-1)!)$ is exactly $v_2\binom{G+c-1}{c-1}\ge 0$. The $c=2$ ad-hoc count ("only $F,F+2$ even") is
the $c=2$ shadow of this identity.

---

## 3. Corollary: c=3 Gap 1 (L3″-A) fully closed

Recall the proved structural identity for the three-row $c=3$ sextic (`2026-06-14-threerow-c3`):
$$ Q_3(a,b,j) \;=\; (a-1)(b-2)\,H(a,b,j)\; -\; 720\,C(j,6), \qquad 720\,C(j,6)=j^{(6)}, $$
$$ H(a,b,j) = (a{+}3)(a{+}4)(b{+}2)(b{+}3) - 6(a{+}3)(b{+}2)\,C(j,1) - 6(ab{+}a{+}2b)\,C(j,2)
+ 36\,C(j,3) + 72\,C(j,4). $$

> **Corollary (L3″-A / Gap 1).** For $a$ **even**, $b\ge 2$, and $0 \le j \le b+3$,
> $$ v_2\binom{b+3}{j} + v_2 Q_3(a,b,j) \;\ge\; v_2(b+3) + 1. $$

**Proof.** Multiply through by $\binom{b+3}{j}$ and split along the two-generator identity:
$$ \binom{b+3}{j}\,Q_3(a,b,j) \;=\; \underbrace{(a-1)(b-2)\,\binom{b+3}{j}H(a,b,j)}_{T_1}
\;-\; \underbrace{\binom{b+3}{j}\,j^{(6)}}_{T_2}. $$
Since $v_2$ of a difference is $\ge$ the minimum of the two valuations, it suffices to show
$v_2 T_1 \ge v_2(b+3)+1$ **and** $v_2 T_2 \ge v_2(b+3)+1$.

**The tip $T_2$ (where NL$_3$ bites).** By the subset identity ($\star$) with $c=3$,
$\binom{b+3}{j} C(j,6) = \binom{b+3}{6}\binom{b-3}{j-6}$, and $720\binom{b+3}{6} = (b+3)^{(6)}$, so
$$ T_2 \;=\; (b+3)^{(6)}\binom{b-3}{j-6} \;=\; (b+3)\cdot\big[(b+2)(b+1)b(b-1)(b-2)\big]\binom{b-3}{j-6}. $$
The bracket is a product of **five consecutive integers**, so $v_2\ge v_2(5!)=3$; hence
$v_2 T_2 \ge v_2(b+3)+3 \ge v_2(b+3)+1$. *(This is precisely the NL$_3$ mechanism — anchor $\beta(3)=3$ —
realised directly, with no parity hypothesis on $F=b-2$.)*

**The heavy part $T_1$ (where "the heavy factor only helps").** Convert $H$ to the falling-factorial
basis,
$$ H = A - 6(a{+}3)(b{+}2)\,j^{(1)} - 3(ab{+}a{+}2b)\,j^{(2)} + 6\,j^{(3)} + 3\,j^{(4)},
\qquad A=(a{+}3)(a{+}4)(b{+}2)(b{+}3), $$
and expand $\binom{b+3}{j}H$ term-by-term with the rearranged subset identity
$\binom{b+3}{j}\,j^{(k)} = (b+3)^{(k)}\binom{b+3-k}{\,j-k\,}$:

| term of $\binom{b+3}{j}H$ | value | factor of $b+3$ | extra factor of $2$ | $v_2 \ge$ |
|---|---|---|---|---|
| $A\binom{b+3}{j}$ | $(a{+}3)(a{+}4)(b{+}2)(b{+}3)\binom{b+3}{j}$ | $b{+}3$ | $a{+}4$ even ($a$ even) | $v_2(b{+}3)+1$ |
| $-6(a{+}3)(b{+}2)\,j^{(1)}\binom{b+3}{j}$ | $-6(a{+}3)(b{+}2)(b{+}3)\binom{b+2}{j-1}$ | $b{+}3$ | coeff $6$ | $v_2(b{+}3)+1$ |
| $-3(ab{+}a{+}2b)\,j^{(2)}\binom{b+3}{j}$ | $-3(ab{+}a{+}2b)(b{+}3)(b{+}2)\binom{b+1}{j-2}$ | $b{+}3$ | $ab{+}a{+}2b$ even ($a$ even) | $v_2(b{+}3)+1$ |
| $6\,j^{(3)}\binom{b+3}{j}$ | $6(b{+}3)(b{+}2)(b{+}1)\binom{b}{j-3}$ | $b{+}3$ | coeff $6$ | $v_2(b{+}3)+1$ |
| $3\,j^{(4)}\binom{b+3}{j}$ | $3(b{+}3)(b{+}2)(b{+}1)b\binom{b-1}{j-4}$ | $b{+}3$ | $\{b,b{+}1,b{+}2\}$ has an even | $v_2(b{+}3)+1$ |

Every term has $v_2 \ge v_2(b+3)+1$, so $v_2\big(\binom{b+3}{j}H\big)\ge v_2(b+3)+1$, whence
$v_2 T_1 = v_2\big((a-1)(b-2)\big) + v_2\big(\binom{b+3}{j}H\big) \ge v_2(b+3)+1$ (the heavy factor
$(a-1)(b-2)$ contributes $\ge 0$; for $a$ even it is in fact odd-times-$(b-2)$, and is genuinely
unneeded — the "$a$ even" hypothesis enters only through $A$ and the cubic-coefficient $ab+a+2b$).

Therefore $v_2\big(\binom{b+3}{j}Q_3\big) = v_2(T_1-T_2) \ge \min(v_2 T_1, v_2 T_2) \ge v_2(b+3)+1$.
$\blacksquare$

This closes Gap 1 of `2026-06-14-threerow-c3-Jstar-even.md` (Compensation Lemma A, $a$ even branch),
unconditionally and for all $b$.

---

## 4. Verification

All checks: 0 failures.

- **NL$_c$** and **sharpness** (min slack $=0$): all even $F\le 198$, all $2c\le j\le F+2c-1$, for
  $c=1,2,3,4,5$.
- **$\beta(c)$ closed form** $(c-1)+v_2((c-1)!) = 2(c-1)-s_2(c-1)$: $c\le 8$.
- **Subset identity ($\star$)** and the **anchor identity (4)** $v_2\binom{F+2c-1}{2c}+v_2((2c)!) =
  \sum_{i=0}^{2c-1}v_2(F+i)$: $c\le 5$, even $F\le 78$.
- **Anchor minimisation (5)** $\min_{F\text{ even}} S(F) = (c-1)+v_2((c-1)!)$: $c\le 6$, $F\le 4998$.
- **Corollary (Gap 1):** the term-by-term identity $\sum(\text{5 terms})=\binom{b+3}{j}H$ and each
  term's bound $v_2\ge v_2(b+3)+1$, plus the full statement $v_2\binom{b+3}{j}+v_2 Q_3 \ge
  v_2(b+3)+1$: all $a$ even $\le 78$, $2\le b\le a$, $0\le j\le b+3$ ($\approx 2.4\times10^5$ triples).

Scripts in `2026-06-14-numberlemma-general-c-code/`.

---

## 5. What this closes, and what remains

**Closed (the Win).** NL$_c$ holds for all $c\ge 1$ with the explicit, sharp constant
$\beta(c)=(c-1)+v_2((c-1)!)$. This is the single falling-factorial tip bound that governs the
**interior** of the three-row even-$|J^*|$ family for *every* $c$ simultaneously: wherever the tip
$-(2c)!\,C(j,2c)$ of $Q_c$ must be 2-adically dominated, NL$_c$ supplies the bound. Gap 1 of the $c=3$
write-up (the $a$-even Compensation Lemma A) is now a corollary, proved unconditionally.

**Remaining (not an NL instance).** Gap 2 of the $c=3$ write-up — **Compensation Lemma B**, the
$a$-odd branch $\tilde\Delta(j)=j+3-2s_2(j)+2U(j)\ge 0$ for $4\le j\le b$ — is the genuine
**two-generator** inequality (simultaneous deficits at $v_2(j)=1$ and $v_2(j)=2$). It is *not* a pure
tip bound: it couples $s_2(j)$ with the full Prop-2 valuation $U$, and the per-$v_2(j)$ minima of $U$
are the irregular sequence $\{0:-6,\,1:-5,\,2:-5,\,3:-4,\,4:-5,\dots\}$ (not a clean $c-v_2(j)$). NL$_3$
controls the tip inside $U$ but does not by itself settle the $s_2(j)$-coupled superposition; that —
and the tie classification S2 (box $=\langle 2,4\rangle$, $|J^*|=4$ iff $a\equiv1,\,b\equiv2 \pmod 4$)
— remains open (verified $m\le 79$, ties $m\le 39$). This is the **$e_2\bmod 2$ wall** proper.

**Methodological note.** As anticipated, *one subset identity plus one Kummer/Legendre count beat the
whole-series machinery* — at $c=2$ and now uniformly in $c$. The lift required no new idea beyond
recognising that the anchor is an exact sum of valuations over consecutive integers, and that the
two-generator $Q_3$ tip and heavy parts *both* dissolve under the same identity ($\star$) read in two
directions (block size $2c$ for the tip, block sizes $1,2,3,4$ for $H$).
