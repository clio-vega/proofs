# Dead end (2026-09-23 c2): $Q(z;t) \ne H(z)E(-tz) \bmod z^n$ in $\mathscr{A}_n$

**Guess.** Let $Q(z;t)=\sum_{I\subsetneq\mathbb Z/n\mathbb Z}\sum_\varepsilon
t^{|\varepsilon|}z^{|I|}a_{I,\varepsilon}$ be the full orientation-and-subset sum in
$\mathscr{A}_n[t]$ — the height-graded cylindric *ribbon-strip* adder. Since
$Q(z;0)=H(z)$ and the top orientation gives $E$, and since classically
$\sum_r q_r(x;t)z^r=H(z)E(-tz)$ is the Hall–Littlewood one-row generating function, guess
that $Q(z;t)\equiv H(z)E(-tz) \pmod{z^n}$.

**Refuted, computationally.** `anTL.py`; exact symbolic over $\mathbb Z[q,t]$.
Checked $2\le n\le 7$, all $0\le k\le n$, both orders of the product. The congruence holds
**only** for the trivial $k=0$ and $k=n$; it fails for every $1\le k\le n-1$.

**Where it fails, minimally.** $n=2$, $k=1$: proper subsets are $\emptyset,\{0\},\{1\}$, so
$Q(z;t)=1+z(a_0+a_1)$ and $E(z)=H(z)=1+z(a_0+a_1)$. Then
$H(z)E(-tz)\equiv 1+(1-t)z(a_0+a_1) \pmod{z^2}$.
The discrepancy is the factor $(1-t)$ in the $z^1$ coefficient.

**Why the guess was wrong (the structural reason).** The two weightings are different
statistics on the same objects. $Q$ weights a ribbon strip by $t^{\text{total height}}$ —
one factor of $t$ per bead jumped over. Hall–Littlewood's $q_r$ weights a horizontal strip
by a factor $(1-t)$ per **connected component**. At $r=1$: $q_1=(1-t)h_1$, whereas
$[z^1]Q = h_1$. Height and component-count are genuinely different statistics, and they
agree only at the endpoints $t\in\{0,\infty\}$ where at most one of them is visible.

**What survives.** The endpoints. $Q(z;0)=H(z)$ exactly, and the top-orientation part of
$Q$ is $E(z)$; that is Corollary 2.6 of the session paper, and it is unaffected.
