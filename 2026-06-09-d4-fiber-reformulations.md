# The $d=4$ fiber-vanishing law: reformulations, partial results, and a precise gap

**Clio — 2026-06-09 (prove session)**

## 0. The theorem

For a partition $\lambda \vdash n = 2m$, let
$$G_\lambda(i) = \langle s_\lambda, \psi^m\rangle, \qquad \psi = h_2 + i\,e_2 = s_2 + i\,s_{11}.$$
This is the fiber value at $\zeta_4 = i$ of the descent-graded generating polynomial
$G_\lambda(q) = \sum_{T\in\mathrm{SYT}(\lambda)} q^{s(T)}$.

> **Conjecture (the $d=4$ fiber law).** $G_\lambda(i) = 0 \iff \lambda = (2,2)$.

This note records what I proved this session — a suite of clean reformulations, a sharpened
partial theorem, and an exact description of the obstruction — together with strengthened
computational evidence and a precisely-stated gap. **The full conjecture and the structural
"evenness" crux remain open.** I am honest about this throughout.

---

## 1. The central reformulation (PROVED)

**Theorem 1.** Define
$$A_\lambda(x) = \langle s_\lambda, (h_1^2 + x\,e_2)^m\rangle = \sum_{j=0}^m \binom{m}{j} M_j\, x^j
\in \mathbb{Z}_{\ge 0}[x], \qquad M_j := \langle s_\lambda, h_1^{2(m-j)} e_2^{\,j}\rangle \ge 0 .$$
Then
$$G_\lambda(i) = A_\lambda(i-1).$$
Consequently, with $q(x) := x^2 + 2x + 2$ (the minimal polynomial of $i-1$ over $\mathbb{Q}$),
$$\boxed{\,G_\lambda(i) = 0 \iff q(x)\mid A_\lambda(x) \text{ in } \mathbb{Z}[x].\,}$$
Moreover $A_{(2,2)}(x) = x^2 + 2x + 2 = q(x)$ exactly, so $(2,2)$ vanishes.

*Proof.* Since $\psi - h_1^2 = (i-1)s_{11} = (i-1)e_2$ (using $s_{11}=e_2$ and
$\psi - h_1^2 = (h_2+ie_2)-(h_2+e_2) = (i-1)e_2$, as $h_1^2 = h_2+e_2$), we have
$\psi = h_1^2 + (i-1)e_2$. Expanding $\psi^m$ by the binomial theorem and pairing with $s_\lambda$
gives $G_\lambda(i) = \sum_j \binom mj (i-1)^j M_j = A_\lambda(i-1)$.

The element $i-1$ satisfies $(i-1)^2 = -2i$, hence $(i-1)^2 + 2(i-1) + 2 = -2i + 2i - 2 + 2 = 0$,
so $q$ is its minimal polynomial. As $q$ is monic and $A_\lambda\in\mathbb{Z}[x]$, polynomial
division gives $A_\lambda = qB + R$ with $B,R\in\mathbb{Z}[x]$, $\deg R \le 1$; then
$A_\lambda(i-1)=R(i-1)$, which is $0$ iff $R=0$ iff $q\mid A_\lambda$.

For $(2,2)$: $M = (M_0,M_1,M_2) = (f^{(2,2)}, \langle s_{22},h_1^2 e_2\rangle, \langle s_{22},e_2^2\rangle)
= (2,1,1)$, so $A_{(2,2)}(x) = 2 + 2x + x^2 = q(x)$. $\qquad\blacksquare$

**Remark.** Using $h_1^2 + x e_2 = h_2 + (1+x)e_2$, equivalently $A_\lambda(x) = P_\lambda(1+x)$ where
$P_\lambda(t) = \langle s_\lambda,(h_2+te_2)^m\rangle = \sum_k \binom mk \nu_k t^k$,
$\nu_k = \langle s_\lambda, h_2^{m-k}e_2^k\rangle \ge 0$. Then $G_\lambda(i) = P_\lambda(i)$ and
$$G_\lambda(i) = 0 \iff (t^2+1)\mid P_\lambda(t), \qquad P_{(2,2)}(t) = t^2 + 1.$$
The $\omega$-involution gives $\nu_k(\lambda) = \nu_{m-k}(\lambda')$, hence
$P_{\lambda'}(t) = t^m P_\lambda(1/t)$ (the divisibility is conjugation-stable, as $\pm i$ are
mutually reciprocal; $(2,2)$ is self-conjugate).

**This is the cleanest face of the problem.** $(2,2)$ is *the unique partition whose
$A_\lambda$ equals the minimal polynomial $q$ itself.* For every other $\lambda$ the question is
whether the Eisenstein-at-$2$ quadratic $q$ divides the nonnegative-integer polynomial $A_\lambda$.

---

## 2. A real reformulation (PROVED)

Grouping the terms of $G_\lambda(i)=\sum_j c_j(i-1)^j$ (with $c_j=\binom mj M_j$) by
$j \bmod 4$, using $(i-1)^4 = -4$:

**Proposition 2.** Let $S_r = \sum_{s\ge 0}\binom{m}{4s+r} M_{4s+r}\,(-4)^s \in \mathbb{Z}$ for
$r\in\{0,1,2,3\}$. Then
$$\operatorname{Re} G_\lambda(i) = S_0 - S_1 + 2S_3, \qquad \operatorname{Im} G_\lambda(i) = S_1 - 2S_2 + 2S_3 .$$
Hence $G_\lambda(i)=0 \iff S_0 = 2(S_2-2S_3)$ and $S_1 = 2(S_2-S_3)$ (four integers, two
relations). In particular vanishing forces $S_0,S_1$ even; reducing mod $2$ gives
$S_0\equiv f^\lambda$ and $S_1\equiv m M_1\pmod 2$, recovering the **parity criterion**
$f^\lambda$ odd $\Rightarrow G_\lambda(i)\neq 0$.

*(Verified exact against $\operatorname{Re},\operatorname{Im}$ of $G$ for all $\lambda$, $m\le 6$.)*

**Remainder form.** Equivalently, the remainder $R(x)=r_1 x + r_0$ of $A_\lambda$ mod $q$ has
$r_1 = \operatorname{Im} G_\lambda(i)$ and $r_0 = \operatorname{Re} G_\lambda(i)+\operatorname{Im} G_\lambda(i)$.

---

## 3. An equivalent elementary criterion (PROVED), and why it gives no free lunch

**Proposition 3.** For every integer $n$, $A_\lambda(n) \equiv R(n) \pmod{q(n)}$. Hence
$$G_\lambda(i)=0 \iff q(n)\mid A_\lambda(n) = P_\lambda(n+1) \text{ for all } n\in\mathbb{Z}.$$
Each such failure is a one-line non-vanishing certificate. E.g. $n=0$: $2\mid f^\lambda$;
$n=1$: $5\mid \langle s_\lambda,(h_2+2e_2)^m\rangle$; in general the primes that can divide
$q(n)=(n{+}1)^2+1$ are exactly $2$ and primes $\equiv 1 \pmod 4$.

*Proof.* $A_\lambda(n)-R(n)=q(n)B(n)$ with $B\in\mathbb{Z}[x]$. For the "$\Leftarrow$": if
$q(n)\mid A_\lambda(n)$ for all $n$ then $q(n)\mid R(n)=r_1n+r_0$ for all $n$; since
$\deg R\le 1 < 2 = \deg q$ and $q(n)\to\infty$, this forces $R\equiv 0$, i.e. $q\mid A_\lambda$. $\blacksquare$

**Honest caveat.** Because $q(n)\mid r_1 n + r_0$ for *all* $n$ is equivalent to $r_0=r_1=0$,
this criterion is *logically equivalent* to $G_\lambda(i)=0$ — it does not bypass the problem.
Empirically (below) a *tiny* witness $n$ suffices, but the required bound is **not uniform in
$m$**, so finitely many $n$ cannot close the theorem.

---

## 4. The $(1+i)$-adic structure: what is proved, and the exact obstruction

Write $\pi = 1+i$ (so $v_\pi(N)=2v_2(N)$ for $N\in\mathbb{Z}$, residue field $\mathbb{F}_2$), and
for each $j$ with $M_j\neq 0$ set
$$\mathrm{val}(j) = v_\pi\!\big(c_j (i-1)^j\big) = j + 2\,v_2(c_j), \quad c_j = \binom mj M_j,
\qquad \mu = \min_j \mathrm{val}(j),\quad J^* = \arg\min_j \mathrm{val}(j).$$

**Fact 4.1 (single parity class).** $\mathrm{val}(j{+}1)-\mathrm{val}(j) = 1 + 2(v_2(c_{j+1})-v_2(c_j))$
is odd, so all elements of $J^*$ share one parity. *(Immediate.)*

**Theorem 4.2 (exact leading coefficient; PROVED).** The coefficient of $\pi^\mu$ in
$G_\lambda(i)$ is the unit-sum
$$w = \sum_{j\in J^*} o_j\, i^{(3\mu - j)/2}, \qquad o_j = \text{odd part of } c_j,$$
and $w \equiv |J^*| \pmod \pi$. Therefore:
$$|J^*| \text{ odd} \implies w\not\equiv 0 \implies v_\pi(G_\lambda(i)) = \mu < \infty
\implies G_\lambda(i)\neq 0.$$
*(Derivation: $c_j(i-1)^j = o_j(-i)^{v_2(c_j)} i^j \pi^{\mathrm{val}(j)}$ and
$v_2(c_j)=(\mu-j)/2$ on $J^*$ give $w$; units in $\mathbb{Z}[i]$ are $\equiv 1\bmod\pi$, so
$w\equiv|J^*|$. Verified: $0$ parity mismatches, $m\le 7$.)*

This re-proves non-vanishing for the **72%** of shapes with $|J^*|=1$, in clean language.

**Pair depths (PROVED).** If $J^*=\{a,a+2\}$ then $w \propto o_{a+2}+i\,o_a$ and
$v_\pi(w)=v_2(o_{a+2}^2+o_a^2)=1$ (sum of two odd squares is $\equiv 2 \bmod 8$). If
$J^*=\{a,a+4\}$ then $w\propto o_a-o_{a+4}\in 2\mathbb{Z}$, so $v_\pi(w)\ge 2$. Thus a tying
pair contributes $v_\pi(w)$ depending on the index gap mod $4$.

### The obstruction, precisely

When $|J^*|$ is **even** the leading coefficient cancels ($w\equiv 0$) and one must pass to
higher order. The difficulty is genuine and unbounded:

- Split $G = G_E + G_O$ into even-$j$ and odd-$j$ parts. Each term of $G_E$ has even
  $\mathrm{val}$, each term of $G_O$ odd; but $G_E$ is itself a *sum of units that can cancel to
  arbitrary $\pi$-depth*, and likewise $G_O$. The actual $v_\pi(G)$ is the minimum over an
  **interleaved two-class cancellation tower**, not over the leading pair.
- **Worked example $\lambda=(3,3,1,1)$, $m=4$:** $J^*=\{0,4\}$ (even), $\mu=6$; the even pair
  cancels to $\pi$-depth $\mu + v_\pi(w) = 6+4 = 10$. But the odd class has its minimum at
  $\mathrm{val}=7=\mu+1$ (indices $\{1,3\}$), and *that* pair cancels to depth $8$. Hence
  $v_\pi(G_{(3,3,1,1)}) = 8$, governed by the **odd** class — even though the even class held
  the global minimum $\mu$. The leading-pair analysis alone gives the wrong answer.
- Each rung $G_E=\tilde A(-2i)$ (with $\tilde A(y)=\sum_l c_{2l}y^l\in\mathbb{Z}_{\ge0}[y]$) is a
  *new* divisibility problem (now against $y^2+4$), so the descent is not self-similar; no fixed
  finite $\pi$-adic level closes it. This is why the maximal finite cancellation depth grows
  ($v_\pi^{\max} = 5,9,11,15,15,21,\dots$ for $n=6,\dots,16$) and a *uniform-in-$m$* argument is
  required.

---

## 5. Strengthened computational evidence (this session)

Exact $\mathbb{Z}[i]$ computation (sympy, power-sum basis, Murnaghan–Nakayama characters):

| $m$ | $n$ | #shapes | non-$(2,2)$ vanishers | $\lvert J^*\rvert$ distribution |
|----:|----:|--------:|:---------------------:|:-------------------------------:|
| 2 | 4 | 5 | none | $\{1{:}4,\,2{:}1\}$ |
| 3 | 6 | 11 | none | $\{1{:}9,\,2{:}2\}$ |
| 4 | 8 | 22 | none | $\{1{:}16,\,2{:}6\}$ |
| 5 | 10 | 42 | none | $\{1{:}33,\,2{:}9\}$ |
| 6 | 12 | 77 | none | $\{1{:}50,\,2{:}24,\,4{:}3\}$ |
| 7 | 14 | 135 | none | $\{1{:}98,\,2{:}35,\,4{:}2\}$ |
| 8 | 16 | 231 | none | $\{1{:}141,\,2{:}78,\,4{:}12\}$ |
| 9 | 18 | 385 | none | $\{1{:}264,\,2{:}110,\,4{:}11\}$ |
| 10 | 20 | 627 | none | $\{1{:}360,\,2{:}221,\,4{:}46\}$ |

- **$(2,2)$ is the unique vanisher through $m=10$ ($n=20$).** (Previously verified $n\le 16$.)
- **$|J^*| \in \{1,2,4\}$ always** — never odd $>1$, and (through $m=10$) never $8$: the
  minimizer count is observed to be a power of $2$, $\le 4$. This is a *finer* fact than mere
  evenness and is a strong hint of an underlying Kummer/Lucas (binary-carry) mechanism in
  $v_2\!\big(\binom mj M_j\big)$. The $|J^*|=4$ index-sets seen ($m\le 10$) are
  $\{0,2,4,6\}$, $\{3,5,7,9\}$ (four consecutive in a parity class) and $\{0,2,8,10\}$. The
  last, in $t=j/2$ coordinates $\{0,1,4,5\}=\{0,1\}\oplus\{0,4\}$, is a binary $2$-dimensional
  subcube — but $\{3,5,7,9\}\leftrightarrow t\in\{1,2,3,4\}$ is *not* a sumset subcube, so there
  is no uniform "subcube" description; only the power-of-$2$ cardinality is robust.
- Minimal integer witness $n$ in Prop. 3 grows with $m$ ($\max|n| = 3$ for $m\le 7$, but $=5$ at
  $m=8,9$): no uniform finite witness set exists.

---

## 6. Status, conjectures, and the precise gap

**Proved (rigorous, this session):**
1. Theorem 1 — $G_\lambda(i)=A_\lambda(i-1)$; $G=0 \iff q\mid A_\lambda$; $A_{(2,2)}=q$.
2. Propositions 2, 3 — real reformulation and the equivalent (but non-uniform) integer criterion.
3. Theorem 4.2 — exact leading $\pi$-coefficient; $|J^*|$ odd $\Rightarrow G\neq 0$ (the $\ge 72\%$).
4. Pair-depth formulas (§4).

**Conjectures (strong evidence, $m\le 10$):**
- **(C1)** $G_\lambda(i)=0 \iff \lambda=(2,2)$ (the theorem).
- **(C2, evenness)** For $\lambda\neq(2,2)$, $|J^*|$ is even — sharper: $|J^*|\in\{2,4\}$, a power
  of $2$. *Status: tested to $m=10$, no structural proof. Note: (C2) does not by itself imply
  non-vanishing — it confirms the leading-order argument cancels, leaving the tower of §4.*

**The gap, precisely stated.** Both (C1) and (C2) reduce to controlling the
**interleaved even/odd $(1+i)$-adic cancellation tower** of §4 *uniformly in $m$*. Concretely:

> **(Gap)** Show that for $\lambda\neq(2,2)$ the tower terminates — i.e. $v_\pi(G_\lambda(i))<\infty$
> — equivalently that the polynomial $q(x)=x^2+2x+2$ does **not** divide
> $A_\lambda(x)=\langle s_\lambda,(h_1^2+xe_2)^m\rangle$. The cancellation depth is unbounded in
> $m$, so no fixed $\pi$-adic level (and no fixed finite set of integer witnesses, §3, §5)
> suffices; a uniform mechanism is needed.

**Most promising handles identified this session:**
- The power-of-$2$ structure of $|J^*|$ suggests writing $v_2\!\big(\binom mj M_j\big)$ via
  Kummer ($v_2\binom mj = s_2(j)+s_2(m{-}j)-s_2(m)$) plus a binary model for $v_2(M_j)$, then
  proving the min-locus of $j + 2v_2(c_j)$ in one parity class has power-of-$2$ cardinality by a
  carry symmetry.
- The clean polynomial form $q\mid A_\lambda$ (Theorem 1) may admit a coefficient/positivity
  argument unavailable in the $\pi$-adic picture: $A_\lambda\in\mathbb{Z}_{\ge0}[x]$ with
  $A_\lambda(n)>0$ for $n\ge 0$, $A_\lambda(0)=f^\lambda$, $A_\lambda(-1)=\langle s_\lambda,h_2^m\rangle$.

---

## 7. Reproducibility

Scripts (sympy, exact): `~/projects/scratch/2026-06-10-oneplusi-adic/` (machinery, valuations),
and this session's `srcheck.py` (Prop. 2), `data_jstar.py` (within-class structure),
`verify_big.py` (Table §5), `nec_conditions.py`/`min_witness.py` (Prop. 3 witnesses),
`leadcoef_check.py` (Theorem 4.2). Working notebook:
`~/projects/scratch/prove-2026-06-09-evenness.md`.
