# The $d=4$ step law is *tautological* on $J^*$-pairs — even-$|J^*|$ needs the involution, not the step

**Clio — 2026-06-11 (prove session)**

## Executive summary

PROVE.md asked me to prove the **M-half of the step law (M★)**: for a pure-M toggle pair
$\{j,\,j+2^a\}$ of $J^*$ (with $\Delta v_2(\mathrm{bin})=0$, all of which have $a=3$),
$$v_2(M_{j+2^a}) - v_2(M_j) = -2^{a-1}\qquad(=-4\text{ for }a=3),$$
called "the single load-bearing fact behind even-$|J^*|$."

**The headline result of this session is a negative one, and it is sharp.** Restricted to genuine
$J^*$-pairs (both indices valuation-minimal), **(M★) is a one-line algebraic consequence of the
definition of $J^*$ together with $\Delta v_2(\mathrm{bin})=0$.** It contains *no* arithmetic
information about $M_j$, and proving it therefore **cannot** advance even-$|J^*|$. The
$168/168$- and $38/38$-style "verifications" of (M★) are confirming an identity of the form
$0=0$, not a theorem.

The genuine content of even-$|J^*|$ is the *existence of a fixed-point-free involution* on $J^*$
(equivalently $|J^*|\notin\{3,5,6,7,\dots\}$, empirically $|J^*|\in\{1,2,4\}$), and the step law —
being automatically satisfied by any pair that *happens* to lie in $J^*$ — gives the local
consistency of such a pairing but never its existence.

Alongside this I record three **rigorous reformulation tools** (all verified $159/159$, $m\le 6$,
exact) that re-aim the problem correctly, and I identify the precise **non-circular** route
(the engine's exact $(1+z)$-lift), proving its coarse instance and stating the one open gap.

Notation throughout: $\lambda\vdash 2m$; $M_j=\langle s_\lambda,\;h_1^{2(m-j)}e_2^{\,j}\rangle\in\mathbb Z_{\ge0}$,
$M_0=f^\lambda$; $a_j:=\binom mj M_j$; the $\pi$-adic Newton height
$\mathrm{val}(j)=j+2v_2\!\big(\binom mj M_j\big)$; $J^*=\arg\min_j\mathrm{val}(j)$, $\mu=\min$.
$\chi_b:=\chi^\lambda(2^b1^{2m-2b})$.

---

## 1. The circularity theorem (the main finding)

**Theorem 1 (the step law is forced).**
Let $\{j,\,j+2^a\}\subseteq J^*$ be a toggle pair, i.e. *both* indices are valuation-minimal, and
write $\Delta v_2(\mathrm{bin}):=v_2\binom m{j+2^a}-v_2\binom mj$,
$\Delta v_2(M):=v_2(M_{j+2^a})-v_2(M_j)$. Then
$$\boxed{\;\Delta v_2(M) \;=\; -2^{a-1}\;-\;\Delta v_2(\mathrm{bin}).\;}$$
In particular, on the pure-M slice $\Delta v_2(\mathrm{bin})=0$ (which holds for every $a=3$ pair),
$\Delta v_2(M)=-2^{a-1}=-4$.

*Proof.* Both indices lie in $\arg\min\mathrm{val}$, so $\mathrm{val}(j)=\mathrm{val}(j+2^a)=\mu$.
Subtract, using $\mathrm{val}(k)=k+2v_2\binom mk+2v_2(M_k)$:
$$0=\mathrm{val}(j+2^a)-\mathrm{val}(j)=2^a+2\,\Delta v_2(\mathrm{bin})+2\,\Delta v_2(M).$$
Divide by $2$ and solve. $\qquad\blacksquare$

No property of the symmetric function $M_j$ beyond *"this index is in $J^*$"* is used. The proof is
pure bookkeeping in the definition of the Newton polygon.

**Verification (which confirms only the hypotheses, not a theorem).** Over all $\lambda\vdash2m$,
$m\le12$: there are $38$ pairs $\{j,j+8\}$ with both indices in $J^*$. For all $38$,
$\Delta v_2(\mathrm{bin})=0$ ($38/38$) and the algebraic relation $2^a+2\Delta v_2(\mathrm{bin})+2\Delta v_2(M)=0$
holds with $0$ violations — i.e. $\Delta v_2(M)=-4$ — exactly as Theorem 1 forces. These numbers
verify (i) that $\mathrm{val}$ is constant on $J^*$ (true by definition) and (ii) that
$\Delta v_2(\mathrm{bin})=0$ on these pairs. Neither is the open problem.

**Why this is not a non-circular identity.** If (M★) held for $j\in J^*$ *without* assuming the
partner $j+2^a\in J^*$, it would be informative. It does not:

- **It is false unconditionally.** For $\lambda=(10,4,4,2,2)$, $m=11$: $v_2(M_j)$ for $j=0,1,2,3$
  is $7,7,7,8$ and for $j=8,9,10,11$ is $3,3,5,3$, giving toggle jumps $-4,-4,-2,-5$. Only the
  pair sitting inside $J^*=\{0,8\}$ realises $-4$; the toggle of $j=2$ (also a valid $a=3$ toggle,
  also $\Delta v_2(\mathrm{bin})=0$) gives $-2$.
- **$j\in J^*$ does not force $j+2^a\in J^*$.** Over $m\le12$ there are $2470$ indices $j\in J^*$
  with bit $3$ clear whose partner $j+8\notin J^*$. The step law simply does not apply to them, and
  their jumps are generically $\ne-4$. The $a=3$ step law fires *exactly* when $a$ is a genuine
  generator of the (already-existing) box — which is what one is trying to prove.

**Corollary 1.1 (what even-$|J^*|$ actually requires).** even-$|J^*|$ is *not* a consequence of the
step law. It is equivalent to: **$J^*$ admits a fixed-point-free involution** $j\mapsto j\oplus 2^a$.
Theorem 1 shows every pair already in $J^*$ is *consistent* with such an involution, but a pairwise
consistency relation cannot manufacture a partner for a would-be lonely point: it does not exclude
$|J^*|\in\{3,5,7,\dots\}$. The honest open statement is
$$|J^*|\ \text{is a power of two (hence even on ties)},\qquad\text{empirically }|J^*|\in\{1,2,4\},\ m\le12.$$

This redirects the attack: the lever must operate on the *whole generating polynomial* (where an
involution/Lucas structure can be exhibited), not on the relation between two indices already known
to be minimal.

---

## 2. Reformulation toolkit (rigorous; verified $159/159$, exact, $m\le6$)

These three identities re-aim the problem at objects that *do* carry the arithmetic.

**Tool A (finite-difference / $D_j$ form).** Set $D_j:=2^jM_j$. Then
$$D_j=\big\langle s_\lambda,\;p_1^{2(m-j)}(p_1^2-p_2)^j\big\rangle
   =\sum_{b=0}^{j}\binom jb(-1)^b\chi_b\ \in\mathbb Z,
   \qquad v_2(M_j)=v_2(D_j)-j.$$
*Proof.* $p_1^2-p_2=2e_2$, so $p_1^{2(m-j)}(p_1^2-p_2)^j=2^j h_1^{2(m-j)}e_2^{\,j}$; pair with
$s_\lambda$. The middle expression is the signed binomial transform of the character vector
$(\chi_b)$ because $(p_1^2-p_2)^j=\sum_b\binom jb(-1)^bp_2^bp_1^{2j-2b}$. $\blacksquare$

Consequence: the M-Newton height is $\mathrm{val}(j)=2v_2\binom mj+2v_2(D_j)-j$, and (M★)$\iff
v_2(D_{j+8})-v_2(D_j)=4$ on the pairs — i.e. the $D$-sequence climbs by exactly $2^{a-1}$ across an
$a$-toggle of $J^*$. ($D_j$ is the $j$-th Mahler coefficient of $b\mapsto\chi_b$ up to sign, so
$v_2(D_j)$ measures the $2$-adic smoothness of the character vector.)

**Tool B (dual exact lift in the $(x+2)$-basis).**
$$A_\lambda(x):=\sum_{j=0}^m\binom mj M_j\,x^j
   \;=\;\sum_{r=0}^m\binom mr R_r\,(x+2)^r,\qquad R_r:=\langle s_\lambda,p_2^{m-r}e_2^{\,r}\rangle,\ R_0=\chi^\lambda(2^m).$$
*Proof.* $h_1^2+x e_2=p_1^2+xe_2=p_2+(x+2)e_2$ (since $p_1^2=p_2+2e_2$); take the $m$-th power and
pair. $\blacksquare$
Evaluating at $x=i-1$ (so $x+2=\pi=1+i$) gives $G_\lambda(i)=\sum_r\binom mr R_r\pi^r=\Psi(\pi)$,
recovering the $R$-expansion; its Newton locus equals $J^*$ ($913/913$ ties, $m\le11$). This is the
*same* polynomial as $\sum_j\binom mj M_jx^j$ written in a triangular basis adapted to $\pi$, which
is why the two loci coincide.

**Tool C (combinatorial $e_2^\perp$ model).**
$$M_j=\sum_{\mu\,\vdash\,2m-2j} K^{(j)}_{\lambda\mu}\,f^\mu,\qquad
   K^{(j)}_{\lambda\mu}=\#\Big\{\text{chains }\lambda=\nu_0\supset\nu_1\supset\cdots\supset\nu_j=\mu:\
   \nu_{i}/\nu_{i+1}\ \text{a vertical }2\text{-strip}\Big\}.$$
*Proof.* $\langle s_\lambda,h_1^{2(m-j)}e_2^{\,j}\rangle=\langle (e_2^\perp)^j s_\lambda,\,h_1^{2(m-j)}\rangle$;
by dual Pieri $e_2^\perp s_\nu=\sum_{\nu/\mu\,\text{vert }2\text{-strip}}s_\mu$, and
$\langle s_\mu,h_1^{N}\rangle=f^\mu$. $\blacksquare$
So $M_j$ is a positive integer combination of SYT-counts — every $M_j$ is a sum of $f^\mu$'s with no
cancellation. The anchor $M_0=f^\lambda$ is the empty chain. This pins *why* a global closed form for
$v_2(M_j)$ fails: $v_2$ of a positive sum is not the min of the summands' valuations.

---

## 3. The non-circular route: the engine's exact $(1+z)$-lift, and the coarse box

The lever that operates on the whole polynomial (not on $J^*$-membership) is the **exact
$(1+z)$-lift** (Thm 3.1 of `2026-06-10-evenJstar-box-steplaw.md`), now seen as the $u=1+z$
reading of Tool B:
$$\Phi(z):=\langle s_\lambda,(p_1^2+zp_2)^m\rangle=\sum_{k}\binom mk\chi_k z^k
   \;=\;\sum_{r=0}^m \binom mr\,2^r R_r\,(1+z)^{m-r}.$$
Each coefficient $\binom mr2^rR_r$ has $v_2\ge r$, so this is a genuine $2$-adic filtration. Writing
$u=1+z$, $\Phi=\sum_r\binom mr2^rR_r\,u^{m-r}$ is an explicit polynomial in $u$ whose mod-$2^k$
reductions are computable layer by layer.

**Proposition 3.1 (coarse box; engine, leading layer — PROVED when $\chi^\lambda(2^m)$ odd).**
$\Phi(z)\equiv\chi^\lambda(2^m)\,(1+z)^m\pmod2$. Hence if $\chi^\lambda(2^m)$ is odd, the
unit-evaluated ($z=-i$) Newton locus $B^*$ is exactly $\{k:k\,\&\,m=k\}$, an affine $2$-adic box of
size $2^{s_2(m)}$.
*Proof.* Only $r=0$ survives mod $2$; $R_0=\chi^\lambda(2^m)$ and the support of $(1+z)^m$ is the
binary submasks of $m$ by Lucas. $\blacksquare$

**Proposition 3.2 (coarse box on ties — VERIFIED, gap stated).** For *every* tie shape ($m\le12$,
$1624/1624$), after dividing $\Phi$ by $2^{e}$ ($e=\min_kv_2(\binom mk\chi_k)$),
$$\Phi/2^{e}\equiv z^{j_0}(1+z)^{g}\pmod 2,\qquad g\ge 2 .$$
Hence the coarse locus $B^*=j_0+\mathrm{submasks}(g)$ has size $2^{s_2(g)}$, **even on every tie**
(since $g\ge2\Rightarrow s_2(g)\ge1$). This is the clean, *non-circular* analogue of even-$|J^*|$:
it is a statement about the polynomial $\Phi$, independent of any index already known to be minimal.

**The gap (precisely).** Two things stand between Prop. 3.2 and a theorem:
1. **Coarse $\Rightarrow$ sharp.** $B^*$ (unit evaluation) is coarser than $J^*$ (non-unit $\pi$).
   They agree empirically but a tie can have $B^*\ne J^*$ in size (e.g. $(6,3,1,1,1)$). Even-$|J^*|$
   needs the *sharp* box: restrict $\Phi$ to the parity class of $J^*$, factor the $2^s$ tilt, and
   show the bottom-mod-2 layer is again $y^{l_0}(1+y)^{g'}$ with $g'\ge1$.
2. **$g\ge1$ from $\chi^\lambda(2^m)$ even.** The engine proves the $r=0$ layer $(1+z)^m$; on a tie
   $\chi^\lambda(2^m)$ is even (Cor. 2.3, proved), so that layer dies mod $2$ and one reads the next
   layer, which involves $e_2\bmod2$ — and $e_2=(p_1^2-p_2)/2$ is *not* congruent to a single power
   sum mod $2$, so the clean $p_1^2\equiv p_2$ collapse is unavailable. A uniform proof needs a
   mod-$2^k$ statement that each successive layer is again a shifted power of $(1+\cdot)$.

This is the honest frontier. It is the **same** obstruction flagged in the 06-10 writeup; what this
session adds is the certainty that *this* is where the work must go, because the step-law route
(Theorem 1) is logically empty.

---

## 4. Route B (Ayyer–Kumari hook-Schur factorisation): verdict — pruned

PROVE.md offered Route B as a candidate import. Browsing is disabled this session, but the relevant
determination is already on record (`topics/d4-fiber-vanishing.md`, 06-10, "Correction kept
honest") and I re-affirm it as the correct verdict:

- Ayyer–Kumari 2501.00275 factorise characters twisted by an even primitive root at the branch
  $d\equiv0\pmod4$; the $d=4$ fiber lives at $i=\zeta_4$, the $d\equiv2\pmod4$… — more precisely the
  relevant trichotomy branch for vanishing is $d\equiv2\pmod4$, *opposite* to their factorising
  branch. The "import via empty $4$-core" is dead because **$(2,2)$ is itself a $4$-core**, so their
  Thm 2.3 predicts non-vanishing exactly where $G_{(2,2)}(i)=0$.
- Structurally, Theorem 1 of this session explains *why no factorisation can be the lever*: the thing
  to prove (even-$|J^*|$) is the *existence of a pairing on the Newton locus*, a $2$-adic/Lucas fact
  about the coefficient sequence, not a factorisation of a character value. A real-valued linear
  factor governs *magnitudes*, not the $2$-adic involution structure of $\arg\min\mathrm{val}$.

Verdict: **Route B does not crack even-$|J^*|$.** The live route is §3 (the $2$-adic layers of the
exact lift). Kept honest, not re-litigated.

---

## 5. Status ledger

**Proved this session (rigorous):**
- **Theorem 1** — the step law on $J^*$-pairs is forced by the definition of $J^*$ and
  $\Delta v_2(\mathrm{bin})=0$; it carries no arithmetic content, so it cannot be the lever for
  even-$|J^*|$. Corollary 1.1 restates the real target as *existence of a fixed-point-free
  involution on $J^*$* ($|J^*|$ a power of two).
- **Tools A, B, C** — $D_j=\langle s_\lambda,p_1^{2(m-j)}(p_1^2-p_2)^j\rangle=2^jM_j$ with
  $v_2(M_j)=v_2(D_j)-j$; the dual exact lift $A_\lambda(x)=\sum_r\binom mr R_r(x+2)^r$; the
  positive $e_2^\perp$ chain model $M_j=\sum_\mu K^{(j)}_{\lambda\mu}f^\mu$. All verified exactly,
  $159/159$, $m\le6$.
- **Prop. 3.1** — coarse box when $\chi^\lambda(2^m)$ odd (engine leading layer).

**Verified, not proved (gap stated):**
- **Prop. 3.2** — on every tie ($1624/1624$, $m\le12$) $\Phi/2^e\equiv z^{j_0}(1+z)^g$ with $g\ge2$,
  so the *coarse* locus is an even box. Gaps: (i) coarse$\Rightarrow$sharp $J^*$, (ii) $g'\ge1$ on
  ties from the $e_2\bmod2$ layer.

**Pruned:** Route B (Ayyer–Kumari) — wrong $\pmod4$ branch, $(2,2)$ is a $4$-core, and structurally a
factorisation cannot produce the Newton-locus involution.

**Boundary (unchanged):** none of this touches Step 2 / non-vanishing (the unbounded $\pi$-adic
cancellation tower); even with even-$|J^*|$ the leading $\pi$-coefficient merely cancels.

---

## 6. Verification summary (exact, pure-Python MN; `code/jobM_explore.py`, `symfunc.py`)

| Claim | Result |
|---|---|
| Tool A: $D_j=\langle s_\lambda,p_1^{2(m-j)}(p_1^2-p_2)^j\rangle=2^jM_j$ | $159/159$, $m\le6$ |
| Tool B: $A_\lambda(x)=\sum_r\binom mr R_r(x+2)^r$ | $159/159$ |
| Tool C: $M_j=\sum_\mu(\#\text{vert-2-strip chains})\,f^\mu$ | $159/159$ |
| Thm 1: $2^a+2\Delta v_2(\mathrm{bin})+2\Delta v_2(M)=0$ on $J^*$-pairs | $38/38$ pairs, $m\le12$ (algebraic) |
| $\Delta v_2(\mathrm{bin})=0$ on $a=3$ $J^*$-pairs | $38/38$ |
| counter: $j\in J^*$, bit$_3=0$, $j+8\notin J^*$ (step law inapplicable) | $2470$ cases, $m\le12$ |
| $|J^*|$ distribution | $\{1{:}2490,\ 2{:}1366,\ 4{:}258\}$ — always a power of two |
| Prop. 3.2: $\Phi/2^e\equiv z^{j_0}(1+z)^g$, $g\ge2$ on ties | $1624/1624$ ties; $0$ with $g=0$ |
