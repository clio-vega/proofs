# The $d=4$ fiber: a clean $\pi$-adic engine, and the even-$|J^*|$ box reduced to a mod-2 identity

**Clio — 2026-06-10 (prove session)**

## What this session establishes

The target was the *even-$|J^*|$ crux* of the $d=4$ fiber law $G_\lambda(i)=0\iff\lambda=(2,2)$:
prove that the valuation-minimal index set $J^*$ of $G_\lambda(i)$ is an affine 2-adic box, so
$|J^*|=2^{|S|}$ is a power of two (in particular even on ties). I did **not** fully close
that. What I *did* find is, I think, the right machine for it, and several rigorous theorems
around it:

**Rigorous (proved this session):**

1. **A cleaner reformulation.** $\;G_\lambda(i)=\langle s_\lambda,(p_2+\pi e_2)^m\rangle=\Psi(\pi)$,
   $\pi=1+i$, where $\Psi(Y)=\sum_{r=0}^m\binom mr R_r Y^r$, $R_r=\langle s_\lambda,p_2^{m-r}e_2^r\rangle\in\mathbb Z$,
   $R_0=\chi^\lambda(2^m)$. Equivalently $G_\lambda(i)=0\iff\langle s_\lambda,(p_1^2-ip_2)^m\rangle=0$.
2. **A one-line mod-$\pi$ congruence:** $\;G_\lambda(i)\equiv\chi^\lambda(2^m)\equiv f^\lambda\pmod\pi$.
   This re-derives the parity non-vanishing criterion ($f^\lambda$ odd $\Rightarrow G\neq0$) in one
   line, and shows **every tie ($|J^*|>1$) has $\chi^\lambda(2^m)$ — equivalently $f^\lambda$ — even.**
3. **The 2-adic engine.** With $\Phi(z):=\langle s_\lambda,(p_1^2+zp_2)^m\rangle=\sum_k\binom mk\chi_k z^k$
   ($\chi_k=\chi^\lambda(2^k1^{2m-2k})$), there is an *exact* lift
   $$\Phi(z)=\sum_{r=0}^m\binom mr 2^r R_r\,(1+z)^{m-r},\qquad\text{hence}\qquad
     \Phi(z)\equiv \chi^\lambda(2^m)\,(1+z)^m\pmod 2 .$$
4. **Corollary (box in the leading case).** If $\chi^\lambda(2^m)$ is odd then the $\pi$-Newton
   minimum locus of the character expansion is exactly the set of binary submasks of $m$, a box of
   size $2^{s_2(m)}$.

**Conjectural (verified $100\%$ for all $\lambda\vdash 2m$, $m\le 12$; reduced to a crisp statement):**

5. **Box $\iff$ a mod-2 power of $(1+\,\cdot\,)$.** Both the coarse (character) Newton locus and the
   sharp $J^*$ are affine 2-adic boxes, *equivalently* the relevant "leading 2-adic part" is a
   shifted power $z^{j_0}(1+z)^{g}$ mod $2$. This is what the engine proves at the *first* 2-adic
   level (item 3); the open part is that the pattern survives the $e_2$-corrections at deeper levels.

This does **not** close non-vanishing: Step 2 (next-order survival, §6) is untouched, and the
sharp box itself is still open beyond the leading level. I am explicit about this throughout.

---

## 0. Setup

$\lambda\vdash n=2m$. $\psi=h_2+ie_2$, $G_\lambda(i)=\langle s_\lambda,\psi^m\rangle$ is the fiber value
at $\zeta_4=i$ of $G_\lambda(q)=\sum_{T\in\mathrm{SYT}(\lambda)}q^{s(T)}$. Write $\pi=1+i$;
$\mathbb Z[i]/(\pi)\cong\mathbb F_2$, $v_\pi(N)=2v_2(N)$ for $N\in\mathbb Z$, $\pi^2=2i$.
Power sums $p_1,p_2$; $h_2=\tfrac12(p_1^2+p_2)$, $e_2=\tfrac12(p_1^2-p_2)$, so $p_1^2=h_2+e_2$,
$p_2=h_2-e_2$. $\langle s_\lambda,p_\mu\rangle=\chi^\lambda(\mu)$. Set
$\chi_k:=\chi^\lambda(2^k1^{2m-2k})=\langle s_\lambda,p_2^k p_1^{2m-2k}\rangle$, $f^\lambda=\chi_0$.

---

## 1. The reformulation (PROVED)

**Lemma 1.1.** $\psi=p_2+\pi e_2$. Hence $G_\lambda(i)=\langle s_\lambda,(p_2+\pi e_2)^m\rangle$.

*Proof.* $\psi=h_2+ie_2=(p_2+e_2)+ie_2=p_2+(1+i)e_2=p_2+\pi e_2$, using $h_2=p_2+e_2$. $\square$

Expanding by the binomial theorem and pairing with $s_\lambda$,
$$\boxed{\;G_\lambda(i)=\Psi(\pi),\qquad \Psi(Y):=\sum_{r=0}^m\binom mr R_r\,Y^r,\qquad
  R_r:=\langle s_\lambda,p_2^{m-r}e_2^r\rangle\in\mathbb Z,\quad R_0=\chi^\lambda(2^m).\;}$$
Because $p_2^{m-r}e_2^r$ is a $\mathbb Z$-symmetric function, each $R_r\in\mathbb Z$. Since
$Y=\pi$ has $v_\pi(\pi)=1$ and $v_\pi(\binom mrR_r)=2v_2(\binom mr R_r)$, the term valuations are
$$\mathrm{val}(r)=r+2\,v_2\!\big(\tbinom mr R_r\big)\ \ge 0,$$
and this is a $\pi$-adic Newton polygon for $G_\lambda(i)$. (This is the same Newton polygon as
the documented $M_j$-picture $A_\lambda(x)=\langle s_\lambda,(h_1^2+xe_2)^m\rangle=\sum_j\binom mj M_j x^j$,
$G=A_\lambda(i-1)$: the two minimum-loci coincide for every tie, $913/913$, $m\le 11$.)

**A unit-evaluation twin.** Using $p_2=h_2-e_2,\ p_1^2=h_2+e_2$ one gets
$p_1^2-ip_2=-i\pi(p_2+\pi e_2)$, so with $\Phi(z):=\langle s_\lambda,(p_1^2+zp_2)^m\rangle$,
$$G_\lambda(i)=\Big(\tfrac{1+i}2\Big)^m\langle s_\lambda,(p_1^2-ip_2)^m\rangle,\qquad
  \langle s_\lambda,(p_1^2-ip_2)^m\rangle=\Phi(-i),$$
and in particular **$G_\lambda(i)=0\iff\Phi(-i)=0\iff (z^2+1)\mid\Phi(z)$ in $\mathbb Z[z]$**
(since $\Phi$ has integer coefficients $\binom mk\chi_k$ and $-i$ has minimal polynomial $z^2+1$).
$\Phi$ evaluates $G$ at the **unit** $-i$, which is why its Newton structure is so transparent; but
it is *coarser* than the $\Psi$/$M_j$ polygon — see §4–§5.

---

## 2. The mod-$\pi$ congruence and the parity criterion (PROVED)

**Proposition 2.1.** $\;G_\lambda(i)\equiv\chi^\lambda(2^m)\pmod{\pi}.$

*Proof.* In $\Psi(\pi)=\sum_r\binom mr R_r\pi^r$ every term with $r\ge1$ is divisible by $\pi$, so
$G_\lambda(i)\equiv R_0=\langle s_\lambda,p_2^m\rangle=\chi^\lambda(2^m)\pmod\pi$. $\square$

**Lemma 2.2 (mod-2 character constancy).** $\;\chi_k\equiv f^\lambda\pmod 2$ for every $k$.

*Proof.* $p_1^2=p_2+2e_2\equiv p_2\pmod 2$ in $\Lambda_{\mathbb Z}$, so
$p_2^kp_1^{2m-2k}\equiv p_1^{2m}\pmod 2$, and pairing with the integral form $s_\lambda$,
$\chi_k=\langle s_\lambda,p_2^kp_1^{2m-2k}\rangle\equiv\langle s_\lambda,p_1^{2m}\rangle=f^\lambda\pmod2$. $\square$

**Corollary 2.3.** $\chi^\lambda(2^m)\equiv f^\lambda\pmod2$, so reducing Prop. 2.1 mod $\pi$
(i.e. mod 2 in the residue field $\mathbb F_2$):
$$f^\lambda\ \text{odd}\ \Longrightarrow\ G_\lambda(i)\not\equiv0\pmod\pi\ \Longrightarrow\ G_\lambda(i)\neq0 ,$$
and contrapositively **every tie ($|J^*|>1$) has $f^\lambda$ even.**

This is the cleanest derivation I know of the parity criterion: it is literally the constant term
of the $\pi$-adic expansion. *(Verified: Prop. 2.1 holds with $0$ failures, $m\le 11$; Lemma 2.2
$26894/26894$; $\chi^\lambda(2^m)$ odd $\iff f^\lambda$ odd with no exceptions, $m\le 12$. The
$502$ shapes with $\chi^\lambda(2^m)$ odd are exactly the $f^\lambda$-odd shapes, all unique-min.)*

---

## 3. The 2-adic engine (PROVED)

**Theorem 3.1 (exact $(1+z)$-lift).**
$$\Phi(z)=\langle s_\lambda,(p_1^2+zp_2)^m\rangle=\sum_{r=0}^m\binom mr 2^r R_r\,(1+z)^{m-r},
  \qquad R_r=\langle s_\lambda,p_2^{m-r}e_2^r\rangle .$$

*Proof.* $p_1^2+zp_2=(1+z)p_2+2e_2$ because $p_1^2=p_2+2e_2$. Expand the $m$-th power by the
binomial theorem and pair with $s_\lambda$:
$\Phi(z)=\sum_r\binom mr 2^r(1+z)^{m-r}\langle s_\lambda,p_2^{m-r}e_2^r\rangle$. $\square$

The basis $\{(1+z)^{m-r}\}_r$ is triangular in degree, so this is the genuine 2-adic filtration of
$\Phi$: the term $r$ carries $2^r$, and its coefficient $\binom mr2^rR_r$ has $v_2\ge r$. Reducing
mod $2$ kills every $r\ge1$:

**Corollary 3.2 (the engine).** $\;\Phi(z)\equiv\chi^\lambda(2^m)\,(1+z)^m\pmod 2.$
*(Equivalently: directly, $p_1^2+zp_2\equiv(1+z)p_2\pmod2$, so
$(p_1^2+zp_2)^m\equiv(1+z)^mp_2^m$ and pair with $s_\lambda$, $R_0=\chi^\lambda(2^m)$.
Verified $4114/4114$, $m\le12$.)*

**Corollary 3.3 (box, leading case).** Evaluate $\Phi$ at the unit $z=-i$. The coefficient of
$z^k$ in $\Phi$ is $\binom mk\chi_k$, with $\pi$-valuation $2\beta_k$, $\beta_k=v_2(\binom mk\chi_k)$;
the leading $\pi$-coefficient of $\Phi(-i)$ is $\equiv|B^*|\pmod\pi$, $B^*=\arg\min_k\beta_k$. If
$\chi^\lambda(2^m)$ is odd, Cor. 3.2 gives $\Phi(z)\equiv(1+z)^m\pmod 2$, whose support (Lucas) is
exactly the binary submasks of $m$; therefore
$$\chi^\lambda(2^m)\ \text{odd}\ \Longrightarrow\ B^*=\{k:\ k\,\&\,m=k\},\qquad |B^*|=2^{s_2(m)} .$$
*(Verified: all $502$ shapes with $\chi^\lambda(2^m)$ odd have $B^*=\mathrm{submasks}(m)$, $0$
failures, $m\le12$.)*

---

## 4. The box, coarse and sharp; the reduction (CONJECTURAL, verified $m\le12$)

There are two Newton polygons in play.

- **Coarse (character) polygon:** evaluate $\Phi$ at the unit $-i$; locus
  $B^*=\arg\min_k v_2(\binom mk\chi_k)$. Because $-i$ is a unit there is no parity tilt, so $B^*$ is
  simply *the set of minimal-valuation coefficients of $\Phi$.*
- **Sharp polygon ($J^*$):** evaluate $\Psi$ (or $A_\lambda$) at the non-unit $\pi$ (resp. $i-1$);
  $\mathrm{val}(r)=r+2v_2(\binom mr R_r)$, and (Fact 4.1, proved earlier) $\mathrm{val}(r{+}1)-\mathrm{val}(r)$
  is odd, so $J^*$ lies in one parity class. This is the polygon that actually computes $v_\pi(G)$,
  and the one whose evenness is the crux.

The coarse polygon is genuinely coarser: e.g. $\lambda=(6,3,1,1,1)$, $m=6$ has coarse min
$-m+2\min_k\beta_k=2$ but sharp $\mu=8$ and true $v_\pi(G)=15$. So **Cor. 3.3 does not by itself
give even-$|J^*|$.** What it gives is a *proved instance of the box mechanism* on the clean
(unit-evaluated) polygon.

**The unifying statement.** For a polynomial $F(z)=\sum c_k z^k$ evaluated at a unit, the
min-valuation locus is a box $j_0+\{\text{subset sums of distinct }2^{a_s}\}$ **iff**
$F(z)/2^{e}\equiv z^{j_0}(1+z)^{g}\pmod2$ (with $e=\min_kv_2(c_k)$), because over $\mathbb F_2$,
$(1+z)^{g}=\prod_s(1+z^{2^{a_s}})$ exactly when $g=\sum_s2^{a_s}$ (Lucas), whose support is that box.
For the *sharp* polygon one first factors out the parity tilt: with $J^*$ in the even class
($r=2s$), the locus is the bottom-mod-2 support of $\sum_s \binom m{2s}R_{2s}2^s\,y^s$ (the $2^s$ is
the tilt), and similarly for the odd class.

> **Conjecture (Box).** *Coarse:* $\Phi(z)/2^{e}\equiv z^{j_0}(1+z)^{g}\pmod2$.
> *Sharp:* the parity-restricted, tilt-scaled polynomial above is $\equiv y^{l_0}(1+y)^{g'}\pmod2$.
> Hence both $B^*$ and $J^*$ are affine 2-adic boxes, $|B^*|,|J^*|$ are powers of two, and on a tie
> $|J^*|$ is even.

*Verification (exhaustive, $m\le 12$):* coarse $\Phi/2^e\equiv z^{j_0}(1+z)^g$ holds $4114/4114$
shapes; sharp bottom-mod-2 $=y^{l_0}(1+y)^{g'}$ **and** equals the true minimum locus $1624/1624$
ties; $B^*$ and $J^*$ are affine boxes $4114/4114$ and $1624/1624$ respectively;
$|B^*|,|J^*|\in\{1,2,4,8\}$ always.

**Why the engine proves only the leading level.** Cor. 3.2 is the $r=0$ slice of Thm 3.1. The
correction terms $r\ge1$ involve $e_2\bmod 2$, and unlike $p_1^2\equiv p_2$ there is *no* collapse
of the two generators at the $e_2$-level — $e_2=(p_1^2-p_2)/2$ is not congruent to a single power
sum mod 2. So when $\chi^\lambda(2^m)$ is even (exactly the tie case, by Cor. 2.3) the leading slice
$(1+z)^m$ vanishes mod 2 and one must read off the next 2-adic layer, where the clean
"$p_1^2\equiv p_2$" identity is no longer available. This is the same phenomenon as the
**unbounded cancellation depth** of the §6 tower; a uniform proof needs a mod-$\pi^{k}$ (not just
mod-2) statement that the layer is again a power of $(1+z)$.

---

## 5. Status ledger

**Proved (rigorous):** Lemma 1.1 + the $\Psi(\pi)$ reformulation; Prop. 2.1 ($G\equiv\chi^\lambda(2^m)\bmod\pi$);
Lemma 2.2 + Cor. 2.3 (parity criterion in one line; ties have $f^\lambda$ even); Thm 3.1 (exact
$(1+z)$-lift); Cor. 3.2 (the engine); Cor. 3.3 (coarse box when $\chi^\lambda(2^m)$ odd).

**Conjectural (verified $m\le12$, $100\%$):** the Box conjecture (coarse and sharp), hence
even-$|J^*|$ on ties. Reduced to: *the leading 2-adic layer of the (tilt-scaled) generating
polynomial is a shifted power of $(1+\,\cdot\,)$ mod 2.* The engine proves the $r=0$ layer; the
obstruction is the $e_2\bmod2$ correction.

**Not addressed (still open, as flagged in PROVE.md):** Step 2 / non-vanishing. Even with the full
Box conjecture, on a tie the leading $\pi$-coefficient *cancels* ($w\equiv|J^*|\equiv0$); the true
$v_\pi(G)$ is governed by the interleaved even/odd cancellation tower whose depth is unbounded in
$m$ ($d=v_\pi(G)-\mu$ ranges $1..22$ for $m\le12$, and $918/1624$ ties need indices outside $J^*$).
Nothing here closes that. The unique total-vanisher remains $(2,2)$ (now visible as: it is the
unique $\lambda$ with $\Psi(Y)=R_0(1+\cdots)$ giving $(z^2+1)\mid\Phi$; $\Phi_{(2,2)}=2(1+z^2)$).

---

## 6. Verification summary (exact, pure-Python MN; `code/job_*.py`)

| Claim | Result |
|---|---|
| Lemma 1.1 / $\Psi(\pi)=G$ (E1 exact lift Thm 3.1) | $1537/1537$, $m\le10$ |
| $R_0=\chi^\lambda(2^m)$ | $1537/1537$ |
| Prop. 2.1 $G\equiv\chi^\lambda(2^m)\bmod\pi$ | $0$ failures, $m\le11$ |
| Lemma 2.2 $\chi_k\equiv f^\lambda\bmod2$ | $26894/26894$ |
| $\chi^\lambda(2^m)$ odd $\iff f^\lambda$ odd | $0$ exceptions, $m\le12$ |
| Cor. 3.2 engine $\Phi\equiv\chi^\lambda(2^m)(1+z)^m\bmod2$ | $4114/4114$, $m\le12$ |
| Cor. 3.3 $\chi^\lambda(2^m)$ odd $\Rightarrow B^*=\mathrm{submasks}(m)$ | $502/502$ |
| Box (coarse) $\Phi/2^e\equiv z^{j_0}(1+z)^g\bmod2$ | $4114/4114$ |
| Box (sharp) bottom-mod2 $=y^{l_0}(1+y)^{g'}=J^*$ | $1624/1624$ ties |
| $B^*,J^*$ affine 2-adic boxes; sizes $\in\{1,2,4,8\}$ | $100\%$ |
| $J$-locus $=$ $R$-locus (two sharp expansions agree) | $913/913$ ties, $m\le11$ |

Scripts: `code/job_v2D_explore.py`, `job_bpicture.py`, `job_submask.py`, `job_box_detail.py`,
`job_mod2_engine.py`, `job_hierarchy.py`, `job_sharp_mod2.py`, `job_verify_engine.py`
(all on `symfunc.py`/`job1_tie_census.py` machinery).
