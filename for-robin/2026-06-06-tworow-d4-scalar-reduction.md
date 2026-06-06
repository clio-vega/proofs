# Two-row d=4 fiber law: a scalar reduction (2026-06-06)

I reduced the two-row $d{=}4$ fiber-vanishing conjecture
($G_{(a,b)}(i)=0 \iff (a,b)=(2,2)$, where $G_{(a,b)}(i)=\langle s_{(a,b)},\psi^m\rangle$,
$\psi=h_2+i\,e_2$) to a genuinely one-dimensional scalar positivity statement. Via the
recessive continued fraction $\xi_{n-1}=A_n/((1+i)n+\xi_n)$, $A_n=(m+n)(m-n+1)$, I pass
to coordinates $p_n=\Re\xi_n-(m-n)$, $q_n=\Im\xi_n+n$, in which the law reads simply
$q_n>0$ (with $q_n=nP_n$). The $q$-recursion couples to $p$ only through the denominator
$D_n=q_n^2+(m+p_n)^2$, so under $p_n\ge0$ it is dominated by a $p$-free scalar comparison
sequence $Q_m=m$, $Q_{n-1}=(n-1)-Q_nA_n/(Q_n^2+m^2)$. The whole law then follows from two
clean lemmas: domination $q_n\ge Q_n$ (F3) and scalar positivity $Q_n>0$ (F4). The sharp
theorem is $q_n\ge0$ for all $m\ge3$ with equality iff $(m,n)=(3,1)$ — exactly the $(2,2)$
fiber. This is the corrected biconditional (the strict ">0" fails at one point). Verified
exhaustively $m=3..400$ plus $m=419$: the only zero of $P_n$ is at $(3,1)$, and the binding
margin $Q_1\to1/4$. Honest residual: neither F3 nor F4 has a closed form yet — the slow
manifold is a Cardano cubic $(\varphi-t)(\varphi^2+1)=-\varphi(1-t^2)$ with no elementary
envelope, the map $g(q)$ is non-monotone, and one-sided bounds fail ($\max_q g<0$ at $n=2$).
The entire margin concentrates in the single step $n=2\to1$, needing a two-sided envelope
$0.6\lesssim Q_2\lesssim0.75$ of accuracy $O(1/m)$ around the cubic root.

GitHub: https://github.com/clio-vega/proofs/blob/2026-06-06-tworow-d4-scalar-reduction/2026-06-06-tworow-d4-scalar-reduction.tex
