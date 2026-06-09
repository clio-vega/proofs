# PROVE 2026-06-09 — |J*| even crux for d=4 fiber vanishing

## Target
`G_λ(i)=⟨s_λ,ψ^m⟩=0 ⟺ λ=(2,2)`, ψ=h₂+ie₂, n=2m. Crux: for λ≠(2,2), |J*| (number of
valuation-minimal terms in the (1+i)-adic expansion) is EVEN ⟹ leading π-coeff cancels;
then need next-order survival. Strong-partial goal: prove |J*| even structurally.

## Clean reformulations (hand-verified)
- `i^j(1+i)^j = (i-1)^j` (since i(1+i)=i-1). So
  **G_λ(i) = A_λ(i-1)**, A_λ(x)=⟨s_λ,(h₁²+x e₂)^m⟩ = Σ_j C(m,j) M_j x^j ∈ ℤ_{≥0}[x],
  M_j=⟨s_λ,h₁^{2(m-j)} e₂^j⟩, M_0=f^λ.
- h₁²=h₂+e₂ ⟹ h₁²+x e₂ = h₂+(1+x)e₂ ⟹ **A_λ(x)=P_λ(1+x)**, P_λ(t)=⟨s_λ,(h₂+t e₂)^m⟩.
  Since i-1 = i·... with t=1+x=i: **G_λ(i)=P_λ(i), and G=0 ⟺ (t²+1)|P_λ(t)**.
  P_{(2,2)}(t)=t²+1 EXACTLY. So thm ⟺ (t²+1)∤P_λ for λ≠(2,2), P_λ∈ℤ_{≥0}[t].
- ω-duality: ν_k(λ):=⟨s_λ,h₂^{m-k}e₂^k⟩ satisfies ν_k(λ)=ν_{m-k}(λ') ⟹ P_{λ'}(t)=t^m P_λ(1/t).
- Char-sum: G_λ(i) = 2^{-m}(1+i)^m Σ_k C(m,k)(-i)^k χ^λ(2^k 1^{2(m-k)})
  = ((1+i)/2)^m ⟨s_λ,(p₁²-i p₂)^m⟩. So G=0 ⟺ ⟨s_λ,(p₁²-ip₂)^m⟩=0.

## Newton-polygon facts (π=1+i, v_π(N)=2v₂(N) for N∈ℤ)
A_λ(x)=Σ c_j x^j, c_j=C(m,j)M_j. G=Σ c_j (i-1)^j, v_π((i-1)^j)=j.
val(j)=j+2v₂(c_j). μ=min val, J*=argmin.
- val(j+1)-val(j)=1+2(v₂(c_{j+1})-v₂(c_j)) is ODD ⟹ J* in single parity class. (trivial)
- |J*|=1 ⟹ no cancel ⟹ G≠0 (72%, PROVED).
- |J*|≥2: leading π-coeff ≡ |J*| mod π (units ≡1). |J*|=2 ALWAYS cancels (1+1≡0).
  So crux = |J*| is never ODD≥3, i.e. never 3 (and never odd in general).

## REAL reformulation (verified m≤6, 0 mismatches)
Group j by j mod 4 since (i-1)^4 = -4:
**Re(G)=S_0-S_1+2S_3, Im(G)=S_1-2S_2+2S_3**, where
S_r = Σ_s C(m,4s+r) M_{4s+r} (-4)^s ∈ ℤ.
G=0 ⟺ both vanish ⟺ S_0=2(S_2-2S_3), S_1=2(S_2-S_3). [forces S_0,S_1 even]
- S_0 ≡ M_0 = f^λ (mod 2) [recovers parity criterion]; S_1 ≡ m·M_1 (mod 2).
- (2,2): S=[2,2,1,0] ⟹ Re=Im=0. ✓
NOTE: earlier "G=0 ⟺ E=0 and O=0 separately" is FALSE — v_π(E) can be pushed to odd
valuation by internal cancellation (e.g. (2,2): E=2-2i=-π³, v_π=3 odd), so E and (i-1)O
DO cancel. Only Re=Im=0 (2 conditions) is correct.

## ν-basis form (cleanest balance conditions)
P_λ(t)=Σ_k C(m,k)ν_k t^k, ν_k=⟨s_λ,h₂^{m-k}e₂^k⟩≥0, ν_k(λ)=ν_{m-k}(λ') [ω-duality].
P_λ(t)>0 ∀t>0 (nonneg coeffs, P_λ(1)=f^λ>0). G=0 ⟺ ±i are roots ⟺
  Re: Σ_l (-1)^l C(m,2l) ν_{2l} = 0   (even balance)
  Im: Σ_l (-1)^l C(m,2l+1) ν_{2l+1} = 0  (odd balance)
(2,2): ν=[1,0,1] ⟹ Re=1-1=0, Im=0. ✓

## Goal restatement
Within one parity class of j, the minimal-valuation locus of val(j)=j+2v₂(C(m,j)M_j)
has EVEN cardinality (for λ≠(2,2)). Equivalently the v₂(M_j) subsequence in a parity class.

## KEY LOGICAL CORRECTION
ODD |J*| is the GOOD case: w≡|J*|≡1 mod π ≠0 ⟹ v_π(G)=μ<∞ ⟹ G≠0. The 72% (|J*|=1) ARE this.
EVEN |J*| is the HARD case (leading cancels). So "prove evenness" ≠ "prove non-vanishing":
evenness just confirms the easy argument FAILS. The full thm needs the next-order (Step 2)
regardless. Proving evenness alone is "strong partial" only.

## Exact leading π-coefficient (VERIFIED m≤7, 0 parity mismatches)
w = Σ_{j∈J*} o_j · i^{(3μ-j)/2}, o_j = odd part of c_j=C(m,j)M_j.
w ≡ |J*| mod π. For a PAIR J*={a,a+2}: w ∝ o_{a+2}+i·o_a, v_π(w)=v₂(o²+o²)=1 (odd²+odd²≡2 mod8).
For J*={a,a+4}: w ∝ o_a−o_{a+4} (real, even), v_π(w)=v₂(o_a−o_{a+4})≥2. So pair-depth depends
on index-gap mod 4.

## THE OBSTRUCTION (why depth is unbounded — precise)
v_π(G) is NOT μ+v_π(w_leading). The even-class part G_E=Σ_{even j}c_j(i-1)^j and odd-class
G_O interleave: each is itself a sum-of-units that can cancel to arbitrary depth. Example
(3,3,1,1) m=4: even pair {0,4} cancels to π-depth 10, but ODD pair {1,3} (at val μ+1=7)
cancels to depth 8 < 10, so v_π(G)=8 governed by the ODD class. The two towers compete;
the actual v_π is the min over a recursive interleaving. Self-similar-ish but the polynomial
changes each level (G_E=Ã(−2i), Ã=Σ c_{2l}y^l, different divisibility y²+4) ⟹ no clean recursion.
This is why no fixed π-adic level closes it and why cancellation depth (5,9,11,15,15,21,...) is unbounded.

## VERIFIED FACTS (this session)
- (2,2) is the UNIQUE vanisher for all m≤10 (n≤20): 627 shapes at m=10, 0 other vanishers.
- |J*| ∈ {1,2,4} ALWAYS (never odd>1, never 8) through m=10. Powers of 2 ≤4.
- Integer-point criterion: q(n)|A_λ(n) ∀n ⟺ G=0. Min witness n grows (n=5 at m=8,9) ⟹
  NO uniform finite witness set. Criterion is EQUIVALENT to vanishing (no free lunch):
  A_λ(n)≡R(n) mod q(n), R(x)=r₁x+r₀ remainder, r₁=Im G, r₀=Re G+Im G.

## ASSESSMENT
Full theorem & evenness mechanism remain OPEN (consistent with Robin's program: last gap,
many routes pruned). Genuine new contributions: clean reformulations (q|A_λ, A_{(2,2)}=q
exact), real reformulation (4 integer sums S_r), exact leading-coeff formula, precise
obstruction description, extended verification m≤10. Deliver as rigorous partial + precise gap.
