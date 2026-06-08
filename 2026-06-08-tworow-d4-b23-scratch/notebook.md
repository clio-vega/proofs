# Prove session 2026-06-08 — close b ≡ 0 (mod 4)

## Target
Identity (V): for b ≡ 0 (mod 4), R=(b-2)/2, every integer m≥b:
   v₂(I_b(m)) = v₂((m)_{R+1}) − v₂(R!).
⟹ Q_b(m) ≠ 0 ⟹ two-row d=4 law for b ≡ 0 (mod 4).

## Key structural facts (b ≡ 0 mod 4 ⟹ b=4t)
- R = (b−2)/2 = 2t−1 is **ODD**; R+1 = b/2 even.
- (E): I_b = Σ_{j=1}^b τ_j, τ_j = C(m,j)·Im((1+i)^j)·(T(b−j)−T(b−1−j)).
- Im((1+i)^j): v₂ = ⌊j/2⌋ for j≢0 mod4; =0 (term vanishes) for j≡0 mod4. Odd iff j=1.
- τ_1 = −(R+1)C(m,R+1), v₂(τ_1) = a := v₂((m)_{R+1}) − v₂(R!).

## The decomposition (rederived for b even; CONFIRMED 24480 cases)
For τ_j ≠ 0, h := ⌊j/2⌋, s_j = R+1+h, β_j = s_j − j:
   v₂(τ_j) − a  =  D0(j) + S_j(m),    S_j = Σ_{t=1}^{h} v₂(m−R−t)
with the m-independent offset
   j odd  : D0 = v₂(C(R, h))                    ≥ 0
   j even : D0 = v₂(C(R, h−1)) − v₂(h)          (can be < 0)

## (P1) no dip:  D0 + S_j ≥ 0 always
- j odd : D0 ≥0, S_j ≥0. done.
- j even: product (m−R−1)…(m−R−h) is h consecutive ints ⟹ divisible by h! ⟹
          S_j ≥ v₂(h!) ≥ v₂(h) ≥ −D0.  done.

## Which j≥2 can tie (D0+S_j = 0)?
- j even, h≥3: S_j ≥ v₂(h!) = v₂(h)+v₂((h−1)!) ≥ v₂(h)+1 (since (h−1)!≥2 even) > −D0 ⟹ strict.
- j even, h=2 ⟹ j=4 ≡0 mod4 ⟹ τ=0 (excluded).
- j even, h=1 ⟹ j=2: D0=0, S_2=v₂(m−R−1). TIE ⟺ m−R−1 odd ⟺ m odd (R odd).
- j odd, h≥2 (j≥5): S_j has ≥2 consecutive ints ⟹ ≥1; D0≥0 ⟹ strict.
- j odd, h=1 ⟹ j=3: D0=v₂(C(R,1))=v₂(R)=0 (R odd), S_3=v₂(m−R−1). TIE ⟺ m odd.
⟹ ONLY j∈{2,3} can tie τ_1, and they BOTH tie ⟺ m odd. (matches data: even m→{1}, odd m→{1,2,3})

## Punchline (V)
Let a = v₂(τ_1). Every term has v₂ ≥ a (P1). I_b/2^a = Σ_j τ_j/2^a integer; mod 2 only ties survive.
- m even: only τ_1 ties ⟹ I_b/2^a ≡ τ_1/2^a ≡ 1 (odd). 
- m odd : ties = {1,2,3}, each τ_j/2^a odd ⟹ sum of THREE odds = odd ⟹ I_b/2^a odd.
Either way v₂(I_b)=a finite ⟹ I_b≠0 for m≥b. ∎

Beautiful contrast with b≡1: there τ_1 strictly dominates (one odd term). Here for odd m three
terms tie, but |ties| is odd (1 or 3) so the leading 2-adic digit always survives. The parity of
the NUMBER of minimal terms is the real invariant.

## Verified
verify.py: 24480 checks, b∈{4,...,64} step4, m∈[b,b+20]∪40 random up to 10^6. 0 problems.
explore.py: per-term valuation tables b=4,8,12,16.
