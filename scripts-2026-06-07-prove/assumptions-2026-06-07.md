## What I'm trying to show
Two-row d=4 fiber law: G_{(2m-b,b)}(i) = 0  ⟺  (m,b)=(2,2), for integers 0≤b≤m.
Reduced (prior cycles) to: I_b(m) := Im G_{(2m-b,b)}(i) ≠ 0 for all integer m≥b, except (2,2).
Sharper (♦): Q_b (deflation of I_b by forced roots {0..⌊(b-1)/2⌋}) has no rational root, save (2b-1)/2 when 4|b.

## What's going wrong
Browse said "Q_b is binary-Krawtchouk-type; import Cullinan/Filaseta irreducibility." But computational
sweep TODAY refutes the identification: I_b/Q_b is NOT Meixner/Hahn/dual-Hahn/Krawtchouk/Chebyshev/MO#286705
up to scalar+affine. Discriminants non-square (Galois ∉ A_n). So no known irreducibility theorem to cite.
2-adic route dead (v₂ unbounded, flat Newton polygon). No reflection symmetry of I_b for b≥4. Stuck on a
uniform-in-b "no integer/rational root" with no certificate and no named family.

## Every Assumption I Am Making

### About the objects
A1. I_b(m) is correctly = Im[u^b]((1-u)(1+su+u²)^m), s=1+i. (matches writeup, verified m≤13)
A2. The two-row law is EQUIVALENT to "I_b(m)≠0 for integer m≥b except (2,2)" — i.e. Im alone decides it.
A3. The vanishing locus is exactly {(2,2)}; no other integer (m,b), m≥b, kills G.
A4. "m≥b" is the right range (partition condition 2m-b≥b). Edge: m=b is shape (b,b) — included.
A5. (2,2) corresponds to m=2,b=2; it is the UNIQUE vanisher among two-row shapes.

### About the maps / the reduction
A6. Reducing the FULL law (Re=Im=0) to "Im≠0" is valid: if Im G≠0 then G≠0. (one-directional, safe)
A7. ...but the CONVERSE leakage: could there be integer (m,b) with Im G=0 yet Re G≠0 (so G≠0, law still holds)?
    If so, "I_b(m)≠0 for m≥b" is SUFFICIENT but maybe NOT NECESSARY. I've been assuming I must prove Im≠0
    for ALL m≥b. But the LAW only needs G≠0, i.e. NOT(Re=0 AND Im=0). Im=0 alone does not violate the law!
A8. The forced roots {0..⌊(b-1)/2⌋} are all < b, so irrelevant to range m≥b. (verified)

### About definitions / classical families
A9. "Binary-Krawtchouk-type" identification — ASSUMED by browse. REFUTED today.
A10. Q_b irreducible for 4∤b — computational only (b≤16). Proof is the open hard part.

### About prior results / tools
A11. Parity criterion G_λ≡f^λ mod(1+i): f^λ odd ⟹ G≠0. So open shapes are f^λ-EVEN ones only.
A12. Structural identity G_λ(i)=Σ_k C(m,k)i^k N_{λ,k}, N_{λ,k}=⟨s_λ,h₂^{m-k}e₂^k⟩≥0.
A13. Palindromicity f(u)=u²f(1/u) ⟹ for INTEGER m, C_l(m)=C_{2m-l}(m). (true, integer m)

### About the problem itself
A14. ASSUMING the right object is Im G = I_b. But maybe |G|² or Re G has a cleaner non-vanishing proof.
A15. ASSUMING uniform-in-b is the only acceptable target. Maybe a clean infinite family + honest gap is the
     real deliverable, and "uniform (♦)" is genuinely a hard open number-theory problem (not closable now).
A16. ASSUMING the two-row law is even the right battle — PROVE.md's PRIMARY target was Gap A (all λ), and the
     fast-path was conditional on a usable OP result that does NOT exist.
A17. ASSUMING I_b is "generic." But it has rich structure (forced roots, 4|b factor, palindromic seed) — maybe
     a NON-OP but still exploitable structure (e.g., a recurrence in b, a contiguous relation) decides roots.

### The ones that feel too obvious
A18. That I must rule out roots of Q_b at all. The LAW (A7!) only needs the JOINT vanishing to be exactly (2,2).
     I have been silently upgrading "G≠0" to "Im G≠0", which is STRICTLY STRONGER and may be FALSE/unprovable.
A19. That Im G=0 at integer m≥b actually OCCURS for some shape other than (2,2). Does it? If Im G=0 never
     happens for m≥b except (2,2), then A2 holds and the strong target is the right one. If Im G=0 DOES happen
     at some integer m≥b (with Re≠0), then the strong target is FALSE and I've been chasing an unprovable claim.
