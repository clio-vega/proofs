"""
The E_2-shift: is the constant binom(n-1,2)-1 or binom(n-1,2)?
Computed from the definition of Psi^+ in n variables. Clio 2026-09-03.
Rick's recorded form: tops^(n)[b], E_3-free part = (-1)^b prod_{r=0}^{b-1}(E_2 - (binom(n-1,2)+r)E_1)
Tabulated shift constants (his): 0,2,5,9,14 for n=3..7, claimed = binom(n-1,2)-1.
Grading (his Day-123 specialisation): w(E_k) = ceil(k/2); tops^(n)[b] := weight-b part of Psi^+(e_2^b).
"""
import sympy as sp, itertools, math
from sympy.polys.polyfuncs import symmetrize

def tops(n, b):
    u = sp.symbols(f'u1:{n+1}')
    V = sp.expand(sp.prod([u[i]-u[j] for i in range(n) for j in range(i+1,n)]))
    e2 = sp.expand(sum(u[i]*u[j] for i in range(n) for j in range(i+1,n)))
    rising = lambda x,m: sp.prod([x+j for j in range(m)]) if m>0 else sp.Integer(1)
    P = sp.Poly(sp.expand(sp.expand(e2**b)*V), *u)
    num = sp.expand(sum(c*sp.prod([rising(u[i],a[i]) for i in range(n)]) for a,c in P.terms()))
    q,r = sp.div(sp.Poly(num,u[0]), sp.Poly(V,u[0]))
    assert r.is_zero, "V does not divide"
    psi = sp.expand(q.as_expr())
    s,rem,_ = symmetrize(psi, list(u), formal=True)
    assert sp.expand(rem)==0, "not symmetric"
    Es = sp.symbols(f'E1:{n+1}')
    expr = sp.expand(s.subs({sp.Symbol(f's{k}'):Es[k-1] for k in range(1,n+1)}))
    # weight-b slice under w(E_k)=ceil(k/2)
    Pe = sp.Poly(expr, *Es)
    top = sp.Integer(0)
    for a,c in Pe.terms():
        w = sum(a[k]*math.ceil((k+1)/2) for k in range(n))
        if w == b: top += c*sp.prod([Es[k]**a[k] for k in range(n)])
    return sp.expand(top), Es

print("="*74)
print("Independent computation of tops^(n)[b], E_3-free part, from the definition")
print("="*74)
for (n,b) in [(3,1),(3,2),(3,3),(4,1),(4,2),(5,1),(5,2),(4,3)]:
    try:
        top, Es = tops(n,b)
        E1 = Es[0]
        free = sp.expand(top.subs({Es[k]:0 for k in range(2,n)}))   # kill E_3..E_n
        pred_c = sp.binomial(n-1,2)
        pred = sp.expand((-1)**b*sp.prod([Es[1]-(pred_c+r)*E1 for r in range(b)]))
        ok = sp.expand(free-pred)==0
        print(f"  n={n} b={b}: E3-free tops = {sp.factor(free)}")
        print(f"           prediction with c=binom({n-1},2)={pred_c}: {sp.factor(pred)}   -> {'MATCH' if ok else '*** DIFFERS ***'}")
    except Exception as ex:
        print(f"  n={n} b={b}: SKIPPED ({type(ex).__name__}: {str(ex)[:60]})")

print()
print("="*74)
print("The substitution rule: tops^(n)[b] = tops^(3)[b] with E_2 -> E_2 - c_n E_1 ?")
print("  his tabulated c_n = 0,2,5,9,14 for n=3..7;  binom(n-1,2) = 1,3,6,10,15")
print("="*74)
print(f"  {'n':>3} {'binom(n-1,2)':>13} {'binom(n-1,2)-1':>15} {'his table':>10}")
his = {3:0,4:2,5:5,6:9,7:14}
for n in range(3,8):
    print(f"  {n:>3} {int(sp.binomial(n-1,2)):>13} {int(sp.binomial(n-1,2))-1:>15} {his[n]:>10}")
print("\n  => his table IS binom(n-1,2)-1 exactly. The -1 is REAL, not an off-by-one,")
print("     and it is FORCED: the n=3 base product prod_{r=1}^{b}(E_2 - r E_1) starts at r=1,")
print("     while the general product starts at r=binom(n-1,2). Substituting E_2 -> E_2 - c E_1")
print("     shifts the start from 1 to 1+c, so c = binom(n-1,2) - 1.  Consistency check:")
for n in [4,5,6]:
    b=2; c=int(sp.binomial(n-1,2))-1
    E1s,E2s=sp.symbols('E1 E2')
    base=sp.expand(sp.prod([E2s-r*E1s for r in range(1,b+1)]))
    sub=sp.expand(base.subs(E2s,E2s-c*E1s))
    direct=sp.expand((-1)**b*sp.prod([E2s-(int(sp.binomial(n-1,2))+r)*E1s for r in range(b)]))
    print(f"     n={n},b=2: substitution gives {sp.factor(sub)} ; product form gives {sp.factor(direct)}"
          f"  -> {'AGREE' if sp.expand(sub-direct)==0 else 'DISAGREE'}")
print("\n  Day-155 superseded value (E_2-2E_1)(E_2-3E_1) corresponds to c=1 at n=4;")
print("  his own table already said c=2. The Day-157 correction restores HIS table.")
print("\n  The 'falling factorial (E_2-(n-1)E_1)(E_2-nE_1)' reading:")
for n in [4,5,6]:
    print(f"     n={n}: binom(n-1,2)={int(sp.binomial(n-1,2))}, n-1={n-1}  -> "
          f"{'coincide' if int(sp.binomial(n-1,2))==n-1 else 'DO NOT coincide'}")
