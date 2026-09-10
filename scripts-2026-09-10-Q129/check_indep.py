"""Independent verification, deliberately NOT reusing the abacus.
(a) direct border-strip enumeration on Young diagrams
(b) numeric Schur polynomials via bialternant -> R_e(-1) == mult by p_e
"""
import sympy as sp, random
from ribbon import partitions, all_parts_upto, add_ribbons, conj

# ---------- (a) direct border-strip enumeration ----------
def cells(lam): return {(i,j) for i,p in enumerate(lam) for j in range(p)}
def is_partition(v):
    v=[x for x in v if x>0]
    return all(v[i]>=v[i+1] for i in range(len(v)-1))

def add_ribbons_direct(lam, e, maxn=99):
    """brute force: all mu >= lam with |mu/lam|=e, mu/lam connected & no 2x2."""
    out=[]
    R=len(lam)+e+1
    base=list(lam)+[0]*(R-len(lam))
    def rec(i, cur, left):
        if i==R:
            if left==0: yield tuple(x for x in cur if x>0)
            return
        lo=base[i]; hi=cur[i-1] if i>0 else lo+left
        for v in range(lo, min(lo+left, hi)+1):
            yield from rec(i+1, cur+[v], left-(v-lo))
    for mu in rec(0,[],e):
        S=cells(mu)-cells(lam)
        if len(S)!=e: continue
        # connected (edge adjacency)
        seen={next(iter(S))}; st=[next(iter(S))]
        while st:
            (a,b)=st.pop()
            for d in ((1,0),(-1,0),(0,1),(0,-1)):
                c=(a+d[0],b+d[1])
                if c in S and c not in seen: seen.add(c); st.append(c)
        if seen!=S: continue
        # no 2x2
        if any((a,b) in S and (a+1,b) in S and (a,b+1) in S and (a+1,b+1) in S for (a,b) in S): continue
        rows={a for (a,b) in S}
        out.append((mu, len(rows)-1))
    return sorted(out)

# ---------- (b) numeric symmetric functions ----------
def schur_num(lam, xs):
    m=len(xs); lam=list(lam)+[0]*(m-len(lam))
    num=sp.Matrix(m,m, lambda i,j: sp.Rational(xs[i])**(lam[j]+m-1-j))
    den=sp.Matrix(m,m, lambda i,j: sp.Rational(xs[i])**(m-1-j))
    return sp.simplify(num.det()/den.det())
def p_num(k, xs): return sum(sp.Rational(x)**k for x in xs)

print("=== (a) abacus height vs direct Young-diagram height ===")
bad=0; tot=0
for e in (2,3,4,5):
    for n in range(0,7):
        for lam in partitions(n):
            A=sorted(add_ribbons(lam,e,len(lam)+e+n+3))
            B=add_ribbons_direct(lam,e)
            tot+=1
            if A!=B:
                bad+=1; print("MISMATCH",e,lam,A,B)
print(f"  {tot-bad}/{tot} shapes agree exactly (shape sets AND heights)")

print("=== (b) R_e(-1) s_mu  ==  p_e * s_mu   (numeric, bialternant) ===")
random.seed(11)
bad=0; tot=0
for e in (1,2,3,4):
    for n in range(0,6):
        for mu in partitions(n):
            m=n+e+1
            xs=[sp.Rational(random.randint(2,40), random.randint(1,7)) for _ in range(m)]
            while len(set(xs))<m: xs=[sp.Rational(random.randint(2,60),random.randint(1,7)) for _ in range(m)]
            lhs=sum((-1)**ht*schur_num(lam,xs) for lam,ht in add_ribbons(mu,e,len(mu)+e+n+3))
            rhs=p_num(e,xs)*schur_num(mu,xs)
            tot+=1
            if sp.simplify(lhs-rhs)!=0: bad+=1; print("  FAIL",e,mu)
print(f"  {tot-bad}/{tot} (e,mu) pairs verified numerically")
