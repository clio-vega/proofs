"""Is the M-convexity failure STRUCTURAL for two-variable-per-row homogeneous weights?

Row-degree lemma (proof, not computation): each row has exactly m vertices, each of
weight homogeneous of degree 1 in (x_i,y_i).  So EVERY monomial of Z has
    deg_{x_i} + deg_{y_i} = m   for every i.
Hence supp(Z) lies on the graph G = {(alpha, m*1 - alpha)} inside Z^{2n}.

Test: over MANY weight families (free-fermion and not), is every non-monomial Z
non-M-convex?  And is the x-support always confined to one total-degree level?
"""
import sympy as sp, itertools, random
from fractions import Fraction
from sixvertex import Z
from lorentzian import is_lorentzian, is_M_convex, multi_factorial

def cdict(expr, vars_):
    p=sp.Poly(sp.expand(expr),*vars_)
    return {tuple(mm):Fraction(int(c)) for mm,c in zip(p.monoms(),p.coeffs())}

random.seed(11)
tot=mono=nontriv=notMconv=xlevel1=0
ffcount=0
for trial in range(60):
    n,m = random.choice([(2,2),(2,3),(3,2),(2,4)])
    xs=sp.symbols('x1:%d'%(n+1)); ys=sp.symbols('y1:%d'%(n+1))
    # random degree-1 nonneg weights, independent per vertex type (shared across rows)
    def lf(xi,yi):
        a,b=random.randint(0,3),random.randint(0,3)
        if a==0 and b==0: a=1
        return a*xi+b*yi
    pat=[(random.randint(0,3),random.randint(0,3)) for _ in range(6)]
    pat=[(1,0) if p==(0,0) else p for p in pat]
    def mk(i):
        xi,yi=xs[i],ys[i]
        nm=['a1','a2','b1','b2','c1','c2']
        return {k:(p[0]*xi+p[1]*yi) for k,p in zip(nm,pat)}
    w=[mk(i) for i in range(n)]
    isff = sp.expand(w[0]['a1']*w[0]['a2']+w[0]['b1']*w[0]['b2']-w[0]['c1']*w[0]['c2'])==0
    if isff: ffcount+=1
    vars_=list(xs)+list(ys)
    for top in itertools.product([0,1],repeat=m):
        for bot in itertools.product([0,1],repeat=m):
            z=Z(n,m,list(top),list(bot),[0]*n,[0]*n,w)
            if z==0: continue
            cd=cdict(z,vars_); tot+=1
            # verify the row-degree lemma
            for a in cd:
                for i in range(n):
                    assert a[i]+a[n+i]==m, ("ROW-DEGREE LEMMA VIOLATED",a,i,m)
            if len(cd)==1: mono+=1; continue
            nontriv+=1
            if not is_M_convex(set(cd)): notMconv+=1
            xs_supp={tuple(a[:n]) for a in cd}
            if len({sum(t) for t in xs_supp})==1: xlevel1+=1
print("random weight families tested (free-fermion among them: %d)"%ffcount)
print("row-degree lemma deg_{x_i}+deg_{y_i}=m : VERIFIED on every monomial of every Z")
print("total nonzero Z      : %d"%tot)
print("  monomial           : %d"%mono)
print("  NON-MONOMIAL       : %d"%nontriv)
print("    not M-convex     : %d   <-- %s"%(notMconv, "ALL of them" if notMconv==nontriv else "NOT all!"))
print("    x-support in a single total-degree level : %d / %d"%(xlevel1,nontriv))
