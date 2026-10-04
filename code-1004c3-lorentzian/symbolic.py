"""Independent symbolic verification of Lemma rawhess and Lemma minor.
These two are MY OWN derivation, so they get an instrument that does not share
the derivation: sympy builds N(f), differentiates, and forms the Hessian from
scratch, with symbolic coefficients."""
import sympy as sp, itertools, random
from math import factorial, prod

print("=== Lemma rawhess, symbolically, with SYMBOLIC coefficients ===")
tot=ok=0
for n in (2,3,4):
    for d in (2,3,4,5):
        if n**d > 700: continue
        v=sp.symbols(f'v0:{n}')
        mons=[a for a in itertools.product(range(d+1),repeat=n) if sum(a)==d]
        c={a: sp.Symbol('c_'+'_'.join(map(str,a))) for a in mons}
        Nf=sum(c[a]/prod(factorial(t) for t in a)*prod(v[i]**a[i] for i in range(n))
               for a in mons)
        for beta in [b for b in itertools.product(range(d+1),repeat=n) if sum(b)==d-2]:
            g=Nf
            for i in range(n):
                for _ in range(beta[i]): g=sp.diff(g,v[i])
            Hs=sp.hessian(sp.expand(g), v)
            pred=sp.Matrix(n,n, lambda i,j: c.get(tuple(
                beta[k]+(1 if k==i else 0)+(1 if k==j else 0) for k in range(n)), sp.Integer(0)))
            tot+=1
            if sp.simplify(Hs-pred)==sp.zeros(n,n): ok+=1
            else: print("   MISMATCH",n,d,beta)
print(f"  (n,d,beta) triples enumerated {tot}, agreements {ok}  "
      f"(n up to 4, degree up to 5, coefficients symbolic)")

print("=== Lemma minor: at most one positive eigenvalue + nonneg diagonal => 2x2 minors <=0 ===")
random.seed(7)
tot=ok=viol=0
for n in (2,3,4,5):
    for _ in range(400):
        A=[[random.randint(-6,6) for _ in range(n)] for _ in range(n)]
        M=[[A[i][j]+A[j][i] for j in range(n)] for i in range(n)]
        if any(M[i][i]<0 for i in range(n)): continue
        from lor import count_pos_eigen_exact
        if count_pos_eigen_exact(M)>1: continue
        tot+=1
        bad=any(M[i][i]*M[j][j]-M[i][j]**2>0 for i in range(n) for j in range(n))
        ok+= (not bad); viol+= bad
print(f"  matrices satisfying the hypotheses enumerated {tot}, conclusion held {ok}, "
      f"violations {viol}  (n=2..5, entries in [-12,12])")

print("=== final (Q) gate: the theorem's conclusion, independent enumerator, wider range ===")
from lor import k_sequence, is_logconcave_pf2
tot=bad=0; Ds=[]; sizes=[]
for m in (1,2,3,4,5):
    hmax={1:9,2:7,3:5,4:4,5:3}[m]
    for h in itertools.product(range(1,hmax+1),repeat=m):
        for l in itertools.product(*[range(0,hi+1) for hi in h]):
            for D in range(sum(l),sum(h)+1):
                ks=k_sequence(list(l),list(h),D)
                if not any(ks): continue
                tot+=1; Ds.append(D); sizes.append(sum(ks))
                if not is_logconcave_pf2(ks): bad+=1; print("   FAILURE",l,h,D,ks)
print(f"  instances enumerated {tot}, PF2 failures {bad}; D up to {max(Ds)}, "
      f"largest total count {max(sizes)}, m up to 5")
