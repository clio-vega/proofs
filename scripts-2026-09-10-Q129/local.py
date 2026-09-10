"""K-L local identity for my t-deformed A_n:  S(lam,L) = L-ribbon removals, wt_A = t^ht.
Take the natural symmetric T(mu,L)=S(mu,L); leave wt_B(mu,gam) FREE over Q(t).
Local identity:  sum_{gam in G(lam,mu)} t^{ht(lam/gam)} w(mu,gam) = delta(lam,mu)."""
import sympy as sp
from ribbon import *
t=sp.Symbol('t')

def preds(lam, N):
    """all (gamma, L, ht) with lam/gamma a rim hook of size L>=1"""
    out=[]
    for L in range(1,sum(lam)+1):
        for g,ht in remove_ribbons(lam,L,N+3):
            out.append((g,L,ht))
    return out

for n in (2,3,4,5,6):
    P=list(partitions(n))
    print(f"--- n={n} ---")
    allok_gen=True; allok_m1=True
    for mu in P:
        pm={g:(L,ht) for g,L,ht in preds(mu,n)}
        gams=sorted(pm)
        w={g:sp.Symbol(f'w_{i}') for i,g in enumerate(gams)}
        eqs=[]
        for lam in P:
            pl={g:(L,ht) for g,L,ht in preds(lam,n)}
            s=0
            for g in gams:
                if g in pl and pl[g][0]==pm[g][0]:
                    s+= t**pl[g][1]*w[g]
            eqs.append(sp.expand(s-(1 if lam==mu else 0)))
        unk=list(w.values())
        M,b=sp.linear_eq_to_matrix(eqs,unk)
        rM=M.rank(); rA=M.row_join(b).rank()
        ok_gen = (rA==rM)
        M1=M.subs(t,-1); b1=b.subs(t,-1)
        ok_m1 = (M1.row_join(b1).rank()==M1.rank())
        allok_gen &= ok_gen; allok_m1 &= ok_m1
        if not ok_gen or n<=3:
            print(f"   mu={mu}: #unk={len(unk)} #eq={len(eqs)}  general t: {'SOLVABLE' if ok_gen else 'INCONSISTENT'} | t=-1: {'SOLVABLE' if ok_m1 else 'INCONSISTENT'}")
    print(f"  => n={n}: local identity solvable for general t? {allok_gen} ;  at t=-1? {allok_m1}")
