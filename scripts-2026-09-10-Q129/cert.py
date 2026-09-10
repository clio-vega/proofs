import sympy as sp
from ribbon import *
t=sp.Symbol('t')
def preds(lam,N):
    out=[]
    for L in range(1,sum(lam)+1):
        for g,ht in remove_ribbons(lam,L,N+3): out.append((g,L,ht))
    return out

n=4; mu=(4,)
P=list(partitions(n))
pm={g:(L,ht) for g,L,ht in preds(mu,n)}
gams=sorted(pm)
print("mu =",mu," rim-hook predecessors (gamma, L, ht(mu/gamma)):",[(g,)+pm[g] for g in gams])
w={g:sp.Symbol('w'+str(sum(g))) for g in gams}
rows=[];rhs=[];labels=[]
for lam in P:
    pl={g:(L,ht) for g,L,ht in preds(lam,n)}
    coeffs=[]
    for g in gams:
        coeffs.append(t**pl[g][1] if (g in pl and pl[g][0]==pm[g][0]) else 0)
    rows.append(coeffs); rhs.append(1 if lam==mu else 0); labels.append(lam)
M=sp.Matrix(rows); b=sp.Matrix(rhs)
print()
print("      unknowns:", [str(w[g]) for g in gams], "   (w_k = wt_B((4), (k)))")
for i,lam in enumerate(labels):
    print(f"  lam={str(lam):14s} " + "  ".join(f"{sp.sstr(M[i,j]):>6s}" for j in range(M.cols)) + f"   =  {b[i]}")
print()
ns=M.T.nullspace()
print("left null space of M (certificates), dim =",len(ns))
for v in ns:
    v=sp.simplify(v*sp.lcm([sp.denom(sp.cancel(x)) for x in v]))
    pair=sp.factor(sp.expand((v.T*b)[0,0]))
    print("   c =",[sp.factor(x) for x in v], " ->  c.b =",pair)

print()
print("verify c^T M = 0 :", sp.simplify(sp.Matrix([[t**2*(t+1),-t**2*(t+1),-t*(3*t-1),-(t-1)**2,2*(t-1)]])*M)==sp.zeros(1,4))
print("c^T b =", sp.factor(sp.expand(sp.Matrix([[t**2*(t+1),-t**2*(t+1),-t*(3*t-1),-(t-1)**2,2*(t-1)]])*b)[0,0]))
print()
print("=== which specialisations t=c are consistent? (rank test on the specialised system) ===")
for val in [-1,0,1,2,-2,sp.Rational(1,2),sp.I,sp.Symbol('gen')]:
    if val==sp.Symbol('gen'):
        Ms,bs=M,b; nm='generic'
    else:
        Ms,bs=M.subs(t,val),b.subs(t,val); nm=f"t={val}"
    print(f"  {nm:10s}: rank(M)={Ms.rank()}  rank([M|b])={Ms.row_join(bs).rank()}  -> {'SOLVABLE' if Ms.rank()==Ms.row_join(bs).rank() else 'INCONSISTENT'}")
