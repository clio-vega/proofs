"""Q430 PART 8-fast: the two-tableau fibre, tableau side only (no Warnaar engine)."""
import sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from korff import paths, J_cyclic, mvec
t=sp.Symbol('t')
def sign_hkko(la,k,w):
    h=2*k
    if all(la[2*i]==la[2*i+1] for i in range(k)): return 1
    if la[0]-la[h-1]==w and all(la[2*i+1]==la[2*i+2] for i in range(k-1)): return -1
    return 0
def readings(P,w):
    out=[]
    for a in range(len(P)-1):
        th=tuple(P[a+1][i]-P[a][i] for i in range(len(P[a])))
        J=J_cyclic(th); assert len(J)==1
        i=J[0]; m=mvec(P[a],w)[i]; assert m>=1
        out.append((i,m))
    return out
f=sp.Function('f'); g=sp.Function('g')
print("PART 8 (tableau side) -- alpha=(1,1), ell=0 => w=2, h=2k")
print("  c_alpha(t) = 1 - t for k=1,2,3 [Warnaar engine, PART 7 grid]")
print()
for k in range(1,8):
    h,w=2*k,2
    P=paths(h,w,(1,1))
    tb=sp.Integer(0); ta=sp.Integer(0); lines=[]
    for Pp in P:
        rd=readings(Pp,w); sg=sign_hkko(Pp[-1],k,w)
        lines.append((Pp[-1],sg,tuple(i+1 for i,_ in rd),tuple(m for _,m in rd)))
        tb+=sg*sp.prod([f(m) for _,m in rd]); ta+=sg*sp.prod([g(i+1,m) for i,m in rd])
    tb=sp.expand(tb); ta=sp.expand(ta)
    print(f" k={k} h={h}: #T={len(P)}  in-domain {h>2}")
    for (ep,sg,sl,ms) in lines:
        print(f"      endpoint {str(ep):<22} c^-={sg:+d}  m-INDEX seq {sl}  m-VALUE seq {ms}")
    print(f"      sum c^- = {sum(l[1] for l in lines)}   (= c_alpha(1) = 0? {sum(l[1] for l in lines)==0})")
    print(f"      MULTIPLICITY-DETERMINED sum = {tb}  === 0 ? {tb==0}")
    print(f"      index-aware sum           = {ta}  === 0 ? {ta==0}")
print()
print("named weights at k=2 (h=4,w=2): all multiplicity-determined, target 1-t")
def run(name,om,k=2,w=2):
    tot=sp.Integer(0)
    for Pp in paths(2*k,w,(1,1)):
        tot+=sign_hkko(Pp[-1],k,w)*sp.prod([om(i,m) for i,m in readings(Pp,w)])
    tot=sp.expand(tot); print(f"   {name:<38} sum = {tot}   match 1-t? {tot==sp.expand(1-t)}")
run("Korff 2019 ribbon    f=1",       lambda i,m: sp.Integer(1))
run("Korff 2011 Psi       f=1-t^m",   lambda i,m: 1-t**m)
run("Gaussian binomial    f=[m+1]_t", lambda i,m: sp.expand(sum(t**j for j in range(m+1))))
run("adversarial          f=t^m+3m",  lambda i,m: t**m+3*m)
run("adversarial          f=(1+t)^m-t",lambda i,m: sp.expand((1+t)**m-t))
run("adversarial          f=m!+t^(2m)",lambda i,m: sp.factorial(m)+t**(2*m))
print()
print("CONTROL (satisfiability): an INDEX-AWARE weight hitting 1-t exactly,")
print("  g(h,2)=1, g(1,1)=1, g(h,1)=t   [T1 reads m_1, T2 reads m_h -> separated]")
k=2
def ga(i,m):
    s=i+1
    if (s,m)==(2*k,2): return sp.Integer(1)
    if (s,m)==(1,1):   return sp.Integer(1)
    if (s,m)==(2*k,1): return t
    return sp.Integer(1)
run("index-aware witness", ga)
print()
print("SENTINEL: q430e.py complete.")
