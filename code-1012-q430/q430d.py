"""Q430 PART 8: the two-tableau fibre alpha=(1,1) at ell=0 (w=2).
CLAIM: for EVERY slot-blind multiplicative weight omega(x,e_i)=f(m_{i-1}(x)),
  sum_T c^-(T) W_T = 0  identically, while c_alpha(t) = 1-t != 0.
This refutes the slot-blind ARCHITECTURE, not one candidate."""
import sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from korff import paths, J_cyclic, mvec
t=sp.Symbol('t')
def c_alpha_M1(k,ell,alpha):
    n=len(alpha); E=Engine(n); F=G_times_prodx(E,k,ell)
    return sp.expand(sp.Poly(F,*E.x).coeff_monomial(sp.prod([E.x[i]**alpha[i] for i in range(n)])) or sp.Integer(0))
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

print("PART 8 -- THE TWO-TABLEAU FIBRE:  alpha=(1,1), ell=0 => w=2, any k>=1")
print()
for k in (1,2,3,4,5):
    h,w=2*k,2
    P=paths(h,w,(1,1))
    c=sp.expand(c_alpha_M1(k,0,(1,1)))
    print(f" k={k}  h={h}  w={w}   in Korff's domain (h>2 and n<=w): {h>2 and 2<=w}")
    print(f"    c_alpha(t) = {c}      c_alpha(1) = {c.subs(t,1)}   deg = {sp.Poly(c,t).degree() if c!=0 else -1}")
    print(f"    #cylindric tableaux = {len(P)}")
    tot_blind = sp.Integer(0); tot_aware = sp.Integer(0)
    f = sp.Function('f'); g = sp.Function('g')
    for Pp in P:
        rd = readings(Pp,w); sg = sign_hkko(Pp[-1],k,w)
        slots = tuple(j[0]+1 for j in rd); ms = tuple(j[1] for j in rd)
        Wb = sp.prod([f(m) for (i,m) in rd])
        Wa = sp.prod([g(i+1,m) for (i,m) in rd])
        print(f"      endpoint {str(Pp[-1]):<16} c^- = {sg:+d}   m-index seq {slots}  m-seq {ms}"
              f"   W_blind = {Wb}")
        tot_blind += sg*Wb; tot_aware += sg*Wa
    tot_blind = sp.expand(tot_blind); tot_aware = sp.expand(tot_aware)
    print(f"    SLOT-BLIND  sum_T c^- W_T = {tot_blind}      identically zero? {tot_blind==0}")
    print(f"    needs to equal c_alpha   = {c}            possible? {tot_blind==c}")
    print(f"    SLOT-AWARE  sum_T c^- W_T = {tot_aware}   identically zero? {tot_aware==0}")
    print()

print("-"*78)
print("CROSS-CHECK: the three named weights, all slot-blind, at k=2, w=2, alpha=(1,1):")
k,w=2,2
def run(name, om):
    tot=sp.Integer(0)
    for Pp in paths(2*k,w,(1,1)):
        sg=sign_hkko(Pp[-1],k,w)
        W=sp.prod([om(i,m) for (i,m) in readings(Pp,w)])
        tot+=sg*W
    tot=sp.expand(tot)
    print(f"   {name:<40} sum = {tot}   target 1-t  -> match {tot==sp.expand(1-t)}")
run("Korff 2019 ribbon   f(m)=1",        lambda i,m: sp.Integer(1))
run("Korff 2011 Psi      f(m)=1-t^m",    lambda i,m: 1-t**m)
run("Gaussian binomial   f(m)=[m+1]_t",  lambda i,m: sp.expand(sum(t**j for j in range(m+1))))
run("adversarial random  f(m)=t^m+3m",   lambda i,m: t**m+3*m)
run("adversarial random  f(m)=(1+t)^m-t",lambda i,m: sp.expand((1+t)**m-t))
print()
print("CONTROL: a SLOT-AWARE weight that DOES hit 1-t (shows the fibre is satisfiable,")
print("         so PART 8 tests slot-blindness and not the fibre):")
# want g(h,2)*g(1,1) - g(h,2)*g(h,1) = 1-t ; take g(h,2)=1, g(1,1)=1, g(h,1)=t
def ga(i,m):
    if (i+1,m)==(2*k,2): return sp.Integer(1)
    if (i+1,m)==(1,1):   return sp.Integer(1)
    if (i+1,m)==(2*k,1): return t
    return sp.Integer(1)
run("slot-aware  g(h,2)=1, g(1,1)=1, g(h,1)=t", ga)
print()
print("SENTINEL: q430d.py part 8 complete.")
