"""Q430 parts 1-6.  See q430.py for the fibre."""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from korff import paths, J_cyclic, mvec
t = sp.Symbol('t')

def c_alpha_M1(k, ell, alpha):
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
        J=J_cyclic(th); assert len(J)==1,(th,J)
        i=J[0]; m=mvec(P[a],w)[i]; assert m>=1
        out.append((i,m))
    return out
def sigma(k,ell,n,omega):
    """sum_T c^-(T) prod_steps omega(idx, m).  omega takes (cyclic index, multiplicity)."""
    h,w=2*k,2*ell+2; tot=sp.Integer(0)
    for P in paths(h,w,tuple([1]*n)):
        c=sign_hkko(P[-1],k,w)
        if c==0: continue
        W=sp.Integer(1)
        for (i,m) in readings(P,w): W*=omega(i,m)
        tot+=c*W
    return sp.expand(tot)

k,ell,n=2,1,4; h,w=2*k,2*ell+2
target=sp.expand(c_alpha_M1(k,ell,tuple([1]*n)))

print("="*78); print("PART 1 -- CONTROLS 1 and 2, as POLYNOMIALS (not labels)")
c1=sigma(k,ell,n, lambda i,m: sp.Integer(1))
print(f" CONTROL 1  omega == 1 (Korff-2019 ribbon arm): sum = {c1}   expected 2 -> {c1==2}")
c2=sigma(k,ell,n, lambda i,m: 1-t**m)
exp2=sp.expand((t-1)**4*(t+1)*(t**2+1)*(t**2+t+2))
print(f" CONTROL 2  omega = 1-t^m (Korff-2011 Psi arm): sum = {sp.factor(c2)}")
print(f"            expected (t-1)^4(t+1)(t^2+1)(t^2+t+2): match -> {sp.expand(c2-exp2)==0}")
print(f"            ord_{{t=1}} = {sp.Poly(c2,t).as_expr().subs(t,1)==0 and min(e for e in range(10) if sp.diff(c2,t,e).subs(t,1)!=0)}")

print(); print("="*78); print("PART 2 -- STRUCTURAL REDUCTION with FREE SYMBOLS (slot-blind weights)")
f1,f2,f3,f4=sp.symbols('f1 f2 f3 f4')
fs={1:f1,2:f2,3:f3,4:f4}
sym=sigma(k,ell,n, lambda i,m: fs[m])
print(f" sum_T c^-(T) W_T  =  {sp.factor(sym)}")
claim=sp.expand(f4*f1**2*(f1+f3))
print(f" CLAIM: = f4*f1^2*(f1+f3).  identity holds -> {sp.expand(sym-claim)==0}")
print(f" NOTE f2 absent: {f2 not in sym.free_symbols}  (the only tableau reading m=2 cancels)")

print(); print("="*78); print("PART 3 -- THE CANDIDATE and the RESIDUAL")
def gb(m): return sp.expand(sum(t**j for j in range(m+1)))     # [m+1]_t
print(f" candidate omega(x,e_i) = [m_{{i-1}}(x)+1]_t ;  [m+1]_t for m=1..4: "
      f"{[sp.expand(gb(m)) for m in range(1,5)]}")
cand=sigma(k,ell,n, lambda i,m: gb(m))
print(f" sum_T c^-(T) W_T = {sp.factor(cand)}")
print(f"                  = {sp.expand(cand)}")
print(f" value at t=1 = {cand.subs(t,1)} ;  deg_t = {sp.Poly(cand,t).degree()}")
print(f" TARGET       = {target} ;  value at t=1 = {target.subs(t,1)} ; deg_t = {sp.Poly(target,t).degree()}")
print(f" EQUAL? {sp.expand(cand-target)==0}")
res=sp.expand(target-cand)
print(f" RESIDUAL  c_alpha - sum = {res}")
print(f"           factored      = {sp.factor(res)}")
print(f"           residual at t=1 = {res.subs(t,1)} ; deg = {sp.Poly(res,t).degree()}")
print(f" via the structural formula: f4*f1^2*(f1+f3) with f1=[2],f3=[4],f4=[5] = "
      f"{sp.expand(claim.subs({f1:gb(1),f3:gb(3),f4:gb(4)}))}  -> agrees {sp.expand(claim.subs({f1:gb(1),f3:gb(3),f4:gb(4)})-cand)==0}")

print(); print("="*78); print("PART 4 -- THEOREM B: FIRST-STEP DIVISIBILITY, as a capacity test")
print(" every path starts theta=e_1 from x=0 and reads m_h(0)=w, so omega(0,e_1) is a")
print(" COMMON FACTOR of every W_T, hence omega(0,e_1) | c_alpha(t) in QQ[t].")
for name,om0 in [("ribbon 2019  omega(0,e_1)=1", sp.Integer(1)),
                 ("Psi 2011     omega(0,e_1)=1-t^w", sp.expand(1-t**w)),
                 ("Gaussian     omega(0,e_1)=[w+1]_t", gb(w))]:
    d=sp.Poly(om0,t).degree() if om0!=1 else 0
    q,r=sp.div(sp.Poly(target,t), sp.Poly(om0,t)) if d>0 else (sp.Poly(target,t), sp.Poly(0,t))
    print(f"   {name:<34} deg={d}  divides c_alpha? {r.as_expr()==0}  (remainder {r.as_expr()})")
print(f"   deg_t c_alpha = {sp.Poly(target,t).degree()} ; so any omega(0,e_1) of degree > "
      f"{sp.Poly(target,t).degree()} is DEAD by degree alone.")
print(f"   integer form: omega(0,e_1)|_{{t=1}} must divide c_alpha(1)={target.subs(t,1)} in ZZ.")
print(f"      Gaussian: [w+1]_t at t=1 = {w+1} ; {w+1} | {target.subs(t,1)} ? "
      f"{int(target.subs(t,1))%(w+1)==0}")

print(); print("="*78); print("PART 5 -- CONTROL 3: a fibre where c_alpha IS constant (k=2, ell=2, n=4)")
k3,ell3,n3=2,2,4; h3,w3=2*k3,2*ell3+2
tg3=sp.expand(c_alpha_M1(k3,ell3,tuple([1]*n3)))
P3=paths(h3,w3,tuple([1]*n3))
nz3=[P for P in P3 if sign_hkko(P[-1],k3,w3)!=0]
allm=[m for P in P3 for (i,m) in readings(P,w3)]
print(f" h={h3} w={w3} n={n3} : n<=w {n3<=w3} ; #tableaux={len(P3)} ; #c^-!=0 = {len(nz3)}")
print(f" c_alpha(t) = {tg3}   [INTEGER: is it constant? {sp.Poly(tg3,t).degree()==0}; deg={sp.Poly(tg3,t).degree()}]")
print(f" c_alpha(1) = {tg3.subs(t,1)}")
print(f" multiplicities actually read: min={min(allm)} max={max(allm)} set={sorted(set(allm))}")
s3a=sigma(k3,ell3,n3, lambda i,m: sp.Integer(1))
s3b=sigma(k3,ell3,n3, lambda i,m: gb(m))
print(f" omega==1  : sum = {s3a}  equals c_alpha? {sp.expand(s3a-tg3)==0}")
print(f" candidate : sum = {sp.expand(s3b)}")
print(f"             equals c_alpha? {sp.expand(s3b-tg3)==0} ; at t=1 = {s3b.subs(t,1)} vs {tg3.subs(t,1)}")
print(f" Theorem B here: omega(0,e_1)=[w+1]_t has deg {w3} ; deg c_alpha = {sp.Poly(tg3,t).degree()} -> dead by degree: {w3>sp.Poly(tg3,t).degree()}")

print(); print("="*78); print("PART 6 -- CONTROL 4, the arm I expect to be redundant: SLOT-AWARE weights")
print(" the cancellation in PART 0 pairs two tableaux with the SAME m-sequence (4,3,2,1)")
print(" but DIFFERENT cyclic-index sequences, so a slot-aware weight can separate them.")
g=sp.symbols('g0:4_1:5')
G={}
gl=sp.symbols(' '.join(f'g{i}_{m}' for i in range(4) for m in range(1,5)))
it=iter(gl)
for i in range(4):
    for m in range(1,5): G[(i,m)]=next(it)
symS=sigma(k,ell,n, lambda i,m: G[(i,m)])
print(f" generic slot-aware sum has {len(symS.free_symbols)} free symbols, "
      f"{len(symS.as_ordered_terms())} terms")
print(f" sum = {symS}")
print(f" does the (4,3,2,1) pair still cancel? substituting slot-blind g(i,m)=f(m): "
      f"{sp.expand(symS.subs({G[(i,m)]: fs[m] for i in range(4) for m in range(1,5)})-claim)==0}")
print(" => the slot-blind reduction is NOT a property of the fibre; it is a property of")
print("    slot-blindness.  Theorem C is correctly scoped; Theorem B is not so scoped.")
print()
print("SENTINEL: q430b.py completed all six parts.")
