"""TASK A: verify and adopt Rick's n=2 domain-legal witness (defect D3 repair).

Two MECHANISMS, not two handles:
  M1  Warnaar Eq_generalise determinant engine (symmetric functions, e_r(x^pm),
      theta-sum over y in Z^k) -> c_alpha(t) as a polynomial.
  M2  HKKO signed lattice-path count (prop:cssyt=path + c^-_{2h,w}) -> c_alpha(1).
They share NO code path: M1 never builds a path, M2 never builds a determinant.
"""
import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from walks import AB_from_walks
from korff import paths, J_cyclic, weight_path
t = sp.Symbol('t')

def ordt1(f):
    f = sp.expand(f)
    if f == 0: return None
    d = 0
    while d <= 30:
        if sp.simplify(sp.diff(f, t, d).subs(t, 1)) != 0: return d
        d += 1
    return ">30"

def c_alpha_M1(k, ell, alpha):
    """Mechanism 1: coefficient of x^alpha in (x_1..x_n)^k * G."""
    n = len(alpha)
    E = Engine(n)
    F = G_times_prodx(E, k, ell)
    mon = sp.prod([E.x[i]**alpha[i] for i in range(n)])
    return sp.expand(sp.Poly(F, *E.x).coeff_monomial(mon) or sp.Integer(0))

def c_alpha_M2(k, ell, alpha):
    """Mechanism 2: A - B over HKKO paths (= c_alpha(1))."""
    A,B,TOT = AB_from_walks(k, 2*ell+2, alpha)
    return A-B, A, B, TOT

def nu_data(k, ell, alpha):
    h = 2*k; w = 2*ell+2
    out = []
    for P in paths(h, w, alpha):
        nu = sum(len(J_cyclic(tuple(P[a+1][i]-P[a][i] for i in range(h))))
                 for a in range(len(P)-1))
        out.append((P[-1], nu, sp.expand(weight_path(P, w))))
    return out

def sign_hkko(la, k, w):
    h = 2*k
    if all(la[2*i]==la[2*i+1] for i in range(k)): return 1
    if la[0]-la[h-1]==w and all(la[2*i+1]==la[2*i+2] for i in range(k-1)): return -1
    return 0

print("="*78)
print("TASK A.1  THE CRITICAL POINT  k=2, ell=1, n=2, alpha=(1,1)")
print("          domain: n_K = h_cyc = 2k = 4 > 2 (OK);  ell_K = n = 2 <= k_K = w = 4 (OK)")
print("="*78)
k,ell,alpha = 2,1,(1,1)
h,w = 2*k, 2*ell+2
nd = nu_data(k,ell,alpha)
print(f"  fibre size (#cylindric tableaux of content {alpha} in h={h}, w={w}) = {len(nd)}")
for la,nu,psi in nd:
    print(f"    endpoint {la}  nu(T)={nu}  c^-={sign_hkko(la,k,w):+d}  Psi_T={sp.factor(psi)}  Psi_T(1)={psi.subs(t,1)}")
print(f"  min nu over fibre = {min(x[1] for x in nd)}   <= the obstruction's lower bound")
c1 = c_alpha_M1(k,ell,alpha)
c2,A,B,TOT = c_alpha_M2(k,ell,alpha)
print(f"  M1 (Warnaar determinant): c_alpha(t) = {c1}   c_alpha(1) = {sp.expand(c1).subs(t,1)}   ord_t=1 = {ordt1(c1)}")
print(f"  M2 (HKKO signed paths)  : A={A} B={B} TOT={TOT}  c_alpha(1) = A-B = {c2}")
print(f"  TWO MECHANISMS AGREE at t=1: {sp.expand(c1).subs(t,1) == c2}")
korffsum = sp.expand(sum(sign_hkko(la,k,w)*psi for la,nu,psi in nd))
print(f"  Korff side  sum_T c^-(T) Psi_T(t) = {sp.factor(korffsum)}  ; at t=1 = {korffsum.subs(t,1)}")
print(f"  IDENTITY HOLDS? {sp.simplify(korffsum-c1)==0}")
