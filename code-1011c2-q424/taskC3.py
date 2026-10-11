"""TASK C consolidated: the ribbon weight of Korff 1906.02565, screened and refuted
by a NEW mechanism (t-blindness), plus the sharpened specification."""
import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from korff import paths, J_cyclic, weight_path
t,q = sp.symbols('t q')

def c_alpha_M1(k,ell,alpha):
    n=len(alpha); E=Engine(n); F=G_times_prodx(E,k,ell)
    return sp.expand(sp.Poly(F,*E.x).coeff_monomial(sp.prod([E.x[i]**alpha[i] for i in range(n)])) or sp.Integer(0))
def sign_hkko(la,k,w):
    h=2*k
    if all(la[2*i]==la[2*i+1] for i in range(k)): return 1
    if la[0]-la[h-1]==w and all(la[2*i+1]==la[2*i+2] for i in range(k-1)): return -1
    return 0

def ribbon_wt_strip(ncomp, ncols_per_comp):
    """Korff 1906.02565 (Hecke_wt) specialised to a HORIZONTAL strip: every connected
    component h is a single row, so r(h)=1 and wt = (t-1)^{p-1} prod_j t^{c_j-1}."""
    p = ncomp
    return sp.expand((t-1)**(p-1) * sp.prod([t**(c-1) for c in ncols_per_comp]))

print("1. RIBBON WEIGHT ON A HORIZONTAL STRIP (r(h)=1 for every component):")
print("   wt = (t-1)^{p-1} * prod_j t^{c_j-1},  p = #components, c_j = #columns of comp j")
for p,cs in [(1,[1]),(1,[2]),(1,[3]),(2,[1,1]),(2,[2,1]),(3,[1,1,1])]:
    W = ribbon_wt_strip(p,cs)
    print(f"   p={p} cols={cs}: wt = {str(sp.factor(W)):<22} wt(1) = {W.subs(t,1)}")
print("   => SINGLE BOX (p=1, c=1): wt = 1 IDENTICALLY, for every t.  t-FREE.")
print("   => PASSES the t=1 screen (weight-must-be-one-at-t1 needs W_T(1)=1): YES, 1.")
print("   => pre-registered kill ('another prod (1-t^m) over a cyclic index set') DOES NOT FIRE.")
print()
print("2. BUT: on alpha=(1^n) every step is a single box, so W_T = 1 for EVERY T,")
print("   hence sum_T c^-(T) W_T = sum_T c^-(T) = c_alpha(1), a CONSTANT in t.")
print()
print("3. THE REFUTATION, in Korff's domain (h=2k>2 and n<=w):")
for (k,ell,n) in [(2,1,4),(3,1,4)]:
    h,w=2*k,2*ell+2; al=tuple([1]*n)
    P=paths(h,w,al)
    sgn = sum(sign_hkko(p[-1],k,w) for p in P)
    c = c_alpha_M1(k,ell,al)
    korff = sp.expand(sum(sign_hkko(p[-1],k,w)*weight_path(p,w) for p in P))
    print(f"   k={k} ell={ell} n={n} (h={h}>2, n={n}<=w={w}):  #T={len(P)}")
    print(f"      ribbon architecture : sum_T c^- W_T = {sgn}  (t-free)")
    print(f"      TRUTH               : c_alpha(t)    = {sp.expand(c)}   [Warnaar engine]")
    print(f"      ribbon FAILS: {sp.simplify(sp.expand(c)-sgn)!=0}   (agrees only at t=1: {sp.expand(c).subs(t,1)==sgn})")
    print(f"      Korff 2011 Psi      : sum_T c^- Psi_T = {sp.factor(korff)} ; at t=1 = {korff.subs(t,1)}  -> ord too HIGH")
print()
print("4. q^d TWIST CANNOT RESCUE IT (q = Novikov/quantum parameter of qH^*(Gr), t-free):")
k,ell,n=2,1,4; h,w=2*k,2*ell+2
P=paths(h,w,tuple([1]*n))
tw=sp.expand(sum(sign_hkko(p[-1],k,w)*q**p[-1][0] for p in P))
c=c_alpha_M1(k,ell,tuple([1]*n))
print(f"   sum_T c^-(T) q^{{d(T)}} = {tw}        target = {sp.expand(c)}")
poly=sp.Poly(sp.expand(tw-c), q)
print(f"   as a polynomial in q: {poly.as_expr()}")
print(f"   irreducible over QQ(t)? {sp.factor(poly.as_expr())}")
print(f"   constant q in QQ solving it: {[s for s in sp.solve(sp.Eq(tw,c), q) if s.is_number]}  (none, since LHS is t-free and RHS is not)")
print()
print("5. THE SHARPENED SPECIFICATION (what BOTH failures jointly say).")
print("   On the fibre alpha=(1^n) the single-box weight omega(x,e_i) must satisfy:")
print("     (a) omega(x,e_i)|_{t=1} != 0      [else sum_T c^- W_T(1) = 0 != c_alpha(1)]")
print("     (b) omega(x,e_i) NOT t-free       [else sum_T c^- W_T is t-free != 3 - t^2]")
print("   Korff 2011 Psi  : omega = (1-t^{m_i}), vanishes at t=1  -> violates (a)")
print("   Korff 2019 wt   : omega = 1,           t-free           -> violates (b)")
print("   A Gaussian binomial [m+1 choose 1]_t = (1-t^{m+1})/(1-t) satisfies BOTH:")
for m in range(1,5):
    gb = sp.expand(sp.cancel((1-t**(m+1))/(1-t)))
    print(f"      m={m}: [m+1 choose 1]_t = {gb}  value at t=1 = {gb.subs(t,1)} != 0, and non-constant")
