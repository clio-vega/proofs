"""T4. THE BRIDGE.  Claim (derived by hand; this is the load-bearing step):

For T in the class (K,J,eps) -- K = pairs touched, J = unpaired elements used,
a_i = eps_i v_i -- the term B(T)A(T) factorises as

    B(T)A(T) = C(K,J) * prod_{i<i' in K} R(a_i,a_i')^{-1}
                      * prod_{i in K, p notin J} R(a_i, -u_p)^{-1},
    R(a,b) = (ab-1)(a+b),

with C(K,J) INDEPENDENT of eps.  Test: the quotient must be constant over the orbit.

Second claim: R(a,b) = a*b*(sigma(a)+sigma(b)),  sigma(x) = x - 1/x.  (identity, exact)
"""
import itertools, random
from fractions import Fraction as F
from eng import admissible, terms

def R(a, b): return (a*b-1)*(a+b)
def sig(x):  return x - 1/x

# --- claim 2 first: a pure identity
bad2 = 0
for _ in range(400):
    a = F(random.randint(-30,30), random.randint(1,7)); b = F(random.randint(-30,30), random.randint(1,7))
    if a == 0 or b == 0: continue
    if R(a,b) != a*b*(sig(a)+sig(b)): bad2 += 1
print(f"T4a  R(a,b) = ab(sigma(a)+sigma(b)) : {'OK' if bad2==0 else f'{bad2} FAILURES'}")

# --- claim 1: the eps-invariance of C(K,J)
def build(v, u):
    M = []; pairidx = []
    for x in v: pairidx.append((len(M), len(M)+1)); M += [x, -x]
    uidx = list(range(len(M), len(M)+len(u))); M += list(u)
    return M, pairidx, uidx

def classify(T, pairidx, uidx):
    Ts = set(T); K = []; eps = []
    for i,(p,q) in enumerate(pairidx):
        if p in Ts and q in Ts: return None
        if p in Ts: K.append(i); eps.append(+1)
        elif q in Ts: K.append(i); eps.append(-1)
    J = tuple(sorted(j for j,ix in enumerate(uidx) if ix in Ts))
    return tuple(K), J, tuple(eps)

random.seed(4242)
ok = bad = 0; fails = []
cells = [([F(2),F(3)],[],2), ([F(2),F(3),F(5)],[],3), ([F(2),F(3),F(5),F(7)],[],4),
         ([F(2),F(3),F(5)],[F(11,3)],3), ([F(2),F(3),F(5)],[F(11,3)],2),
         ([F(2),F(3),F(5),F(7)],[F(11,3)],3), ([F(2),F(3),F(5)],[F(11,3),F(13,5)],4),
         ([F(2),F(3),F(5)],[F(11,3),F(13,5)],3), ([F(5,2),F(7,3),F(9,4),F(11,5)],[F(8,3),F(13,6)],4),
         ([F(2),F(3)],[F(11,3),F(13,5),F(17,7)],3)]
for v, u, r in cells:
    M, pairidx, uidx = build(v, u)
    if not admissible(M, require_recip_free=True): print("   skip (not recip-free)", v, u); continue
    theta = list(range(r))
    cls = {}
    for T, B, A in terms(M, theta):
        c = classify(T, pairidx, uidx)
        if c is None: continue
        K, J, eps = c
        a = {K[idx]: (eps[idx]*v[K[idx]]) for idx in range(len(K))}
        pred = F(1)
        for i, ip in itertools.combinations(K, 2): pred *= R(a[i], a[ip])
        for i in K:
            for p in range(len(u)):
                if p not in J: pred *= R(a[i], -u[p])
        C = (B*A)*pred            # should be eps-independent
        cls.setdefault((K,J), set()).add(C)
    for (K,J), s in cls.items():
        if len(s) == 1: ok += 1
        else:
            bad += 1
            if len(fails) < 3: fails.append((v,u,r,K,J,sorted(s)[:3]))
print(f"T4b  C(K,J) is eps-independent : {ok}/{ok+bad} classes   (10 cells, lam<=4, f<=3)")
for f_ in fails: print("   FAIL", f_)
