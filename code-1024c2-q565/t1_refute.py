"""T1. Candidate counterexamples to conj:transversal (Psi=0 <=> r >= f+2).

Reasoning: thm:r1 gives Psi^C_{(1),(0)}(M) = 1 - e^{-2}, e = prod_t w_t.  So at r=1,
theta=(0) the zero locus is e = +-1 -- which has nothing to do with negation pairs.
If a reciprocal-free M of lambda negation pairs and f unpaired elements can have
prod w = +-1, conj:transversal's "=>" direction dies at r=1, inside its own swept
domain (lambda<=4, f<=3, n<=7, r in {f+1,f+2,f+3}).
"""
from fractions import Fraction as F
from eng import psi_pair, psi_head, admissible, terms

def show(name, w, theta, lam, f):
    rf = admissible(w, require_recip_free=True)
    a = psi_pair(w, theta); b = psi_head(w, theta)
    e = F(1)
    for x in w: e *= x
    nz = sum(1 for _, B, A in terms(w, theta) if B*A != 0)
    tot = sum(1 for _ in terms(w, theta))
    r = len(theta)
    print(f"{name}")
    print(f"   M = {[str(x) for x in w]}   n={len(w)}  lambda={lam}  f={f}  r={r}  theta={theta}")
    print(f"   admissible={admissible(w)}  reciprocal-free={rf}   prod w = {e}")
    print(f"   psi_pair = {a}   psi_head = {b}   agree={a==b}")
    print(f"   nonzero terms B(T)A(T): {nz}/{tot}")
    print(f"   conj:transversal predicts zero iff r>={f+2}  ->  predicts {'ZERO' if r>=f+2 else 'NONZERO'}; measured {'ZERO' if a==0 else 'NONZERO'}")
    print()

# (1) lambda=3, f=0, r=1 : prod v_i = -1 with no pairwise product = +-1
v = [F(2), F(3), F(1,6)]
M = []
for x in v: M += [x, -x]
show("C1  lambda=3,f=0,r=1  (v = 2,3,1/6 ; prod v = 1)", M, [0], 3, 0)

# (2) lambda=0, f=2, r=1 : w_1 w_2 = -1  (an ANTI-reciprocal pair)
show("C2  lambda=0,f=2,r=1  (w_1 w_2 = -1)", [F(2), F(-1,2)], [0], 0, 2)

# (3) lambda=1, f=1, r=1 : v,-v,u with v^2 u = +-1
#     prod = -v^2 u ; want = +-1 -> u = +-1/v^2
v1 = F(3)
show("C3  lambda=1,f=1,r=1  (u = 1/v^2)", [v1, -v1, F(1,9)], [0], 1, 1)

# (4) lambda=2, f=1, r=1
v1, v2 = F(2), F(5)
u = F(1, 100)   # prod = (v1 v2)^2 u = 100 u = 1
show("C4  lambda=2,f=1,r=1  (u = 1/(v1 v2)^2)", [v1, -v1, v2, -v2, u], [0], 2, 1)
