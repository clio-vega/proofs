# Rigorous certificate: for b=2,3 mod4, b<=N, the ONLY rational roots of I_b(m) in Z[m]
# are the forced {0,1,...,R}, R=floor((b-1)/2).  -> Q_b has no rational root -> (star).
# Uses exact rational-root test (Poly.ground_roots / nroots-free).  No tail bound needed.
import sympy as sp
from qb import I_b
m=sp.symbols('m')
import sys
N=int(sys.argv[1]) if len(sys.argv)>1 else 50
bad=[]
for b in range(3, N+1):
    if b%4 not in (2,3): continue
    R=(b-1)//2
    P=sp.Poly(I_b(b), m)
    # exact rational roots with multiplicity
    rr=P.ground_roots()  # dict {root: mult} of rational roots
    rat=set(sp.nsimplify(r) for r in rr.keys())
    forced=set(sp.Integer(r) for r in range(0,R+1))
    extra = rat - forced
    if extra:
        bad.append((b, extra))
        print(f"b={b}: EXTRA RATIONAL ROOTS {extra}  (forced {sorted(forced)})")
    else:
        if b%30 in (2,3,6,7,10,11) or b>N-5:
            print(f"b={b:3d}(mod4={b%4}): rational roots = forced {sorted(int(x) for x in rat)} only.  OK")
print("CERTIFIED no extra rational root for all checked b<=%d: %s"%(N, "YES" if not bad else f"NO {bad}"))
