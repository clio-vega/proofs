"""
Q460 -- the order-vector table.  Orders of Delta_{U,V}(x,y;z) along the three
divisors  {y=x}, {y=qx}, {x=qy}, for each tie-free type.
If the map type -> order-vector is INJECTIVE, a pole order reads off the type.
"""
import sympy as sp, sys
sys.path.insert(0, '/home/clio/projects/code/q460')
from ikdelta import Delta, pattern, q

x, t = sp.symbols('x t')

def order_along(u, v, sub):
    e = sp.cancel(sp.together(Delta(u, v, x, sub)))
    num, den = sp.fraction(e)
    def lo(p):
        p = sp.expand(p)
        if p == 0: return None
        return min(m[0] for m in sp.Poly(p, t).monoms())
    a, b = lo(num), lo(den)
    return (a - b) if (a is not None and b is not None) else None

# representatives of each tie-free type on 6 columns
REP = {
 'A'  : ((1,2),(3,4)),    # U entirely before V
 'A^T': ((3,4),(1,2)),
 'B'  : ((1,4),(2,3)),    # V nested in U
 'B^T': ((2,3),(1,4)),
 'C'  : ((1,3),(2,4)),    # crossing
 'C^T': ((2,4),(1,3)),
}
DIV = [("y=x   ", x + t), ("y=qx  ", q*x + t), ("x=qy  ", x/q + t)]
print("order of Delta along each divisor (negative = pole):")
print("  type   " + "  ".join(d[0] for d in DIV))
vecs = {}
for nm, (u, v) in REP.items():
    row = tuple(order_along(u, v, s) for _, s in DIV)
    vecs[nm] = row
    print("  %-6s %s" % (nm, "   ".join("%+d" % r if r is not None else "  ." for r in row)))
inv = len(set(vecs.values())) == len(vecs)
print("\n  injective on the 6 labels? %s" % inv)
print("  distinct order-vectors: %d  (labels: %d)" % (len(set(vecs.values())), len(vecs)))
# B vs B^T: is case B literally a symmetric function?
dB1 = Delta((1,4),(2,3), sp.Symbol('X'), sp.Symbol('Y'))
dB2 = Delta((1,4),(2,3), sp.Symbol('Y'), sp.Symbol('X'))
print("  is case B symmetric, Delta(x,y)=Delta(y,x)?  %s"
      % (sp.simplify(dB1 - dB2) == 0))
for nm, vv in sorted(vecs.items()):
    print("    %-4s -> %s" % (nm, vv))
