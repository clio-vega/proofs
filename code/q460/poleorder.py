"""
Q460 -- the pole-order engine.

CLAIM: for arcs u,v on disjoint column sets (tie-free),
   ord_{y=x} Delta_{U,V}(x,y;z)  =  -1  if the closed intervals [u1,u2],[v1,v2] are DISJOINT
                                 =  -2  if they OVERLAP (nested or crossing).
This is what lets a pole order read off the interval-intersection graph of sigma.
We also tabulate the collision cases, where the rule is EXPECTED to break.
"""
import sympy as sp, sys
sys.path.insert(0, '/home/clio/projects/code/q460')
from ikdelta import Delta, pattern, q

x, t = sp.symbols('x t')

def ord_at_y_eq_x(u, v):
    """Order of Delta_{u,v}(x, x+t) in t at t=0 (negative = pole)."""
    e = Delta(u, v, x, x + t)
    e = sp.cancel(sp.together(e))
    num, den = sp.fraction(e)
    on = sp.Poly(sp.expand(num), t).monoms()
    od = sp.Poly(sp.expand(den), t).monoms()
    return min(m[0] for m in on) - min(m[0] for m in od)

def overlap(u, v):
    return not (u[1] < v[0] or v[1] < u[0])

print("TIE-FREE pairs (no shared column): is ord = -(1 + [overlap]) ?")
cases = []
MAXC = 6
arcs = [(i, j) for i in range(1, MAXC+1) for j in range(i, MAXC+1)]
checked = agree = 0
by_pat = {}
for u in arcs:
    for v in arcs:
        su = set([u[0], u[1]]); sv = set([v[0], v[1]])
        if u[0] != u[1] and v[0] != v[1] and (su & sv):  continue   # collision
        if u[0] == u[1] and v[0] == v[1] and u[0] == v[0]: continue # impossible col sum
        if u[0] == u[1] and v[0] != v[1] and u[0] in sv:  continue  # col sum 3
        if v[0] == v[1] and u[0] != u[1] and v[0] in su:  continue
        if u == v: continue
        try:
            o = ord_at_y_eq_x(u, v)
        except ValueError:
            continue
        pred = -(1 + (1 if overlap(u, v) else 0))
        checked += 1
        if o == pred: agree += 1
        else: print("   MISMATCH u=%s v=%s ord=%d pred=%d" % (u, v, o, pred))
        p = pattern(u, v) or pattern(v, u)
        nm = (p[0] + ('' if pattern(u,v) else '^T'))
        by_pat.setdefault(nm, set()).add(o)
print("   %d/%d agree   (stratum: all tie-free ordered arc pairs on <=%d columns)"
      % (agree, checked, MAXC))
print("   pole order by pattern: %s" % {k: sorted(v) for k, v in sorted(by_pat.items())})

print("\nCOLLISION pairs (shared column) -- the rule is expected to BREAK:")
seen = {}
for u in arcs:
    for v in arcs:
        if u[0] == u[1] or v[0] == v[1]: continue
        su = set(u); sv = set(v)
        if not (su & sv): continue
        if u == v: continue
        # column sums must be <= 2: a shared column holds one 1 from each row: ok
        try:
            o = ord_at_y_eq_x(u, v)
        except ValueError:
            continue
        p = pattern(u, v)
        nm = (p[0] if p else (pattern(v, u)[0] + '^T'))
        seen.setdefault(nm, set()).add(o)
for k in sorted(seen):
    print("   %-6s ord_{y=x} = %s   (overlap always TRUE here; rule would predict -2)"
          % (k, sorted(seen[k])))
