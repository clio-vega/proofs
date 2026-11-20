"""PROVE 11-20 c2: the abstract single-lattice moves.  From the board-free table (localenum), for every GOOD configuration
record what the step does in Lambda_0: s0 -> s0', and the 0-cells swept (pi_0 of the removed 0-triangles); and dually
in Lambda_1.  usage: python3 absmoves.py"""
import sys
sys.path.insert(0, '../code-1120-Q597')
from fractions import Fraction as Fr
import localenum as LE
def kind(t):
    if t[0] != 'tri': return t[0]
    return 'tri0' if len({P[2:] for P in t[1]}) == 1 else 'tri1'
def moves(T):
    out = []
    for u in LE.UNITS:
        if LE.ETA[T](u) == 0: continue
        R, Hs = LE.classify(T, u)
        for typ, H in Hs:
            r = LE.analyse(T, R, typ, set(H))
            if not r['xrule'] or r['k'] is None or r['across'] is False: continue
            V = set().union(*[t[1] for t in H]); zs = LE.convex_hull([LE.phi(P) for P in V])
            corners = [P for P in V if LE.is_corner(LE.phi(P), zs)]
            c = tuple(sum(Fr(P[i]) for P in corners)/6 for i in range(4))
            Hn = frozenset((t[0], frozenset(LE.rot_about(P, c, r['k']) for P in t[1])) for t in H)
            Rn = next(t for t in Hn if t[0] == 'rh')
            rec = {'u': u, 'typ': typ, 'rot': r['k']}
            for i in (0, 1):
                rec['s%d' % i] = frozenset(tuple(int(x) for x in P[2*i:2*i+2]) for P in R[1])
                rec['s%d_new' % i] = frozenset(tuple(int(x) for x in P[2*i:2*i+2]) for P in Rn[1])
                rec['cells%d' % i] = [frozenset(tuple(int(x) for x in P[2*i:2*i+2]) for P in t[1]) for t in H if kind(t) == 'tri%d' % i]
            out.append(rec)
    return out
if __name__ == '__main__':
    for T in 'BC':
        for m in moves(T):
            print(T, m['u'], m['typ'], m['rot'], 's0', sorted(m['s0']), '->', sorted(m['s0_new']), 'cells0', [sorted(c) for c in m['cells0']],
                  '| s1', sorted(m['s1']), '->', sorted(m['s1_new']), 'cells1', [sorted(c) for c in m['cells1']])
