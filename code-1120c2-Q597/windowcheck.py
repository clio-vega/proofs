"""PROVE 11-20 c2: window-lemma ingredients from ENGINE configurations (words3.py config extraction): Q_T, Delta(mu,lambda), stationary steps, and the L1 bound 5.\n(config extraction copied from words3.py:)
Configurations are taken from the ENGINE (the dumped real journeys), not from localenum: for each (T, class, step) the
removed non-rh tiles and the new rhombus, relative to the base point of R (asserted identical across all occurrences).
The class is the unordered pair {+-u} of R's Lambda_0 edge; the step is identified by the engine type plus the new
rhombus offset (a0 vs a1 distinguished by which projection moves).  Ascent is recomputed here from first principles:
v in Lambda_0 is ascending iff eta_T(v) = +1; v in Lambda_1 iff eta_T(v) = -1.  Iterative DFS, own data structures.
usage: python3 words3.py L"""
import sys, pickle, glob
from collections import Counter, defaultdict
ETA = {'B': lambda a, b: b, 'C': lambda a, b: a + b}
def pts(t): return tuple(sorted(tuple(int(x) for x in P) for P in t[1]))
def rhom(S): return next(t for t in S if t[0] == 'rh')
def cls(R):
    P = pts(R); s0 = sorted({p[:2] for p in P}); d = (s0[1][0]-s0[0][0], s0[1][1]-s0[0][1])
    return frozenset({d, (-d[0], -d[1])})
def shift(P, b): return tuple(x - y for x, y in zip(P, b))
CONF = defaultdict(set)
for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for T, M0, J in D[fam]:
            for ty, g, n in J:
                R = rhom(g); b = min(pts(R)); Rn = rhom(n)
                rem = frozenset((t[0], tuple(sorted(shift(P, b) for P in pts(t)))) for t in g if t[0] != 'rh')
                nb = shift(min(pts(Rn)), b)
                s0moves = {p[:2] for p in pts(R)} != {p[:2] for p in pts(Rn)}
                s1moves = {p[2:] for p in pts(R)} != {p[2:] for p in pts(Rn)}
                step = ty + ('0' if ty == 'a' and s0moves and not s1moves else ('1' if ty == 'a' else ''))
                CONF[(T, cls(R), step)].add((rem, nb, cls(Rn)))
bad = [k for k, v in CONF.items() if len(v) != 1]
print('engine configurations', len(CONF), 'non-unique', bad)
assert not bad and len(CONF) == 12
CONF = {k: next(iter(v)) for k, v in CONF.items()}

def psi(P, T): return (min(-ETA[T](p[2], p[3]) for p in P), min(ETA[T](p[0], p[1]) for p in P))
def A(P): return (min(p[0] for p in P), min(p[2] for p in P))
def om(P, R, T): a, b = psi(P, T), psi(R, T); c, d = A(P), A(R); return (a[0]-b[0], a[1]-b[1], c[0]-d[0], c[1]-d[1])
def canonR(c):
    d = sorted(c)[0]; return [(0,0,0,0), (d[0],d[1],0,0), (0,0,d[0],d[1]), (d[0],d[1],d[0],d[1])]
QS = {}
for (T, c, s), (rem, nb, cn) in sorted(CONF.items(), key=str):
    R = canonR(c); R = [tuple(x - y for x, y in zip(p, min(R))) for p in R]
    Rn = [tuple(x + y for x, y in zip(p, nb)) for p in canonR(cn)]; Rn = [tuple(x - y for x, y in zip(p, min(canonR(cn)))) for p in canonR(cn)]
    Rn = [tuple(x + y for x, y in zip(p, nb)) for p in Rn]
    dl = om(Rn, R, T)
    print(T, sorted(c), s, 'Delta(mu,lambda)', dl, 'stationary' if dl == (0,0,0,0) else '')
    for k, P in rem: QS.setdefault(T, set()).add(om(P, R, T))
for T in 'BC': print(T, 'Q_T from engine', sorted(QS[T]))
