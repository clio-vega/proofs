"""PROVE 11-20 c2: SECOND, independent implementation of the Lemma AW enumeration.
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
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
def ascending_translates(rem, Rbase_cls, T):
    """ascending R-edge vectors: the class vector placed in Lambda_0 with eta=+1, and in Lambda_1 with eta=-1"""
    out = set()
    for d in Rbase_cls:
        if ETA[T](*d) == 1: v = (d[0], d[1], 0, 0); out |= {(k, tuple(sorted(tuple(x+y for x, y in zip(P, v)) for P in Q))) for k, Q in rem}
        if ETA[T](*d) == -1: v = (0, 0, d[0], d[1]); out |= {(k, tuple(sorted(tuple(x+y for x, y in zip(P, v)) for P in Q))) for k, Q in rem}
    return frozenset(out)
def place(S, b): return frozenset((k, tuple(sorted(tuple(x+y for x, y in zip(P, b)) for P in Q))) for k, Q in S)
L = int(sys.argv[1]); st = Counter()
for T in 'BC':
    classes = {k[1] for k in CONF if k[0] == T}; steps = sorted({k[2] for k in CONF if k[0] == T})
    for c0 in classes:
        stack = [((0, 0, 0, 0), c0, ())]          # (base of R, class, history of (rem, asc))
        while stack:
            b, c, hist = stack.pop()
            for s in steps:
                rem, nb, cn = CONF[(T, c, s)]
                R_rem = place(rem, b); R_asc = place(ascending_translates(rem, c, T), b)
                st['nodes'] += 1
                if R_asc & R_rem: st['SAME-STEP VIOLATION'] += 1
                for (r, a) in hist:
                    st['pairs'] += 1
                    if R_asc & r or a & R_rem: st['AW VIOLATION'] += 1
                if len(hist) + 1 < L:
                    stack.append((tuple(x+y for x, y in zip(b, nb)), cn, hist + ((R_rem, R_asc),)))
print('L', L)
for k in sorted(st, key=str): print(st[k], k)
