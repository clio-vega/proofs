"""PROVE 11-20: replay each journey; at each step record (T, kind, rhombus class u, #candidates, tie-break used,
#candidates of each 'direction' (sign of eta-change of moved triangles)).  Tests whether the plane rule alone forces the direction law."""
import sys, pickle, glob, cmath
from collections import Counter
sys.argv = [sys.argv[0], '5', '1', '1']
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590'); from common import unit_dir
from stepset import candidates, step as sstep
import os; os.chdir('/home/clio/projects/proofs/code-1120-Q597')
DVEC = {'B': unit_dir((0, 0, 0, 1)), 'C': unit_dir((0, -1, 0, 0))}
eta = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
def ttype(t):
    a = {v[:2] for v in t[1]}; b = {v[2:] for v in t[1]}
    return 0 if len(b) == 1 else 1
def rclass(R):
    pts = sorted({v[:2] for v in R[1]}); d = tuple(b-a for a, b in zip(*pts)); return max(d, tuple(-x for x in d))
def tridelta(S, c, k):
    # direction of the moved triangles under rotation of S about c by k*60: returns tuple of (type, eta-sign of displacement of the 'other' coordinate)
    from migrate import prot
    out = []
    for t in S:
        if t[0] != 'tri': continue
        i = ttype(t); v = next(iter(t[1]))
        w = tuple(a+b for a, b in zip(c, prot(tuple(Fr(p)-q for p, q in zip(v, c)), k)))
        dv = tuple(b-a for a, b in zip(v, w)); dd = dv[2:] if i == 0 else dv[:2]
        out.append((i, dd))
    return out
st = Counter()
for fn in sorted(glob.glob('dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for TT, M0, J in D['BC']:
        Bd = set(M0); d = DVEC[TT]
        for (typ, g, n) in J:
            R = next(t for t in g if t[0] == 'rh')
            cs = candidates(Bd, R, d)
            kind = typ + ''.join(sorted(str(ttype(t)) for t in g if t[0] == 'tri'))
            # classify each candidate of type a by the h-sign of its triangle moves
            sig = []
            for (ty, S, c) in cs:
                if ty == 'a':
                    mv = tridelta(S, c, 3)
                    i, dd = mv[0]
                    hs = -eta[TT](dd) if i == 0 else eta[TT](dd)
                    sig.append(('a%d' % i, hs))
                else: sig.append(('b', None))
            chosen = next((ty, S, c) for (ty, S, c) in cs if S == g)
            st[(TT, kind, 'u', rclass(R), '#cands', len(cs), 'sigs', tuple(sorted(sig, key=str)))] += 1
            Bd -= g; Bd |= n
for k in sorted(st, key=str): print(st[k], k)
