"""Lemma side as a FINITE check on the configuration table of Theorem AW (code-1120c2-Q597/words.configs):
for every good configuration (target, class, type), the readout change (pure-triangle keys, as bs6.sdelta) is computed
from the removed and created tiles; report the sides touched and whether keys are preserved.  Type a0/a1 from the
rhombus displacement (as diverge.typed).  Control: a deliberately mislabelled table (created tiles translated by a unit
vector of the other block) must break key preservation.  usage: python3 sidecfg.py"""
import sys, os
from fractions import Fraction as Fr
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597')
import words as Wd
def pure(S):
    out = {}
    for k, vs in S:
        if k != 'tri': continue
        p0 = frozenset(v[:2] for v in vs); p1 = frozenset(v[2:] for v in vs)
        if len(p1) == 1: out[(0, p0)] = next(iter(p1))
        else: out[(1, p1)] = next(iter(p0))
    return out
def sdelta(g, n):
    a, b = pure(g), pure(n)
    if set(a) != set(b): return None
    return {k: tuple(y-x for x, y in zip(a[k], b[k])) for k in a if a[k] != b[k]}
def cen(t): return tuple(sum(Fr(P[i]) for P in t[1])/4 for i in range(4))
def shift(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
st = Counter()
for T in ('B', 'C'):
    C = Wd.configs(T)
    for (ck, typ), (R, rem, R2, cre) in sorted(C.items(), key=str):
        al = tuple(b-a for a, b in zip(cen(R), cen(R2)))
        lab = 'b' if typ == 'b' else ('a1' if not any(al[:2]) else ('a0' if not any(al[2:]) else 'a?'))
        d = sdelta(rem, cre)
        sides = None if d is None else tuple(sorted({k[0] for k in d}))
        st[(T, 'engine type', typ, 'label', lab, 'sides', sides, 'nkeys', None if d is None else len(d))] += 1
        # control: translate the created triangles by a unit vector -> keys must change
        bad = frozenset(shift(t, (0, 0, 1, 0)) if t[0] == 'tri' else t for t in cre)
        st[('CONTROL', T, 'keys preserved after shift', sdelta(rem, bad) is not None)] += 1
for k in sorted(st, key=str): print(st[k], k)
