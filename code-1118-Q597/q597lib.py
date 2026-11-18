"""Shared helpers for PROVE 11-18 Q597 analyses of align.pkl (exact Fractions)."""
import pickle
from fractions import Fraction as Fr
R = pickle.load(open('/home/clio/projects/proofs/code-1118-Q597/align.pkl', 'rb'))
def pre(A, B):
    k = 0
    while k < min(len(A), len(B)) and A[k] == B[k]: k += 1
    return k
def suf(A, B, p):
    k = 0
    while k < min(len(A), len(B)) - p and A[-1-k] == B[-1-k]: k += 1
    return k
def cls(r):
    if r['kind'] == 'BB': return 'BB com' if r['com'] else 'BB noncom'
    if r['Us']: return 'BC U=' + '+'.join(f"{s}:{n}" for s, n in r['Us'])
    return 'BC I=0 overlap' if set(r['X']) & set(r['Y']) else 'BC disjoint'
def rh(S):
    for k, vs in S:
        if k == 'rh': return vs
def cen(vs): return tuple(sum(Fr(v[i]) for v in vs)/len(vs) for i in range(4))
def tr(step, v):
    typ, g, n = step
    f = lambda S: frozenset((k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S)
    return (typ, f(g), f(n))
def tvec(s1, s2):
    if s1[0] != s2[0]: return None
    v = tuple(b-a for a, b in zip(cen(rh(s1[1])), cen(rh(s2[1]))))
    return v if tr(s1, v) == s2 else None
def pure(S):
    out = {}
    for k, vs in S:
        if k != 'tri': continue
        p0 = frozenset(v[:2] for v in vs); p1 = frozenset(v[2:] for v in vs)
        if len(p1) == 1: out[(0, p0)] = next(iter(p1))
        else: out[(1, p1)] = next(iter(p0))
    return out
def sdelta(step):
    _, g, n = step; a, b = pure(g), pure(n)
    if set(a) != set(b): return None
    return {k: tuple(y-x for x, y in zip(a[k], b[k])) for k in a if a[k] != b[k]}
def addd(D, E):
    out = dict(D)
    for k, v in E.items():
        w = tuple(a+b for a, b in zip(out.get(k, (0, 0)), v)); out[k] = w
        if w == (0, 0): del out[k]
    return out
def neg(d): return {k: tuple(-x for x in v) for k, v in d.items()}
def middles(A, B):
    p = pre(A, B); q = suf(A, B, p); return A[p:len(A)-q], B[p:len(B)-q], p, q
