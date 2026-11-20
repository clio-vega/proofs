"""PROVE 11-20 Step 0: candidate potential for no-repeat, in the product picture (Lemma hex, 11-16).
For a journey to nest T in {B,C}, the 11-16 direction law says each moved 0-triangle raises h0(q) by 1 and each
moved 1-triangle raises h1(p) by 1, where (barycentric x=(a,b,-a-b)):  B: h0=-x2, h1=y2 ;  C: h0=x3, h1=-y3.
Candidate potential P0: m0(k) = min of h0 over pi_1(V(H_k)) on steps k that move 0-triangles; claim m0 weakly DEcreasing.
(similarly m1 over pi_0(V(H_k)) on steps moving 1-triangles).  If so, no-repeat follows: a triangle moved at j sits at m0(j)+1 > m0(k).
Also checks: moved triangle starts at min h over pi(H) and ends at min+1 (the 'width 1' fact).
usage: python3 potential.py 5 [lesion]   (lesion = reverse the step order of each journey)"""
import sys, pickle, glob
from collections import Counter
from fractions import Fraction as Fr
SIZE = sys.argv[1]; LESION = len(sys.argv) > 2 and sys.argv[2] == 'lesion'
def hfun(T, i):
    # returns height on a lattice point (2-tuple) of Lambda_{1-i} for i-triangles
    if T == 'B': return (lambda q: -q[1]) if i == 0 else (lambda p: p[1])
    return (lambda q: -q[0]-q[1]) if i == 0 else (lambda p: q_neg3(p))
def q_neg3(p): return p[0]+p[1]          # -y3 = -(-a-b) = a+b
def tri_type(t):
    vs = t[1]; a = {v[:2] for v in vs}; b = {v[2:] for v in vs}
    return 0 if len(b) == 1 else (1 if len(a) == 1 else None)
def other(t, i):  # the constant projection of an i-triangle
    v = next(iter(t[1])); return v[2:] if i == 0 else v[:2]
def proj(vs, i): return {(v[2:] if i == 0 else v[:2]) for v in vs}   # projection to the OTHER lattice of i-triangles
st = Counter(); bad = []
for fn in sorted(glob.glob('dumpj%s-*.pkl' % SIZE)):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in D[fam]:
            if LESION: J = [(typ, n, g) for typ, g, n in reversed(J)]
            ms = {0: [], 1: []}; viol = False
            for k, (typ, g, n) in enumerate(J):
                V = set().union(*[t[1] for t in g])
                for i in (0, 1):
                    h = hfun(TT, i)
                    mg = [t for t in g if t[0] == 'tri' and tri_type(t) == i]
                    mn = [t for t in n if t[0] == 'tri' and tri_type(t) == i]
                    if not mg: continue
                    hs = [h(x) for x in proj(V, i)]; lo, hi = min(hs), max(hs)
                    st[(fam, 'width of pi(H) in h', i, hi-lo)] += 1
                    st[(fam, 'moved start at min', i, all(h(other(t, i)) == lo for t in mg))] += 1
                    st[(fam, 'moved end at min+1', i, all(h(other(t, i)) == lo+1 for t in mn))] += 1
                    ms[i].append(lo)
            for i in (0, 1):
                s = ms[i]; mono = all(s[a] >= s[a+1] for a in range(len(s)-1))
                st[(fam, 'm%d weakly decreasing' % i, mono)] += 1
                if not mono: viol = True
            st[(fam, 'journeys')] += 1
for k in sorted(st, key=str): print(st[k], k)
