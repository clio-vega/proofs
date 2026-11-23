"""DREAM probe (<= 15 min): is Pi_C - Pi_B invariant along a common prefix?  By Lemma shift (share tex) it is a function
of the rhombus SHAPE only.  Check on every rhombus of every journey in roles5.pkl: is (shape class, Pi_C - Pi_B) a function?
And along each common prefix of (A, re-run) of L: does it change, and only at b-steps?   usage: python3 dinvprobe.py"""
import sys, pickle
from collections import Counter, defaultdict
exec(open('/home/clio/projects/proofs/code-1123-Q597/commonpot.py').read().split('\ndef sub(a, b)')[0])
exec(open('/home/clio/projects/proofs/code-1118-Q597/q597lib.py').read().split('\ndef cen(vs)')[0].split('def rh(S):')[1].join(['def rh(S):', '']) if False else '')
def rh(S):
    for k, vs in S:
        if k == 'rh': return ('rh', vs)
def middles(A, B):
    p = 0
    while p < min(len(A), len(B)) and A[p] == B[p]: p += 1
    return None, None, p, None
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles5.pkl', 'rb'))
def shape(R):
    vs = sorted(R[1]); return frozenset(tuple(b - a for a, b in zip(vs[0], P)) for P in vs)
val = defaultdict(set); st = Counter()
for P in PAIRS:
    for (t, A, Bq, Jo) in P['roles']:
        for J in (A, Bq, Jo):
            for s in J:
                R = rh(s[1]); d = tuple(c - b for b, c in zip(Pi(R, 'B'), Pi(R, 'C'))); val[shape(R)].add(d)
        if P['kind'] != 'BC': continue
        mA, mB, p, q = middles(A, Bq)
        ds = [tuple(c - b for b, c in zip(Pi(rh(s[1]), 'B'), Pi(rh(s[1]), 'C'))) for s in A[:p + 1]]
        chg = [i for i in range(1, len(ds)) if ds[i] != ds[i - 1]]
        st[('prefix: changes', len(chg) > 0, 'all changes after a b-step', all(A[i - 1][0] == 'b' for i in chg))] += 1
print('shapes', len(val), {len(v) for v in val.values()}, 'values per shape:', {str(sorted(k))[:40]: sorted(v) for k, v in val.items()})
for k in sorted(st, key=str): print(st[k], k)
