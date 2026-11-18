"""Summarise align.pkl: for y's journey in M (Jy) vs in Mx (Jyx): common prefix p, common suffix q, middle lengths.
Classes: BC by SCEL |U| (one-star = Us==[(side,1)]), BC overlapping with I empty, BC disjoint, BB commuting / non-commuting."""
import pickle
from collections import Counter
R = pickle.load(open('align.pkl', 'rb'))
def pre(A, B):
    k = 0
    while k < min(len(A), len(B)) and A[k] == B[k]: k += 1
    return k
def suf(A, B, p):
    k = 0
    while k < min(len(A), len(B)) - p and A[-1-k] == B[-1-k]: k += 1
    return k
def tiles(J):
    out = set()
    for _, g, n in J: out |= g | n
    return out
def cls(r):
    if r['kind'] == 'BB': return 'BB com' if r['com'] else 'BB noncom'
    overl = bool(set(r['X']) & set(r['Y']))
    if r['Us']: return 'BC U=' + '+'.join(f"{s}:{n}" for s, n in r['Us'])
    return 'BC I=0 overlap' if overl else 'BC disjoint'
T = Counter(); tab = {}
for r in R:
    c = cls(r)
    for name, A, B in (('y', r['Jy'], r['Jyx']), ('x', r['Jx'], r['Jxy'])):
        p = pre(A, B); q = suf(A, B, p)
        T[(c, name, 'mid', len(A)-p-q, len(B)-p-q)] += 1
for k in sorted(T, key=str): print(T[k], k)
