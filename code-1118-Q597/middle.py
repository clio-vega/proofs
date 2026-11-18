"""Inspect the divergent middles of y's journey (M vs Mx) on one-star pairs. Tiles shown as (kind, centroid in Z^4 as Fractions*3)."""
import sys, pickle
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
from fractions import Fraction as Fr
from collections import Counter
exec(open('analyse_align.py').read().split('T = Counter()')[0])
def cen(t):
    k, vs = t; n = len(vs); return k[0] + str(tuple(sum(Fr(v[i]) for v in vs)*12/n for i in range(4))).replace('Fraction', '')
def show(J, a, b):
    for typ, g, n in J[a:b]:
        print('   ', typ, 'gone', sorted(cen(t) for t in g), '\n        new ', sorted(cen(t) for t in n))
one = [r for r in R if r['kind'] == 'BC' and len(r['Us']) == 1 and r['Us'][0][1] == 1]
print(len(one), Counter(r['Us'][0][0] for r in one))
shown = Counter()
for r in one:
    side = r['Us'][0][0]
    if shown[side] >= 1: continue
    shown[side] += 1
    A, B = r['Jy'], r['Jyx']; p = pre(A, B); q = suf(A, B, p)
    print('=== side', side, 'len', len(A), len(B), 'prefix', p, 'suffix', q)
    print('  x-journey types', [t for t, _, _ in r['Jx']])
    print('  y in M  middle:'); show(A, p, len(A)-q)
    print('  y in Mx middle:'); show(B, p, len(B)-q)
