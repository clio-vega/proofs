"""PROVE 11-20: per step type, how the rhombus projections s0=pi_0(R), s1=pi_1(R) move in the heights h1 (on Lambda_0) and h0 (on Lambda_1).
Reports (fam, step kind, lattice, (min h before, max h before) - min_before, (min after,max after) - min_before)."""
import sys, pickle, glob
from collections import Counter
sys.argv = [sys.argv[0], '5']; exec(open('potential.py').read().split('st = Counter()')[0])
def rhom(S): return next(t for t in S if t[0] == 'rh')
def kind(g):
    ts = [tri_type(t) for t in g if t[0] == 'tri']
    return ''.join(str(x) for x in sorted(ts))
st = Counter(); sq = Counter()
for fn in sorted(glob.glob('dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in D[fam]:
            for k, (typ, g, n) in enumerate(J):
                Rb, Ra = rhom(g), rhom(n)
                for i in (0, 1):   # i: triangles type whose height lives on Lambda_{1-i}; rhombus projection to Lambda_{1-i}
                    h = hfun(TT, i)
                    hb = [h(x) for x in proj(Rb[1], i)]; ha = [h(x) for x in proj(Ra[1], i)]
                    m = min(hb)
                    st[(fam, TT, typ + kind(g), 'h%d on rh-proj' % i, (min(hb)-m, max(hb)-m), (min(ha)-m, max(ha)-m))] += 1
                # squares: how does the square move in the lift
                Sb = [t for t in g if t[0] == 'sq']; Sa = [t for t in n if t[0] == 'sq']
                for i in (0, 1):
                    h = hfun(TT, i)
                    hs = lambda t: tuple(sorted({h(x) for x in proj(t[1], i)}))
                    sq[(fam, TT, typ+kind(g), 'h%d sq-proj before/after rel. rh min' % i,
                        tuple(x - min(h(y) for y in proj(Rb[1], i)) for x in hs(Sb[0])),
                        tuple(x - min(h(y) for y in proj(Rb[1], i)) for x in hs(Sa[0])))] += 1
for k in sorted(st, key=str): print(st[k], k)
print('---- squares')
for k in sorted(sq, key=str): print(sq[k], k)
