"""PROVE 11-20 c2: Lemma word on real journeys.  For every step of every dumped journey (3881 B/C + 898 B/B-second),
the removed set and the new rhombus must equal those of exactly one abstract configuration (class of R_k, type in the
engine's a/b, a split into a0/a1), translated to R_k.   usage: python3 wordreal.py"""
import pickle, glob
from collections import Counter
import words2 as W
from psilib import rh
st = Counter(); D = {T: W.prep(T) for T in 'BC'}
import words as Wd
for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
    Dm = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in Dm[fam]:
            for ty, g, n in J:
                R = W.it(rh(g)); b = min(R[1]); ck = Wd.classkey(rh(g))
                remreal = frozenset(W.it(t) for t in g if t[0] != 'rh'); Rn = W.it(rh(n))
                hits = []
                for typ in ('a0', 'a1', 'b'):
                    if typ[0] != ty: continue
                    rem, asc, ctl, bn, ckn = D[TT][(ck, typ)]
                    if frozenset(W.tr(t, b) for t in rem) == remreal and min(Rn[1]) == tuple(x+y for x, y in zip(bn, b)): hits.append(typ)
                st[(fam, TT, ty, 'matching abstract configs', tuple(hits))] += 1
                st[('coverage', TT, tuple(sorted(ck)), tuple(hits))] += 1
for k in sorted(st, key=str): print(st[k], k)
