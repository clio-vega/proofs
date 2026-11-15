"""Locality census: for each diamond, compare the x-journey alone (h1x) with x-after-y (h2x), and y alone (h2y) with y-after-x (h1y).
Report the step-sequence of rhombus displacements, LCP/LCS lengths and whether the displacement MULTISET is preserved."""
import pickle
from fractions import Fraction as Fr
from collections import Counter
D = pickle.load(open('/tmp/c2-diamonds.pkl','rb'))
def ctr(vs): return tuple(sum(Fr(P[k]) for P in vs)/4 for k in range(4))
def disp(h): return [tuple(a-b for a,b in zip(ctr(s[2]), ctr(s[3]))) for s in h]
def lcp(a,b):
    k=0
    while k<min(len(a),len(b)) and a[k]==b[k]: k+=1
    return k
st = Counter()
for d in D:
    for alone, after in (('h1x','h2x'), ('h2y','h1y')):
        A, B = d[alone], d[after]
        dA, dB = disp(A), disp(B)
        same_seq = dA == dB
        st[(d['kind'], 'same disp seq', same_seq)] += 1
        st[(d['kind'], 'same disp multiset', Counter(dA)==Counter(dB))] += 1
        if not same_seq:
            p = lcp(dA,dB); s = lcp(dA[::-1], dB[::-1])
            st[(d['kind'], 'middle lens', len(dA)-p-s, len(dB)-p-s)] += 1
for k in sorted(st, key=str): print(k, st[k])
