"""Second handshake: do E_1 and L_1 remove the same M0 tiles outside both contact hexagons? (reads /tmp/step0_5_*.pkl from step0.py)"""
import pickle, glob
from collections import Counter
st = Counter()
for fn in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(fn, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, _, _ = info[E]; TL, AL, ArL, pL, _, _ = info[L]
        HE0, HL0 = set(AE[pE][1]), set(AL[pL][1])
        e1, l1 = set(AE[pE+1][1]), set(AL[pL+1][1])
        made_e = set(t for s in AE[:pE+1] for t in s[2]); made_l = set(t for s in AL[:pL+1] for t in s[2])
        oth = lambda S, made: {t for t in S if t not in made and t not in HE0 and t not in HL0}
        oe, ol = oth(e1, made_e), oth(l1, made_l)
        st[(case, '|other(E1)|', len(oe), '|other(L1)|', len(ol), '|other(E1) & other(L1)|', len(oe & ol),
            '|created(E1) & created(L1)|', len(set(AE[pE+1][2]) & set(AL[pL+1][2])))] += 1
for k in sorted(st, key=str): print(st[k], k)
