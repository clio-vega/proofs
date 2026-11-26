"""Provenance for the corridor and back window (case b), relative indices; L's labels shifted so that the BS(L) partner is visible."""
import pickle, glob
from collections import Counter
from prov import trace
st = Counter()
for fn in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(fn, 'rb')):
        if case != 'b': continue
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        k = laL - 3
        HE0, HL0 = set(AE[pE][1]), set(AL[pL][1])
        def lab(t, made):
            if t in made: return made[t]
            return 'M0'
        for j in list(range(2, k+2))[:1] + [k+2, k+3]:
            made = trace([('L', AL, pL), ('EL', ArE[:pE+j], pE)])
            nm = 'corr' if j < k+2 else 'k+%d' % (j - k)
            labs = tuple(sorted(lab(t, made).replace('L+%d' % (k+2), 'L+(k+2)').replace('L+%d' % (k+1), 'L+(k+1)') if k > 0 else lab(t, made) for t in ArE[pE+j][1]))
            st[('E^L_' + nm, 'k>0' if k else 'k=0', labs)] += 1
        for j in [k+2, k+3]:
            made = trace([('E', AE[:pE+j], pE)])
            st[('E_k+%d' % (j-k), tuple(sorted('E%+d' % (int(lab(t, made)[1:]) - k) if lab(t, made).startswith('E') else 'M0' for t in AE[pE+j][1])))] += 1
for k in sorted(st, key=str): print(st[k], k)
