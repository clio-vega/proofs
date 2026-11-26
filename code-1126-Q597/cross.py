"""For each step of E's re-run middle: which steps of (E journey, L journey, L re-run) it is a rigid translate of, and by what vector."""
import pickle, glob
from collections import Counter
from match import canon, vsub
st = Counter(); shown = 0
for f in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(f, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        srcs = {'E': AE, 'L': AL, 'Lr': ArL}
        mids = {'E': range(pE, pE+laE), 'L': range(pL, pL+laL), 'Lr': range(pL, pL+lbL)}
        row = []
        for j in range(pE, pE+lbE):
            shr, mr = canon(ArE[j]); hits = []
            for nm, J in srcs.items():
                for i in mids[nm]:
                    sh, m = canon(J[i])
                    if sh == shr: hits.append('%s%d:%s' % (nm, i - (pE if nm == 'E' else pL), ''.join('%+d' % c for c in vsub(mr, m))))
            row.append('|'.join(hits) or 'NONE')
        st[(case, 'any NONE', 'NONE' in row)] += 1
        if shown < 12 and case == 'b':
            shown += 1; print(case, TE, row)
for k in sorted(st, key=str): print(st[k], k)
