"""Analyse bs6.pkl vs bs5.pkl: catalogue per (side,|U|), subset relation, new configs, grammar a.b.w.b with w in {A,B}* (A=(1,-1;0,0), B=(0,0;1,0))."""
import pickle, re, sys
sys.path.insert(0, '../code-1118-Q597')
from q597lib import cen, rh
from collections import Counter
D5 = pickle.load(open('bs5.pkl', 'rb')); D6 = pickle.load(open('bs6.pkl', 'rb'))
def disp(s): return tuple(b-a for a, b in zip(cen(rh(s[1])), cen(rh(s[2]))))
RUN = {(1, -1, 0, 0), (0, 0, 1, 0)}
st = Counter()
allk = sorted(set(D5['Ks']) | set(D6['Ks']), key=str)
for k in allk:
    c5 = set(D5['Ks'].get(k, {})); c6 = set(D6['Ks'].get(k, {}))
    print(k, '<=5:', len(c5), ' size6:', len(c6), ' 5-subset-of-6:', c5 <= c6, ' new at 6:', len(c6 - c5), ' union:', len(c5 | c6))
for D, nm in ((D5, '<=5'), (D6, '6')):
    for k, C in D['Ks'].items():
        for (mA, mB), cnt in C.items():
            w = ''.join(s[0] for s in mA)
            gram = re.fullmatch(r'aba*b', w) is not None and all(disp(s) in RUN for s in mA[2:-1])
            st[(nm, 'grammar a.b.w.b with run vectors', gram)] += 1
            lens = len(mA)
            st[(nm, 'max middle length', k[0])] = max(st[(nm, 'max middle length', k[0])], lens)
for k in sorted(st, key=str): print(st[k], k)
