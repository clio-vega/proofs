"""3-factors around (0,b) in every journey of the census (A, re-run, J_o of every role; target from the role tag):
(prev step, (0,b), next) with next classified straight / other / END.  Is '(0,a?) (0,b) straight' absent?
usage: python3 factor3.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE = int(sys.argv[1]); sys.argv = ['x', str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
STR = {'B': (1, 'a1'), 'C': (1, 'a0')}
st = Counter(); seen = set()
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
for P in PAIRS:
    if P['kind'] != 'BC': continue
    for (tag, A, Bq, Jo) in P['roles']:
        T = 'B' if tag.startswith('x') else 'C'; Tp = 'C' if T == 'B' else 'B'
        for J, TT in ((A, T), (Bq, T), (Jo, Tp)):
            key = (TT, tuple(tuple(sorted(map(str, s[1]))) for s in J))
            if key in seen: continue          # count each distinct journey once
            seen.add(key); w = typed(J)
            for i in range(len(w)):
                if w[i] != (0, 'b'): continue
                prev = w[i-1] if i > 0 else 'START'
                nxt = 'END' if i + 1 == len(w) else ('STRAIGHT' if w[i+1] == STR[TT] else w[i+1])
                if prev != 'START' and prev[0] == 0 and prev[1] != 'b': prev = (0, 'a', 'same-as-straight-letter' if prev[1] == STR[TT][1] else 'other-letter')
                st[(TT, prev, nxt)] += 1
for k in sorted(st, key=str): print(st[k], k)
