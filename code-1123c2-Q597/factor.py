"""Is the factor (0,b)(1,a_T) [straight step of the journey's own target right after a class-0 b] forbidden in journeys?
Census over every journey in roles5.pkl / roles6.pkl (A, re-run, J_o of every role; target from the role tag), and in
the abstract word model (all words, words.py configs: is the factor even placeable?).   usage: python3 factor.py"""
import sys, pickle
from collections import Counter
sys.argv = ['x', '5', '1', '1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
STR = {'B': (1, 'a1'), 'C': (1, 'a0')}
st = Counter()
for SIZE in (5, 6):
    PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
    for P in PAIRS:
        if P['kind'] != 'BC': continue
        for (tag, A, Bq, Jo) in P['roles']:
            T = 'B' if tag.startswith('x') else 'C'; Tp = 'C' if T == 'B' else 'B'
            for nm, J, TT in (('journey', A, T), ('rerun', Bq, T), ('J_o', Jo, Tp)):
                w = typed(J)
                for i in range(len(w) - 1):
                    if w[i] == (0, 'b'): st[('size', SIZE, nm, 'after (0,b):', 'own straight' if w[i+1] == STR[TT] else w[i+1])] += 1
                if w[-1] == (0, 'b'): st[('size', SIZE, nm, 'ends right after (0,b)')] += 1
for k in sorted(st, key=str): print(st[k], k)
