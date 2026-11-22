"""Orientation by TYPE (not by step index): exactly one role at a class-1 a-step of its target; that role's re-run type;
dg; side of I.  Role-swap lesion: call the OTHER role L.  usage: python3 orient.py SIZE"""
import sys, pickle
from collections import Counter
rows = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/signrows%s.pkl' % sys.argv[1], 'rb'))
aT = {'BC': ('a1', 'a0'), 'BB': ('a1', 'a1')}
sideof = {'BC': (1, 0), 'BB': (None, None)}
st = Counter()
for kind, key, x, y, side, r in rows:
    c1a = [i for i in (0, 1) if r[i]['tA'] == (1, aT[kind][i])]
    st[('PRE #class-1-a roles', kind, len(c1a))] += 1
    if len(c1a) != 1: continue
    L = c1a[0]
    for lesion in (False, True):
        Lx = 1 - L if lesion else L
        rb = r[Lx]['tB'][1]
        pred = (sideof[kind][Lx] if rb == 'b' else 'empty') if kind == 'BC' else ('BB-good' if rb == 'b' else 'BB-?')
        obs = side if kind == 'BC' else ('BB-good' if r[Lx]['good'] else 'BB-notgood')
        st[('LESION' if lesion else 'ORIENT', kind, 'L=role%d' % Lx, 'rerun', rb, 'dg', r[0]['g'] - r[1]['g'], 'pred', pred, 'obs', obs, 'agree', pred == obs)] += 1
for k in sorted(st, key=str): print(st[k], k)
