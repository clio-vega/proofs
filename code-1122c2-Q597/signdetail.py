"""detail on signrows: (p0,p1) values, types per role, ordered by later/earlier, vs side. usage: python3 signdetail.py SIZE"""
import sys, os, pickle
from collections import Counter
rows = pickle.load(open('signrows%s%s.pkl' % (sys.argv[1], os.environ.get('TAG', '')), 'rb'))
st = Counter()
for kind, key, x, y, side, r in rows:
    st[('P-DIFF', kind, side, r[0]['p'] - r[1]['p'])] += 1
    st[('P-VALUES', kind, side, r[0]['p'], r[1]['p'])] += 1
    L = 0 if r[0]['p'] > r[1]['p'] else 1
    lat, ear = r[L], r[1 - L]
    st[('LATER role', kind, side, 'later=role%d' % L, 'later', lat['tA'], lat['tB'], 'S', lat['S'], 'good', lat['good'], '| earlier', ear['tA'], ear['tB'], 'S', ear['S'], 'good', ear['good'])] += 1
for k in sorted(st, key=str): print(st[k], k)
