import pickle, sys
from fractions import Fraction as Fr
D = pickle.load(open('/tmp/c2-diamonds.pkl','rb'))
def ctr(vs): return tuple(sum(Fr(P[k]) for P in vs)/4 for k in range(4))
def disp(h): return [(s[0], tuple(str(a-b) for a,b in zip(ctr(s[2]), ctr(s[3])))) for s in h]
inter = [d for d in D if d['kind']=='BC' and [s[0] for s in d['h1x']]!=[s[0] for s in d['h2x']] or [s[0] for s in d['h1y']]!=[s[0] for s in d['h2y']]]
inter = [d for d in D if d['kind']=='BC' and ([s[1:] for s in d['h1x']]!=[s[1:] for s in d['h2x']] or [s[1:] for s in d['h1y']]!=[s[1:] for s in d['h2y']])]
print(len(inter))
d = min(inter, key=lambda d: len(d['h1x'])+len(d['h1y']))
print(d['key'], d['x'], d['y'], sorted(d['A']), d['l'])
for nm in ('h1x','h2x','h1y','h2y'): print(nm, [s[0] for s in d[nm]]); print('   ', disp(d[nm]))
