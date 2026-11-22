import sys, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
seen=Counter()
for P in PAIRS:
    if P['kind']!='BC': continue
    sides = {k[0] for k in P['I']}
    side = 'empty' if not sides else int(next(iter(sides)))
    if seen[side]>=3: continue
    seen[side]+=1
    print('=== side', side, 'key', P['key'], 'x', P['x'], 'y', P['y'], 'I', sorted((k[0], sorted(k[1])) for k in P['I']))
    for t, A, Bq, Jo in P['roles']:
        mA, mB, p, q = middles(A, Bq)
        print('  ', t, 'p', p, 'A  :', ws(typed(A)))
        print('  ', t, '      rerun:', ws(typed(Bq)), ' |mA|', len(mA), '|mB|', len(mB), 'q', q)
