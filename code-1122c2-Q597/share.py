"""(Share): H^L_{pL} cap H^E_{pE} cap M0.  Kinds and count of shared tiles, by L's local type.  usage: python3 share.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
aT = {'BC': ('a1', 'a0'), 'BB': ('a1', 'a1')}; st = Counter()
for P in PAIRS:
    kind = P['kind']; M0 = set(P['M0']); H = []; T = []
    for i, (t, A, Bq, Jo) in enumerate(P['roles']):
        mA, mB, p, q = middles(A, Bq); H.append(set(A[p][1])); T.append((typed(A)[p], typed(Bq)[p]))
    L = [i for i in (0, 1) if T[i][0] == (1, aT[kind][i])][0]
    sh = H[0] & H[1] & M0
    st[(kind, 'L', T[L], 'E', T[1-L], 'shared kinds', tuple(sorted(u[0] for u in sh)), 'shared nonM0', len((H[0] & H[1]) - M0))] += 1
for k in sorted(st, key=str): print(st[k], k)
