"""For the class-1-a role L of every diverging B/C pair: BS(L); alpha in one block Lambda_s; delta_a on side s;
the 11-18 Thm red formula  Delta_formula = sum_C (T_{-alpha} d_sigma - d_sigma) + T_gamma d_a - d_a  versus the observed
mixed difference Delta = Xy - X; case a (re-run a_T'): types of C;  case b: Delta != 0.
Lesion: the same formula on E (the other role).   usage: python3 lformula.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
def sub(D, E):
    out = dict(D)
    for k, v in E.items():
        w = tuple(a-b for a, b in zip(out.get(k, (0, 0)), v)); out[k] = w
        if w == (0, 0): del out[k]
    return out
def formula(mA):
    a, C = mA[0], mA[1:]
    alpha = disp(a, a); gamma = disp(C[0], C[-1]); ma = tuple(-x for x in alpha)
    F = {}
    for s in C:
        d = sdelta(s); F = addd(F, sub(Tsh(d, ma), d))
    da = sdelta(a); F = addd(F, sub(Tsh(da, gamma), da))
    return F, alpha, da
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    sides = {k[0] for k in P['I']}
    side = 'empty' if not sides else ('two' if len(sides) > 1 else int(next(iter(sides))))
    Delta = sub(P['Xy'], P['X'])
    info = []
    for i, (t, A, Bq, Jo) in enumerate(P['roles']):
        mA, mB, p, q = middles(A, Bq); info.append((typed(A)[p], typed(Bq)[p], mA, mB, typed(A)[p:len(A)-q]))
    L = [i for i in (0, 1) if info[i][0] == (1, aT[i])][0]
    for i in (0, 1):
        tA, tB, mA, mB, TT = info[i]; rl = 'L' if i == L else 'E(lesion)'
        F, alpha, da = formula(mA)
        s = L  # side index of L: role0 x->B a1 -> side 1 ; role1 y->C a0 -> side 0
        sL = 1 if L == 0 else 0
        blk = 'L1' if not any(alpha[:2]) else ('L0' if not any(alpha[2:]) else 'mixed')
        st[(rl, side, 'rerun', tB[1], 'BS', bs_test(mA, mB), 'formula==Delta', F == Delta, 'alpha block', blk, 'd_a sides', tuple(sorted({k[0] for k in da})))] += 1
        if rl == 'L':
            st[('L', side, 'rerun', tB[1], 'Delta!=0', bool(Delta), 'supp side', tuple(sorted({k[0] for k in Delta})), 'side(L)', sL)] += 1
            if tB[1] != 'b': st[('L case a: C types', tuple(TT[1:]))] += 1
for k in sorted(st, key=str): print(st[k], k)
