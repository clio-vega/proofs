"""PROVE 11-18 c2 Step 1: the block swap BS on the size-6 census (|lambda| or |mu| = 6, gl4), plus the <=5 census as regression.
usage: python3 bs6.py {6|5} NWORKERS.  Signatures FROZEN (rolesigs.pkl from 11-17c2).  Streams; stores only aggregates + catalogue.
BS(A,B): middles mA=a.C, mB = (C - alpha).(a + gamma) EXACTLY (rigid Z^4 translates), alpha = rh displacement of a, gamma = of C.
Per pair: class, clean-journey BS, formula R2, catalogue K (middles up to translation), resync, lesions (alpha -> alpha +- e_i),
control (BS on both journeys of every B/B non-commuting pair; generic BS on the 'other' journey of B/C)."""
import sys, os, copy, pickle
from collections import Counter, defaultdict
from multiprocessing import Pool
os.chdir('/home/clio/projects/proofs/code-1118c2-Q597')
SIZE = int(sys.argv[1]); NW = int(sys.argv[2]); NCH = int(sys.argv[3]) if len(sys.argv) > 3 else NW
sys.argv = [sys.argv[0], '100000']
exec(open('/home/clio/projects/proofs/code-1118-Q597/widen.py').read().split('st = Counter()')[0])
from fractions import Fraction as Fr
# helpers copied from code-1118-Q597/q597lib.py WITHOUT its module-level align.pkl load (4 GB cgroup; 2 workers were OOM-killed)
def pre(A, B):
    k = 0
    while k < min(len(A), len(B)) and A[k] == B[k]: k += 1
    return k
def suf(A, B, p):
    k = 0
    while k < min(len(A), len(B)) - p and A[-1-k] == B[-1-k]: k += 1
    return k
def rh(S):
    for k, vs in S:
        if k == 'rh': return vs
def cen(vs): return tuple(sum(Fr(v[i]) for v in vs)/len(vs) for i in range(4))
def tr(step, v):
    typ, g, n = step
    f = lambda S: frozenset((k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S)
    return (typ, f(g), f(n))
def pure(S):
    out = {}
    for k, vs in S:
        if k != 'tri': continue
        p0 = frozenset(v[:2] for v in vs); p1 = frozenset(v[2:] for v in vs)
        if len(p1) == 1: out[(0, p0)] = next(iter(p1))
        else: out[(1, p1)] = next(iter(p0))
    return out
def sdelta(step):
    _, g, n = step; a, b = pure(g), pure(n)
    if set(a) != set(b): return None
    return {k: tuple(y-x for x, y in zip(a[k], b[k])) for k in a if a[k] != b[k]}
def addd(D, E):
    out = dict(D)
    for k, v in E.items():
        w = tuple(a+b for a, b in zip(out.get(k, (0, 0)), v)); out[k] = w
        if w == (0, 0): del out[k]
    return out
def neg(d): return {k: tuple(-x for x in v) for k, v in d.items()}
import migrate
_orig = migrate.Board.step
LOG = []
def step_log(self, R, d):
    before = dict(self.tiles); R2 = _orig(self, R, d)
    if R2 is not None:
        gone = {before[t][0:1] + (frozenset(before[t][1]),) for t in set(before) - set(self.tiles)}
        new = {self.tiles[t][0:1] + (frozenset(self.tiles[t][1]),) for t in set(self.tiles) - set(before)}
        LOG.append((self.last, frozenset(gone), frozenset(new)))
    return R2
migrate.Board.step = step_log
def mvl(b, x, T, DV):
    bb = copy.deepcopy(b); LOG.clear(); mig1(bb, 0, x, DV[T]); return bb, list(LOG)
def census_sized():
    for a in range(1, 7):
     for c_ in range(1, 7):
      if SIZE == 6 and a <= 5 and c_ <= 5: continue
      if SIZE == 5 and (a > 5 or c_ > 5): continue
      for lam in parts(a, 4):
       for mu in parts(c_, 4):
        for nu, bs in LR(lam, mu, 4).items():
          if len(bs) < 2: continue
          nt = trim(nu)
          for bb in bs:
            T = beta_cm(bb, lam, mu, nt)
            b0 = canonical(nt, max(nt[0], 1)+1, len(nt)+1); B = b0.nest[1]
            migrate_flock(b0, 1, {B['pos'][(c, r)]: v for (c, r), v in T.items()}, B['u'], B['w'], 'E')
            yield (lam, mu, nu, bb), b0
def disp(s1, s2): return tuple(b-a for a, b in zip(cen(rh(s1[1])), cen(rh(s2[2]))))
def Tsh(d, v):
    out = {}
    for (side, key), val in d.items():
        w = v[:2] if side == 0 else v[2:]
        out[(side, frozenset((a+w[0], b+w[1]) for a, b in key))] = val
    return out
def bs_test(mA, mB, alpha=None):
    if len(mA) < 2 or len(mA) != len(mB): return False
    a, C = mA[0], mA[1:]
    if alpha is None: alpha = disp(a, a)
    gamma = disp(C[0], C[-1]); ma = tuple(-x for x in alpha)
    return all(tr(s, ma) == t for s, t in zip(C, mB[:-1])) and tr(a, gamma) == mB[-1]
def middles(A, B):
    p = pre(A, B); q = suf(A, B, p); return A[p:len(A)-q], B[p:len(B)-q], p, q
E = [tuple((1 if i == j else 0)*sg for i in range(4)) for j in range(4) for sg in (1, -1)]
def work(w):
    st = Counter(); Ks = defaultdict(Counter); bad = []; nstates = 0
    for idx, (key, b0) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b0, cap=CAP):
            nstates += 1
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    for kind in ('BC', 'BB'):
                        if kind == 'BB' and not x < y: continue
                        tx, ty = kind
                        b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                        b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                        P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                        X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                        com = tiles_of(b12) == tiles_of(b21)
                        if kind == 'BB':
                            if com: st['BB com'] += 1; continue
                            st['BB noncom'] += 1
                            for nm, A, B in (('y', Jy, Jyx), ('x', Jx, Jxy)):
                                mA, mB, _, _ = middles(A, B)
                                st[('BB noncom BS', nm, bs_test(mA, mB))] += 1
                                st[('BB noncom resync', nm, bool(mA) and bool(mB) and rh(mA[-1][2]) == rh(mB[-1][2]))] += 1
                            continue
                        K = set(X) | set(Y) | set(Yx) | set(Xy)
                        I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                        sides = {k[0] for k in I}
                        if not I:
                            c = 'BC overlap I=0' if set(X) & set(Y) else 'BC disjoint'
                            st[c] += 1
                            if c == 'BC overlap I=0':
                                for nm, A, B in (('y', Jy, Jyx), ('x', Jx, Jxy)):
                                    mA, mB, _, _ = middles(A, B); st[('I=0 overlap BS', nm, bs_test(mA, mB))] += 1
                            continue
                        Ip, Us = predict(P, X, Y)
                        if len(sides) != 1: st[('BC two-sided I', tuple(sorted(sides)))] += 1; bad.append(('twosided', key, x, y)); continue
                        side = sides.pop(); nU = Us[0][1] if len(Us) == 1 else tuple(Us)
                        st[('BC nontriv', side, nU)] += 1
                        A, B, D1, D2 = (Jy, Jyx, Y, Yx) if side == 0 else (Jx, Jxy, X, Xy)
                        Ao, Bo = (Jx, Jxy) if side == 0 else (Jy, Jyx)
                        mA, mB, p, q = middles(A, B)
                        ok = bs_test(mA, mB)
                        st[('BS clean', side, ok)] += 1
                        mAo, mBo, _, _ = middles(Ao, Bo)
                        st[('BS other (generic)', side, bs_test(mAo, mBo))] += 1
                        st[('resync clean', bool(mA) and rh(mA[-1][2]) == rh(mB[-1][2]))] += 1
                        st[('resync other', bool(mAo) and rh(mAo[-1][2]) == rh(mBo[-1][2]))] += 1
                        if not ok:
                            bad.append(('BSfail', key, x, y, side, nU, mA, mB)); continue
                        a, C = mA[0], mA[1:]; alpha = disp(a, a); gamma = disp(C[0], C[-1])
                        st[('alpha', side, tuple(int(v) if v.denominator == 1 else str(v) for v in alpha))] += 1
                        les = sum(bs_test(mA, mB, tuple(u+e for u, e in zip(alpha, ee))) for ee in E)
                        st[('lesion passes (must be 0)', les)] += 1
                        pred = {}
                        for s in C: pred = addd(addd(pred, Tsh(sdelta(s), tuple(-v for v in alpha))), neg(sdelta(s)))
                        pred = addd(addd(pred, Tsh(sdelta(a), gamma)), neg(sdelta(a)))
                        st[('formula R2', pred == addd(D2, neg(D1)))] += 1
                        o = tuple(-v for v in cen(rh(a[1])))
                        Ks[(side, nU)][(tuple(tr(s, o) for s in mA), tuple(tr(s, o) for s in mB))] += 1
    return st, dict(Ks), bad[:200], nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); Ks = defaultdict(Counter); bad = []; ns = 0
    for s, k, b, n in res:
        st.update(s); bad += b; ns += n
        for kk, v in k.items(): Ks[kk].update(v)
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
    for k in sorted(Ks, key=str): print('catalogue', k, 'configs', len(Ks[k]), 'records', sum(Ks[k].values()))
    pickle.dump(dict(st=st, Ks=dict(Ks), bad=bad), open(f'bs{SIZE}.pkl', 'wb'))
    print('bad examples stored', len(bad))
