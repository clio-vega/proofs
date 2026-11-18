"""PROVE 11-18 step 2. For each B/C pair (x->B, y->C) at M (cap 60 sampler), record y's journey in M and in Mx step by step:
each step = (removed tiles, added tiles) as frozensets of (kind, vertex-tuple-set).  Align: common prefix / common suffix / middle.
Also classify pairs by SCEL |U| (from widen.predict). Pickle per-pair journeys for |U|=1 and for B/B control."""
import sys, copy, pickle
sys.argv = [sys.argv[0], '60'] if len(sys.argv) < 2 else sys.argv
exec(open('widen.py').read().split('st = Counter()')[0])
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
out = []
for key, b0 in census():
    for b, DV in states(b0, cap=CAP):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind, (tx, ty) in (('BC', 'BC'), ('BB', 'BB')):
                    if kind == 'BB' and not x < y: continue
                    b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                    b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                    P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                    X, Y = dl(P, P1), dl(P, P2)
                    Ip, Us = predict(P, X, Y)
                    out.append(dict(kind=kind, key=key, x=x, y=y, Us=Us, Ip=Ip, Jx=Jx, Jy=Jy, Jyx=Jyx, Jxy=Jxy,
                                    M=dict(b.tiles), com=tiles_of(b12) == tiles_of(b21), P=P, X=X, Y=Y,
                                    Yx=dl(P1, P12), Xy=dl(P2, P21)))
pickle.dump(out, open('align.pkl', 'wb'))
print(len(out))
