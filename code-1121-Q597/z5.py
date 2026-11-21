"""PROVE 11-21 retreat: (z5), the Romik-Sniady-type rim test, made operational.
slide k = J_o (the other journey, on M0); slide k+1 = the clean journey re-run on Mx_P (= M0.J_o after the common prefix).
  Ywake := non-rhombus tiles of Mx_P removed by the re-run's middle mB (not created by mB itself)  -- what slide k+1 reads.
  route := rhombus positions R^o_j of J_o (j = 0..end).  rim(S) := tiles of board S sharing an EDGE with some route position.
Tests, per record:  (z5-start) Ywake + alpha  ⊆ rim(M0)        (translate back to the start board)
                    (z5-end)   Ywake          ⊆ rim(M0.J_o)    (no translation)
                    (z5-wake)  Ywake + alpha  ⊆ W ∪ Xadd       (inside slide k's wake, the WS(X) mirror)
Side: for t in rim, psi_{T'}(t) - psi_{T'}(R^o_j) at an adjacent route position (REM side = componentwise <= 0).
Control: B/B pairs (clean = BBy and BBx roles), same recipe.   ONLYC=1 (env): slide k+1 wake = the C-alpha steps only (mB[:-1]); the last step a+gamma is translated by gamma (CT2).
usage: python3 z5.py 5 NW NCH"""
import os
ONLYC = os.environ.get("ONLYC", "0") == "1"
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
from stepset import tedges
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
def psi(t, TT):
    eta = ETA[TT]
    return (min(-eta(P[2:]) for P in t[1]), min(eta(P[:2]) for P in t[1]))
def rtile(s, which=1):
    return ('rh', frozenset(rh(s[which])))
def rim_of(board, route):
    RE = {}
    for r in route:
        for e in tedges(r): RE.setdefault(e, []).append(r)
    out = {}
    for t in board:
        if t[0] == 'rh': continue
        hit = [r for e in tedges(t) for r in RE.get(e, ())]
        if hit: out[t] = hit
    return out
def T(S, v): return {(k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S}
def work(w):
    st = Counter(); nstates = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            nstates += 1; M0 = bare(b)
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    for kind in ('BC', 'BB'):
                        if kind == 'BB' and not x < y: continue
                        tx, ty = kind
                        b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                        b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                        if kind == 'BB':
                            if tiles_of(b12) == tiles_of(b21): continue
                            roles = (('BBy', Jy, Jyx, Jx, 'B'), ('BBx', Jx, Jxy, Jy, 'B'))
                        else:
                            P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                            X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                            K = set(X) | set(Y) | set(Yx) | set(Xy)
                            I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                            sides = {k[0] for k in I}
                            if len(sides) != 1: continue
                            side = sides.pop()
                            roles = (('clean->C', Jy, Jyx, Jx, 'B'),) if side == 0 else (('clean->B', Jx, Jxy, Jy, 'C'),)
                        for tag, A, B, Jo, To in roles:
                            mA, mB, p, q = middles(A, B)
                            if kind == 'BC' and not bs_test(mA, mB): continue
                            if len(mA) < 2 or len(mB) < 1: continue
                            al = disp(mA[0], mA[0])
                            Mx0 = apply(M0, Jo); MxP = apply(Mx0, B[:p])
                            made = set(); Yw = set()
                            for s in (mB[:-1] if ONLYC else mB):
                                for u in s[1]:
                                    if u[0] != 'rh' and u not in made: Yw.add(u)
                                made = (made - set(s[1])) | set(s[2])
                            Yw &= set(MxP)
                            route = [rtile(s) for s in Jo] + [rtile(Jo[-1], 2)]
                            rimS, rimE = rim_of(set(M0), route), rim_of(set(Mx0), route)
                            Wset = set(M0) - set(Mx0); Xadd = set(Mx0) - set(M0)
                            MP = apply(M0, A[:p]); wakeP = set(MP) - set(MxP)
                            st[(tag, 'z5-WS(X): Ywake+alpha in M_P minus Mx_P', T(Yw, al) <= wakeP)] += 1
                            Yt = T(Yw, al)
                            st[(tag, 'records')] += 1
                            st[(tag, 'z5-start: Ywake+alpha in rim(M0)', Yt <= set(rimS))] += 1
                            st[(tag, 'z5-end: Ywake in rim(M0.Jo)', Yw <= set(rimE))] += 1
                            st[(tag, 'z5-wake: Ywake+alpha in W u Xadd', Yt <= (Wset | Xadd))] += 1
                            st[(tag, 'Ywake+alpha in W', Yt <= Wset)] += 1
                            for t in Yt & set(rimS):
                                sd = set()
                                for r in rimS[t]:
                                    d = tuple(a - c for a, c in zip(psi(t, To), psi(r, To))); sd.add(d[0] <= 0 and d[1] <= 0)
                                st[(tag, 'side of translated wake tile in rim(M0): REM side', tuple(sorted(sd)))] += 1
                            st[(tag, 'frac of Ywake+alpha in rim(M0)', round(len(Yt & set(rimS)) / max(1, len(Yt)), 1))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
