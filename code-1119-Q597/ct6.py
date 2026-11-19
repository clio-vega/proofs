"""PROVE 11-19 Steps 0+1, size 6 (or 5 as regression): usage python3 ct6.py {5|6} NW [NCH].
(0) I=0 overlap: WHY bs_test is False (empty middles / len<2 / length mismatch / translation mismatch / identical journeys).
(1) CT on every nontrivial B/C pair (probe5 logic): C-candidates in M_i translate by -alpha to candidates in Mx_i; last cands = a cands + gamma.
    Kinds of non-rh candidate tiles (H / rem-notH / y-made / COMMON -- common would break 'wake or y-made').
    Lesion: alpha +- e_i must break CT1.  Shared track (probe7): run w = translate of a block of x's journey.
(2) Control: CT1 with the 'BS alpha' on B/B non-commuting pairs (expected to fail; record at which step)."""
import sys, os
sys.path.insert(0, '/home/clio/projects/proofs/code-1118c2-Q597')
SZ = sys.argv[1]
src = open('/home/clio/projects/proofs/code-1118c2-Q597/bs6.py').read().split('def work(w):')[0]
exec(src)          # sets SIZE, NW, NCH, census, states, mvl, middles, bs_test, disp, tr, cen, rh ...
from stepset import candidates, step as sstep
import cmath
DIRS = [cmath.exp(1j*cmath.pi*k/12) for k in range(24)]
def T(S, v): return frozenset((k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S)
def bare(b): return {(v[0], frozenset(v[1])) for v in b.tiles.values()}
def apply(Bd, steps):
    Bd = set(Bd)
    for _, g, n in steps:
        assert g <= Bd, 'step not applicable'
        Bd -= g; Bd |= n
    return Bd
def tvec(s1, s2):
    if s1[0] != s2[0]: return None
    v = tuple(b-a for a, b in zip(cen(rh(s1[1])), cen(rh(s2[1]))))
    return v if tr(s1, v) == s2 else None
def rhof(s): return next(u for u in s[1] if u[0] == 'rh')
def finddir(MP, mA):
    for d in DIRS:
        Bm = set(MP); good = True
        for s in mA:
            got, _ = sstep(Bm, rhof(s), d)
            if got is None or tuple(got) != tuple(s): good = False; break
            Bm = apply(Bm, [s])
        if good: return d
    return None
def ct(MP, MxP, mA, mB, d, alpha, H=None, rem=None, kinds=None):
    """CT1 per C-step; returns index of first failing step or -1"""
    a, C = mA[0], mA[1:]; ma = tuple(-v for v in alpha)
    Bm = apply(MP, [a]); Bx = set(MxP)
    for i, (s, t) in enumerate(zip(C, mB[:-1])):
        cs = candidates(Bm, rhof(s), d)
        if kinds is not None:
            for _, S, _ in cs:
                for u in S:
                    if u[0] == 'rh': continue
                    kinds['H' if u in H else ('rem-notH' if u in rem else ('y-made' if u not in MP else 'COMMON'))] += 1
        if {T(S, ma) for _, S, _ in cs} != {S for _, S, _ in candidates(Bx, rhof(t), d)}: return i
        Bm = apply(Bm, [s]); Bx = apply(Bx, [t])
    return -1
def work(w):
    st = Counter(); ex = []; nstates = 0
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
                        P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                        X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                        com = tiles_of(b12) == tiles_of(b21)
                        if kind == 'BB':
                            if com: continue
                            # control: treat y as 'clean' with x as other; CT1 with alpha from y's first middle step
                            for nm, A, B, Jo in (('y', Jy, Jyx, Jx), ('x', Jx, Jxy, Jy)):
                                mA, mB, p, q = middles(A, B)
                                if len(mA) < 2 or len(mB) < 1: st[('BB ctrl', nm, 'middle too short', len(mA), len(mB))] += 1; continue
                                MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
                                d = finddir(MP, mA)
                                if d is None: st[('BB ctrl', nm, 'no dir')] += 1; continue
                                alpha = disp(mA[0], mA[0])
                                f = ct(MP, MxP, mA, mB, d, alpha)
                                lenok = len(mA) == len(mB)
                                st[('BB ctrl', nm, 'CT1 holds' if f < 0 else 'CT1 fails', 'len equal' if lenok else 'len differ')] += 1
                            continue
                        K = set(X) | set(Y) | set(Yx) | set(Xy)
                        I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                        sides = {k[0] for k in I}
                        if not I:
                            if not (set(X) & set(Y)): continue
                            for nm, A, B in (('y', Jy, Jyx), ('x', Jx, Jxy)):
                                mA, mB, p, q = middles(A, B)
                                if bs_test(mA, mB): st[('I0', nm, 'BS True', len(mA), len(mB))] += 1; continue
                                if A == B: why = 'journeys identical (no middle)'
                                elif len(mA) < 2: why = f'mA len {len(mA)} < 2'
                                elif len(mA) != len(mB): why = 'len mismatch'
                                else: why = 'translation mismatch'
                                st[('I0', nm, 'BS False', why, len(mA), len(mB))] += 1
                                if len(ex) < 30 and why == 'translation mismatch': ex.append((key, x, y, nm, mA, mB))
                            continue
                        if len(sides) != 1: st['two-sided'] += 1; continue
                        side = sides.pop()
                        A, B, Jo = (Jy, Jyx, Jx) if side == 0 else (Jx, Jxy, Jy)
                        mA, mB, p, q = middles(A, B)
                        if not bs_test(mA, mB): st[('BS fail', side)] += 1; continue
                        MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
                        rem, add = MP - MxP, MxP - MP
                        a, C = mA[0], mA[1:]
                        alpha = disp(a, a); ma = tuple(-v for v in alpha)
                        gamma = disp(C[0], C[-1])
                        H = {t for t in rem if T([t], ma) <= add}
                        d = finddir(MP, mA)
                        st[(side, 'direction found', d is not None)] += 1
                        if d is None: continue
                        kinds = Counter()
                        f = ct(MP, MxP, mA, mB, d, alpha, H, rem, kinds)
                        st[(side, 'CT1 (C cands translate)', f < 0)] += 1
                        if f >= 0 and len(ex) < 30: ex.append(('CT1fail', key, x, y, side, f))
                        st[(side, 'COMMON candidate tiles', kinds['COMMON'] > 0)] += 1
                        st[(side, 'kinds', tuple(sorted(kinds)))] += 1
                        ok2 = {T(S, gamma) for _, S, _ in candidates(MP, rhof(a), d)} == {S for _, S, _ in candidates(MxP, rhof(mB[-1]), d)}
                        st[(side, 'CT2 (last = a + gamma)', ok2)] += 1
                        les = sum(ct(MP, MxP, mA, mB, d, tuple(u+e for u, e in zip(alpha, ee))) < 0 for ee in E)
                        st[(side, 'CT1 lesion passes (must be 0)', les)] += 1
                        st[(side, 'len C', len(C))] += 1
                        # shared track (probe7)
                        wrun = mA[2:-1]
                        if not wrun: st[(side, 'run empty')] += 1
                        else:
                            hit = None
                            for i in range(len(Jo)-len(wrun)+1):
                                v = tvec(Jo[i], wrun[0])
                                if v is not None and all(tr(Jo[i+j], v) == wrun[j] for j in range(len(wrun))): hit = (i, v); break
                            st[(side, 'run = translate of x-block', hit is not None)] += 1
                            if hit: st[(side, 'delta', tuple(str(c) for c in hit[1]))] += 1
                            st[(side, 'middle shape a.b.w.b', mA[0][0] == 'a' and mA[1][0] == 'b' and mA[-1][0] == 'b' and all(s[0] == 'a' for s in wrun))] += 1
    return st, ex, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ex = []; ns = 0
    for s, e, n in res: st.update(s); ex += e; ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
    pickle.dump(dict(st=st, ex=ex), open(f'ct{SIZE}.pkl', 'wb'))
