"""PROVE 11-20 c2 Step 0c (4-offsets om of u-alpha vs the removed-offset set Q; from threat.py) Step 0b (threat sets K={k: psi(u-alpha)<=mu_k}) derived from psisep.py.\nPROVE 11-20 c2 Step 0 (G1c): does the two-height potential psi separate every (A)-tile neighbour u-alpha from
x's journey J_o?  Test: NOT (psi(u-alpha) <= mu_0 componentwise), i.e. exists i with psi_i(u-alpha) > (mu_0)_i >= (mu_k)_i.
By REM (every removed non-rh tile has psi <= mu_k) and MU this forbids removal at every step k.
Preconditions asserted per A-tile: u-alpha in M0; at the step j of J_o removing u, -alpha is an edge vector of R_j.
Control: B/B pairs (WS fails there), same recipe on read tiles with -alpha an edge vector of u.
Lesion: psi computed with eta_B and eta_C swapped.   usage: python3 psisep.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
SWAP = {'B': 'C', 'C': 'B'}
def psi(t, TT):
    eta = ETA[TT]
    return (min(-eta(P[2:]) for P in t[1]), min(eta(P[:2]) for P in t[1]))
def le(a, b): return a[0] <= b[0] and a[1] <= b[1]
def mus(J, TT): return [psi(rhof(s), TT) for s in J]
def sep(t, J, TT):
    """(exists i forall k psi_i(t) > mu_k_i,  t removed by some step of J,  first k with psi(t)<=mu_k or None)"""
    M = mus(J, TT); pt = psi(t, TT)
    ex = any(all(pt[i] > m[i] for m in M) for i in (0, 1))
    rem = any(t in g for typ, g, n in J)
    fk = next((k for k, m in enumerate(M) if le(pt, m)), None)
    return ex, rem, fk
from offsets import om
QS = {'B': {(-1,0,0,0),(0,-1,0,0),(0,0,0,0),(0,0,1,0),(-1,0,1,0)}, 'C': {(-1,0,0,0),(0,-1,0,0),(0,0,0,0),(0,0,0,1),(0,-1,0,1)}}
def edgevec(t, v): return mult_of_edge(t, v) is not None
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
                    if x < y:   # B/B control
                        c1, Kx = mvl(b, x, 'B', DV); c2, Ky = mvl(b, y, 'B', DV)
                        c12, Kyx = mvl(c1, y, 'B', DV); c21, Kxy = mvl(c2, x, 'B', DV)
                        if tiles_of(c12) != tiles_of(c21):
                            for nm, A_, B_, Jo_ in (('y', Ky, Kyx, Kx), ('x', Kx, Kxy, Ky)):
                                mA_, mB_, p_, q_ = middles(A_, B_)
                                if len(mA_) < 2: continue
                                MP_ = apply(M0, A_[:p_])
                                d_ = finddir(MP_, mA_)
                                if d_ is None: continue
                                al = disp(mA_[0], mA_[0]); mal = tuple(-v for v in al)
                                seen = set()
                                for i in range(len(mA_) - 1):
                                    Mi, whoM = made_by(MP_, list(mA_[:i+1]))
                                    for _, S, _ in candidates(Mi, rhof(mA_[i+1]), d_):
                                        for u in S:
                                            if u[0] == 'rh' or u in whoM or u in seen: continue
                                            seen.add(u)
                                            if not edgevec(u, mal): continue
                                            u1 = next(iter(T([u], mal)))
                                            if u1 not in M0: st[('BB', 'A-shaped, u-alpha not in M0')] += 1; continue
                                            ex, rem, fk = sep(u1, Jo_, 'B')
                                            st[('BB', 'A-shaped tile', u[0], 'psi-sep', ex, 'removed by Jo', rem)] += 1
                                            st[('BB', 'A-shaped tile', 'some k om in Q', any(om(u1, rhof(s), 'B') in QS['B'] for s in Jo_), 'removed by Jo', rem)] += 1
                    b1, Jx = mvl(b, x, 'B', DV); b2, Jy = mvl(b, y, 'C', DV)
                    b12, Jyx = mvl(b1, y, 'C', DV); b21, Jxy = mvl(b2, x, 'B', DV)
                    P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                    X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                    K = set(X) | set(Y) | set(Yx) | set(Xy)
                    I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                    sides = {k[0] for k in I}
                    if len(sides) != 1: continue
                    side = sides.pop()
                    A, B, Jo = (Jy, Jyx, Jx) if side == 0 else (Jx, Jxy, Jy)
                    To = 'B' if side == 0 else 'C'; Tp = SWAP[To]
                    mA, mB, p, q = middles(A, B)
                    if not bs_test(mA, mB): continue
                    st['records'] += 1
                    MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
                    a, C = mA[0], mA[1:]
                    alpha = disp(a, a); ma = tuple(-v for v in alpha)
                    d = finddir(MP, mA)
                    if d is None: st['d None'] += 1; continue
                    seen = set()
                    for i in range(len(C)):
                        Mi, whoM = made_by(MP, [a] + list(C[:i]))
                        for _, S, _ in candidates(Mi, rhof(C[i]), d):
                            for u in S:
                                if u[0] == 'rh' or u in whoM or u in MxP or u in seen: continue
                                seen.add(u)
                                kind = 'A' if edgevec(u, ma) else 'BT'
                                u1 = next(iter(T([u], ma)))
                                if kind == 'BT':
                                    # control-ish: B/T neighbours are made by x; psi-sep at k=0 is not expected
                                    ex, rem, fk = sep(u1, Jo, To) if u1 in M0 else (None, None, None)
                                    st[(side, 'BT', u[0], 'u-alpha in M0', u1 in M0)] += 1
                                    continue
                                # preconditions
                                st[(side, 'A', 'PRE u-alpha in M0', u1 in M0)] += 1
                                jr = next((j for j, s in enumerate(Jo) if u in s[1]), None)
                                pre2 = jr is not None and edgevec(rhof(Jo[jr]), ma)
                                st[(side, 'A', 'PRE removed by Jo at j, -alpha edge vector of R_j', pre2, 'type', Jo[jr][0] if jr is not None else None)] += 1
                                if u1 not in M0: continue
                                tk = u[0] if u[0] != 'tri' else ('tri0' if len({P[2:] for P in u[1]}) == 1 else 'tri1')
                                lo = jr
                                while lo > 0 and Jo[lo-1][0] == 'a' and Jo[jr][0] == 'a': lo -= 1
                                Om = [om(u1, rhof(s), To) for s in Jo]
                                inQ = [k for k, o in enumerate(Om) if o in QS[To]]
                                st[(side, tk, 'rm', Jo[jr][0], 'om(u,R_jr)', om(u, rhof(Jo[jr]), To), 'om(u-a,R_jr)', Om[jr], 'om(u-a,R_0)', Om[0])] += 1
                                st[(side, tk, 'rm', Jo[jr][0], 'steps k with om(u-a,R_k) in Q', 'none' if not inQ else ('within run' if min(inQ) >= lo and max(inQ) <= jr else tuple(k - jr for k in inQ)))] += 1
                                st[('ALL A', 'some k has om(u-a,R_k) in Q', bool(inQ))] += 1
                                # B/T control: u-alpha made by x; is it in Q at the step it is made? (should be: it is near R)
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
