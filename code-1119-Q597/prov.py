"""PROVE 11-19 Step 2 instrument: provenance dictionary for CT1/CT2.  usage: python3 prov.py {5|6} NW NCH
For each clean nontrivial B/C record, each crossing step i (state M_i, rhombus R_i; partner Mx_i, R_i - alpha),
each candidate S on either side, each non-rh tile u: provenance of u on its side and of its translate on the other side.
Provenance in M_i:  'a' (made by y's step a), 'c<j>' made by c_j (j<i) -> 'c', 'wake' (M_P \\ Mx_P), 'common' (M_P & Mx_P).
Provenance in Mx_i: 'c' (made by c_j - alpha), 'xadd' (Mx_P \\ M_P), 'common'.  Also: is the tile x-displaced in a type-a step (zone e0/e1)?"""
import sys
exec(open('/home/clio/projects/proofs/code-1119-Q597/ct6.py').read().split('\ndef work(w):')[0])
def made_by(Bstart, steps):
    """map tile -> index of step that created it (last creation), for tiles present at the end"""
    Bd = set(Bstart); who = {}
    for j, (_, g, n) in enumerate(steps):
        Bd -= g; Bd |= n
        for t in n: who[t] = j
        for t in g: who.pop(t, None)
    return Bd, who
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
                    if x < y:   # B/B control: same recipe, a = first middle step, C = rest
                        c1, Kx = mvl(b, x, 'B', DV); c2, Ky = mvl(b, y, 'B', DV)
                        c12, Kyx = mvl(c1, y, 'B', DV); c21, Kxy = mvl(c2, x, 'B', DV)
                        if tiles_of(c12) != tiles_of(c21):
                            for nm, A_, B_, Jo_ in (('y', Ky, Kyx, Kx), ('x', Kx, Kxy, Ky)):
                                mA_, mB_, p_, q_ = middles(A_, B_)
                                if len(mA_) < 2: st[('BB', nm, 'middle<2')] += 1; continue
                                MP_ = apply(M0, A_[:p_]); MxP_ = apply(apply(M0, Jo_), B_[:p_])
                                d_ = finddir(MP_, mA_)
                                if d_ is None: st[('BB', nm, 'nodir')] += 1; continue
                                al = disp(mA_[0], mA_[0]); mal = tuple(-v for v in al)
                                first = None
                                for i in range(min(len(mA_)-1, len(mB_))):
                                    Mi, whoM = made_by(MP_, list(mA_[:i+1])); Xi, whoX = made_by(MxP_, list(mB_[:i]))
                                    for _, S, _ in candidates(Mi, rhof(mA_[i+1]), d_):
                                        for u in S:
                                            if u[0] == 'rh': continue
                                            u2 = next(iter(T([u], mal)))
                                            pm = ('a' if whoM[u] == 0 else 'c') if u in whoM else ('wake' if u not in MxP_ else 'common')
                                            if u2 not in Xi and first is None: first = (i, pm)
                                    if first is None:
                                        for _, S, _ in candidates(Xi, rhof(mB_[i]), d_):
                                            for u in S:
                                                if u[0] == 'rh': continue
                                                u2 = next(iter(T([u], al)))
                                                if u2 not in Mi and first is None: first = (i, 'X-side')
                                    if first is not None: break
                                if first is None: st[('BB', nm, 'M+X hold on compared steps; lengths', len(mA_), len(mB_))] += 1
                                st[('BB', nm, 'dictionary', 'holds' if first is None else ('fails at C-step', first[0], 'tile prov', first[1]))] += 1
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
                    mA, mB, p, q = middles(A, B)
                    if not bs_test(mA, mB): st['BS fail'] += 1; continue
                    st['records'] += 1
                    MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
                    # x's journey Jo runs from M0; tiles it displaced in type-a steps, by zone vector
                    disp_a = {}
                    for s in Jo:
                        if s[0] != 'a': continue
                        for t in s[1]:
                            if t[0] == 'rh': continue
                            for u in s[2]:
                                if u[0] == t[0]:
                                    v = tuple(bb-aa for aa, bb in zip(cen(t[1]), cen(u[1])))
                                    if T([t], v) == frozenset([u]): disp_a[t] = v
                    a, C = mA[0], mA[1:]
                    alpha = disp(a, a); ma = tuple(-v for v in alpha)
                    d = finddir(MP, mA)
                    if d is None: st['no dir'] += 1; continue
                    def provM(u, who):
                        if u in who: return 'a' if who[u] == 0 else 'c'
                        return 'wake' if u not in MxP else 'common'
                    def provX(u, who):
                        if u in who: return 'c'
                        return 'xadd' if u not in MP else 'common'
                    for i in range(len(C) + 1):
                        Mi, whoM = made_by(MP, [a] + list(C[:i]))
                        Xi, whoX = made_by(MxP, list(mB[:i]))
                        if i < len(C): RM, RX, sh = rhof(C[i]), rhof(mB[i]), ma; lab = 'C'
                        else: RM, RX, sh = rhof(a), rhof(mB[-1]), disp(C[0], C[-1]); lab = 'last'
                        if lab == 'last':
                            Mi, whoM = MP, {}   # CT2 compares Cand(M_P, R_a) + gamma with Cand(Mx_end, R_a + gamma)
                        cM = candidates(Mi, RM, d); cX = candidates(Xi, RX, d)
                        nchosen = len(cM)
                        st[(side, lab, 'ncand', nchosen)] += 1
                        for S0, Bself, Bother, whoS, whoO, pS, pO, v, nm in (
                                (cM, Mi, Xi, whoM, whoX, provM, provX, sh, 'M->Mx'),
                                (cX, Xi, Mi, whoX, whoM, provX, provM, tuple(-z for z in sh), 'Mx->M')):
                            for _, S, _ in S0:
                                for u in S:
                                    if u[0] == 'rh': continue
                                    u2 = next(iter(T([u], v)))
                                    zone = 'x-a-disp' if u in disp_a else ''
                                    st[(side, lab, nm, pS(u, whoS), '->', pO(u2, whoO) if u2 in Bother else 'ABSENT', zone)] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: KEEP = int(os.environ.get('KEEP', NCH)); res = list(pool.imap_unordered(work, range(KEEP)))   # KEEP<NCH: fixed subsample of chunks 0..KEEP-1
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns, 'chunks', KEEP, 'of', NCH)
    for k in sorted(st, key=str): print(st[k], k)
