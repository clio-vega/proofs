"""PROVE 11-22 retreat (a2read.py; a2mid.py adds the middle-end test): is every (A)-tile a READ wake tile of the clean journey?
(A)-tiles exactly as in code-1121-Q597/bbsplit.py / g1aparity.py: non-rh tiles u of any CANDIDATE at the rhombus of
mA[i+1] on M_i = M_P . mA[:i+1] (i >= 0), not made by the middle, with -alpha an edge vector (A-shaped), and u - alpha in M0.
Per A-tile: u in W?  u removed by the clean journey A on M0 (a READ wake tile)? at which step relative to p?
J_o's removing step: (class, type), b-parity before it ((A2) 'even'), and (Sheet)-class (sheet_pre, sheet_post).
Also the brief's sub-question: is the H_p square (removed at p) ever an A-tile?  (expected never: removed at p, Lemma nr.)
Control: B/B roles.   usage: python3 a2read.py SIZE NW NCH"""
import os
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def parity(J, j): return sum(1 for s in J[:j] if s[0] == 'b') % 2
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
                    for kind in ('BC', 'BB'):
                        if kind == 'BB' and not x < y: continue
                        tx, ty = kind
                        b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                        b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                        if kind == 'BB':
                            if tiles_of(b12) == tiles_of(b21): continue
                            roles = (('BBy', Jy, Jyx, Jx), ('BBx', Jx, Jxy, Jy))
                        else:
                            P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                            X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                            K = set(X) | set(Y) | set(Yx) | set(Xy)
                            I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                            sides = {k[0] for k in I}
                            if len(sides) != 1: continue
                            side = sides.pop()
                            roles = (('clean->C', Jy, Jyx, Jx),) if side == 0 else (('clean->B', Jx, Jxy, Jy),)
                        for tag, A, B, Jo in roles:
                            mA, mB, p, q = middles(A, B)
                            if kind == 'BC' and not bs_test(mA, mB): continue
                            if len(mA) < 2 or len(mB) < 1: continue
                            st[(tag, 'records')] += 1
                            MP = apply(M0, A[:p]); d = finddir(MP, mA)
                            if d is None: st[(tag, 'finddir None')] += 1; continue
                            al = disp(mA[0], mA[0]); mal = tuple(-v for v in al)
                            W = set(M0) - set(apply(M0, Jo))
                            TO = typed(Jo)
                            remJ = {}
                            for j, s in enumerate(Jo):
                                for u in s[1]: remJ.setdefault(u, j)
                            remA = {}
                            for k, s in enumerate(A):
                                for u in s[1]: remA.setdefault(u, k)
                            sqp = [u for u in A[p][1] if u[0] == 'sq'][0]
                            seen = set(); nA = 0
                            for i in range(len(mA) - 1):
                                Mi, whoM = made_by(MP, list(mA[:i+1]))
                                for _, S, _ in candidates(Mi, rhof(mA[i+1]), d):
                                    for u in S:
                                        if u[0] == 'rh' or u in whoM or u in seen: continue
                                        seen.add(u)
                                        if not edgevec(u, mal): continue
                                        u1 = next(iter(T([u], mal)))
                                        if u1 not in M0: continue
                                        nA += 1
                                        if u == sqp: st[(tag, 'H_p square is an A-tile (expect never)')] += 1
                                        jr = remJ.get(u); ka = remA.get(u)
                                        cl = None if jr is None else (TO[jr], 'b-parity before', parity(Jo, jr))
                                        st[(tag, 'A-tile', 'in W', u in W, 'read by A', ka is not None,
                                            'A step - p', None if ka is None else ka - p, 'Jo removing step', cl)] += 1
                            st[(tag, 'records with >=1 A-tile', nA > 0)] += 1
                            reads = [sorted({remJ[u] for u in s_[1] if u in W}) for s_ in A]
                            c0 = reads[p][0]; rd = [k for k, r in enumerate(reads) if r]
                            e = max(j for r in reads for j in r); m = e - c0
                            st[(tag, 'MIDDLE: last middle step - p', len(mA) - 1, 'split offset m = e - c0', m, 'last read step - p', rd[-1] - p)] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns, 'KEEP', os.environ.get('KEEP', NCH), 'of', NCH)
    for k in sorted(st, key=str): print(st[k], k)
