"""PROVE 11-22 pairdich.py (from ccdirty.py): pair-level dichotomy over ALL non-commuting B/C pairs (any I), both roles.
Per pair: (side of I: 0 / 1 / two-sided / empty), and per role (x->B, y->C): diverges?, Sheet_p, |S|, local type.
Original: PROVE 11-22 ccdirty.py (from 11-21 c2 handshake.py): the contact configuration on the DIRTY role of B/C pairs
(the other journey as A; its re-run on the board after the clean journey), versus the clean role. Every role is
kept if it has a divergence; middles/bs_test are reported, not filtered.
Original: PROVE 11-21 c2: handshake geometry (Sh1) + contact configuration (CC) offsets and the Sheet law (for size 5 regression and the size-6 subsample).
Per record (clean B/C) and per role (BBy, BBx controls), A = journey of the role's clean corner on M0, Jo = the other.
 SHEET   every W tile removed by A was removed by a Jo step that is NOT a class-1 a-step;
 SH1     contact step p and step p+1 each read exactly the Jo steps {c0, c0+1}; A's step p+1 is (1,b);
 SH2     for i = 2..m (m = e - c0, e = last Jo step read): R_{p+i} = R^o_{c0+i} and A[p+i], Jo[c0+i] remove the same tiles;
 BC      (hypothesis of Lemma lock-step) every tile of Jo's hexagon at c0+i is on A's board before step p+i, i = 2..m;
 SH3     Jo[e] is (0,b), A[p+m] is (0,b) with the same removed tiles, and A reads no W tile after step p+m;
 TAIL    number of A steps after p+m, and Jo steps after e (is the post-split window non-trivial?).
Preconditions: class == #b parity (A and Jo); fc == p.   usage: KEEP=k python3 shadowlaw.py SIZE NW NCH"""
import os
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def parity(J, j): return sum(1 for s in J[:j] if s[0] == 'b') % 2
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
                            st[('BC final boards equal (commute as boards)', tiles_of(b12) == tiles_of(b21))] += 1
                            sides = {k[0] for k in I}
                            sideI = 'empty' if not sides else ('two-sided' if len(sides) > 1 else ('side0 (clean->C)' if 0 in sides else 'side1 (clean->B)'))
                            roles = (('x->B', Jx, Jxy, Jy), ('y->C', Jy, Jyx, Jx))
                            PAIR = {}
                        for tag, A, B, Jo in roles:
                            mA, mB, p, q = middles(A, B)
                            if kind == 'BB': continue
                            bsok = bs_test(mA, mB)
                            if tag.startswith('clean') and not bsok: continue
                            if len(mA) < 1 or len(mB) < 1: st[(tag, 'no divergence (re-run identical / trivial)')] += 1; PAIR[tag] = 'nodiv'; continue
                            st[(tag, 'BS', bsok)] += 1
                            PAIR[tag] = None
                            st[(tag, 'records')] += 1
                            TO = typed(Jo); TA = typed(A)
                            st[(tag, 'PRE class==#b parity (A and Jo)', all(TO[j][0] == parity(Jo, j) for j in range(len(Jo))) and all(TA[j][0] == parity(A, j) for j in range(len(A))))] += 1
                            W = set(M0) - set(apply(M0, Jo))
                            rem = {}
                            for j, s in enumerate(Jo):
                                for u in s[1]: rem.setdefault(u, j)
                            reads = [sorted({rem[u] for u in s[1] if u in W}) for s in A]
                            fc = next((k for k, r in enumerate(reads) if r), None)
                            st[(tag, 'PRE fc == p', fc == p)] += 1
                            Hp_ = {u for u in A[p][1] if u[0] != 'rh'}
                            sheetp = not any(TO[rem[u]][0] == 1 and TO[rem[u]][1] in ('a0', 'a1') for u in Hp_ if u in W)
                            PAIR[tag] = (sheetp, len(Hp_ & W))
                            st[(tag, 'LOCAL TYPE: A', TA[p], 're-run', typed(B)[p], 'J reads', tuple(TO[j] for j in sorted({rem[u] for u in Hp_ if u in W})), '|S|', len(Hp_ & W), 'Sheet_p', sheetp)] += 1
                        if kind == 'BC':
                            st[('PAIR', sideI, 'x->B', PAIR.get('x->B'), 'y->C', PAIR.get('y->C'))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns, 'KEEP', os.environ.get('KEEP', NCH), 'of', NCH)
    for k in sorted(st, key=str): print(st[k], k)
