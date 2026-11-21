"""PROVE 11-21 c2: one consolidated test of the Shadow law and the Sheet law (for size 5 regression and the size-6 subsample).
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
                            sides = {k[0] for k in I}
                            if len(sides) != 1: continue
                            side = sides.pop()
                            roles = (('clean->C', Jy, Jyx, Jx),) if side == 0 else (('clean->B', Jx, Jxy, Jy),)
                        for tag, A, B, Jo in roles:
                            mA, mB, p, q = middles(A, B)
                            if kind == 'BC' and not bs_test(mA, mB): continue
                            if len(mA) < 2 or len(mB) < 1: continue
                            st[(tag, 'records')] += 1
                            TO = typed(Jo); TA = typed(A)
                            st[(tag, 'PRE class==#b parity (A and Jo)', all(TO[j][0] == parity(Jo, j) for j in range(len(Jo))) and all(TA[j][0] == parity(A, j) for j in range(len(A))))] += 1
                            W = set(M0) - set(apply(M0, Jo))
                            rem = {}
                            for j, s in enumerate(Jo):
                                for u in s[1]: rem.setdefault(u, j)
                            reads = [sorted({rem[u] for u in s[1] if u in W}) for s in A]
                            fc = next(k for k, r in enumerate(reads) if r)
                            st[(tag, 'PRE fc == p', fc == p)] += 1
                            sheet = all(not (TO[j][0] == 1 and TO[j][1] != 'b') for r in reads for j in r)
                            st[(tag, 'SHEET', sheet)] += 1
                            c0 = reads[p][0]; e = max(j for r in reads for j in r); m = e - c0
                            sh1 = reads[p] == [c0, c0 + 1] and p + 1 < len(A) and reads[p + 1] == [c0, c0 + 1] and TA[p + 1] == (1, 'b')
                            sh2 = sh1 and p + m < len(A) and all(rhof(A[p+i]) == rhof(Jo[c0+i]) and A[p+i][1] == Jo[c0+i][1] for i in range(2, m + 1))
                            bc = True
                            if sh1:
                                for i in range(2, m + 1):
                                    if p + i >= len(A): bc = False; break
                                    Bd = apply(M0, A[:p+i])
                                    if not set(Jo[c0+i][1]) <= Bd: bc = False
                            sh3 = sh2 and TO[e] == (0, 'b') and TA[p+m] == (0, 'b') and all(not r for r in reads[p+m+1:])
                            st[(tag, 'SH1', sh1, 'SH2', sh2, 'BC', bc, 'SH3', sh3)] += 1
                            st[(tag, 'TAIL A after split', len(A) - 1 - (p + m), 'Jo after e', len(Jo) - 1 - e)] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns, 'KEEP', os.environ.get('KEEP', NCH), 'of', NCH)
    for k in sorted(st, key=str): print(st[k], k)
