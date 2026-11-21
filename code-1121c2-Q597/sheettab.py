"""PROVE 11-21 c2 Step 0: the sheet falsifier.  For every step k of the clean journey A on M0 that removes a tile u of
W = M0 minus M0.Jo, record J_o's sheet at the step j that removed u:  sheet_pre = #b in Jo[:j] mod 2 (= class of R^o_j),
sheet_post = #b in Jo[:j+1] mod 2 (class after the step).  Split: contact (k = fc) / after contact, kind, Jo (class,type).
Preconditions: record counts; typed class == word parity for Jo and A (Lemma word); Jo removes each W tile exactly once.
Lesions: #a1 in Jo[:j] mod 2, j mod 2 (with a differs-from-sheet count).  Controls: BBy, BBx.
usage: python3 sheettab.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def parity(J, j, typ): return sum(1 for s in J[:j] if s[0] == typ) % 2
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
                            st[(tag, 'PRE records')] += 1
                            TO = typed(Jo); TA = typed(A)
                            st[(tag, 'PRE class==#b parity, Jo', all(TO[j][0] == parity(Jo, j, 'b') for j in range(len(Jo))))] += 1
                            st[(tag, 'PRE class==#b parity, A', all(TA[j][0] == parity(A, j, 'b') for j in range(len(A))))] += 1
                            Mx0 = apply(M0, Jo); W = set(M0) - set(Mx0)
                            rem = {}
                            for j, s in enumerate(Jo):
                                for u in s[1]: rem.setdefault(u, []).append(j)
                            st[(tag, 'PRE every W tile removed exactly once by Jo', all(len(rem.get(u, [])) == 1 for u in W))] += 1
                            fc = next((k for k, s in enumerate(A) if any(u in W for u in s[1])), None)
                            st[(tag, 'PRE fc == p', fc == p)] += 1
                            nread = 0
                            for k, s in enumerate(A):
                                for u in s[1]:
                                    if u not in W: continue
                                    nread += 1
                                    j = rem[u][0]
                                    sp, sq = parity(Jo, j, 'b'), parity(Jo, j + 1, 'b')
                                    where = 'contact' if k == fc else ('after' if k > fc else 'BEFORE?!')
                                    st[(tag, 'TAB', where, 'sheet_pre', sp)] += 1
                                    st[(tag, 'TAB', where, 'sheet_post', sq)] += 1
                                    st[(tag, 'CFG', where, u[0], TO[j])] += 1
                                    st[(tag, 'REL', where, 'Jo sheet_pre', sp, 'A sheet at k', TA[k][0])] += 1
                                    la, lj = parity(Jo, j, 'a1'), j % 2
                                    st[(tag, 'LES #a1', where, la)] += 1
                                    st[(tag, 'LES j', where, lj)] += 1
                                    st[(tag, 'LES differs: #a1 vs sheet_pre', la != sp)] += 1
                                    st[(tag, 'LES differs: j vs sheet_pre', lj != sp)] += 1
                                    st[(tag, 'BOTH sheets', where, (sp, sq))] += 1
                            st[(tag, 'records with >=1 read W tile after contact', any(k > fc for k, s in enumerate(A) for u in s[1] if u in W))] += 1
                            st[(tag, 'read W tiles per record total', nread)] += 0
                            st[(tag, 'NREAD sum')] += nread
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str):
        if st[k]: print(st[k], k)
