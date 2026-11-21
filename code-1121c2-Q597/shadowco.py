"""PROVE 11-21 c2 Step 0d: rhombus coincidence in the shadow.  For every step k of the clean journey A on M0 that removes a tile u of
W = M0 minus M0.Jo, record J_o's sheet at the step j that removed u:  sheet_pre = #b in Jo[:j] mod 2 (= class of R^o_j),
sheet_post = #b in Jo[:j+1] mod 2 (class after the step).  Split: contact (k = fc) / after contact, kind, Jo (class,type).
Preconditions: record counts; typed class == word parity for Jo and A (Lemma word); Jo removes each W tile exactly once.
Lesions: #a1 in Jo[:j] mod 2, j mod 2 (with a differs-from-sheet count).  Controls: BBy, BBx.
usage: python3 sheettab.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def parity(J, j, typ): return sum(1 for s in J[:j] if s[0] == typ) % 2

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
                            TO = typed(Jo); TA = typed(A)
                            Mx0 = apply(M0, Jo); W = set(M0) - set(Mx0)
                            rem_ = {}
                            for j, s in enumerate(Jo):
                                for u in s[1]: rem_.setdefault(u, j)
                            readby = {}
                            for k, s in enumerate(A):
                                for u in s[1]:
                                    if u in W: readby.setdefault(u, k)
                            pairs = []
                            for k, s in enumerate(A):
                                js = sorted({rem_[u] for u in s[1] if u in W})
                                pairs.append(js)
                            allj = sorted({j for js in pairs for j in js})
                            e = max(allj)
                            seq = ' '.join(('%d%s' % TA[k]) + ('[' + ','.join(str(j - e) for j in js) + ']' if js else '') + ('|' if k == p else '') for k, js in enumerate(pairs))

                            c0 = pairs[p][0] if pairs[p] else None
                            m = e - c0
                            # rhombus coincidence: R_{p+1+i} vs R^o_{c0+1+i}, i = 0..m; also hexagon (removed set) equality A[p+1+i] vs Jo[c0+1+i]
                            RA = [rhof(s) for s in A]; RO = [rhof(s) for s in Jo]
                            co = []
                            for i in range(0, m + 2):
                                ka, jo = p + i, c0 + i
                                if ka < len(A) and jo < len(Jo): co.append('R' if RA[ka] == RO[jo] else ('h' if set(A[ka][1]) == set(Jo[jo][1]) else '.'))
                            st[(tag, 'COINC R_{p+i}==R^o_{c0+i}, i=0..m+1', ''.join(co))] += 1
                            # same removed set (hexagon) for shadow steps p+2..p+m vs Jo c0+2..e
                            hx = [set(A[p+1+i][1]) == set(Jo[c0+1+i][1]) for i in range(1, m) if p+1+i < len(A)]
                            st[(tag, 'shadow steps remove identical tile sets', all(hx), 'n', len(hx))] += 1
                            fin = Jo[e]; fa = A[p + m] if p + m < len(A) else None
                            st[(tag, 'split step: shared tiles of A[p+m] and Jo[e]', len(set(fin[1]) & set(fa[1])) if fa else None, 'R equal before split', RA[p+m] == RO[e] if fa else None)] += 1

                            # order law: are J_o step indices read in weakly decreasing order along A?
                            flat = [j for js in pairs for j in js]
                            st[(tag, 'reads in J_o order: nonincreasing', all(flat[i] >= flat[i+1] for i in range(len(flat)-1)), 'nondecreasing', all(flat[i] <= flat[i+1] for i in range(len(flat)-1)))] += 1
                            st[(tag, 'x contact reads Jo steps (rel. to e)', tuple(j - e for j in pairs[p]))] += 1
                            st[(tag, 'Jo has class-1 a-step', any(c == 1 and t != 'b' for c, t in TO))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str):
        if st[k]: print(st[k], k)
