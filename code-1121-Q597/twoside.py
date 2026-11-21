"""PROVE 11-21: the two-sided contact.  At the divergence rhombus R_a: the M-step g_a meets W (Jo-removed), and the Mx-step
g'_a meets Xadd = Mx0 minus M0 (Jo-created).  Record: kinds, Jo step(s) (class,type) removing the W tiles and CREATING the
Xadd tiles, whether they are the same Jo step, and the offset base(R_a) - base(R^o_j).   usage: python3 twoside.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def base(vs): return tuple(min(Fr(v[i]) for v in vs) for i in range(4))
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
                            Mx0 = apply(M0, Jo); W = set(M0) - set(Mx0); Xadd = set(Mx0) - set(M0)
                            TO = typed(Jo)
                            rem, cre = {}, {}
                            for j, s in enumerate(Jo):
                                for u in s[1]: rem.setdefault(u, j)
                                for u in s[2]: cre[u] = j
                            wt = [(u[0], rem[u]) for u in A[p][1] if u in W]
                            xt = [(u[0], cre[u]) for u in B[p][1] if u in Xadd]
                            js = {j for _, j in wt} | {j for _, j in xt}
                            st[(tag, 'W tiles', tuple(sorted((k,) + TO[j] for k, j in wt)), 'Xadd tiles', tuple(sorted((k,) + TO[j] for k, j in xt)))] += 1
                            st[(tag, 'Jo steps involved', len(js), 'span', (max(js) - min(js)) if js else None)] += 1
                            st[(tag, 'Xadd tile created by the step removing a W tile', bool({j for _, j in xt} & {j for _, j in wt}))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
