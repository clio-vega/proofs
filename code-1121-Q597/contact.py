"""PROVE 11-21 Step 0c: is the divergence step p the FIRST CONTACT of the clean journey (in M) with Jo's removed set?
Jo-removed set W = M0 \\ (M0 . Jo) (tiles of the start board that Jo removes and does not re-create).
Per record: first k with g_k (non-rh) meeting W; compare with p.  Also: does any prefix step k<p have an
alternative good config (other type, same class) whose availability differs between M_k and Mx_k?
Also nest positions of the two corners (i,j) and the relative offset.   usage: python3 contact.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
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
                            roles = (('BBy', Jy, Jyx, Jx, y, x), ('BBx', Jx, Jxy, Jy, x, y))
                        else:
                            P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                            X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                            K = set(X) | set(Y) | set(Yx) | set(Xy)
                            I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                            sides = {k[0] for k in I}
                            if len(sides) != 1: continue
                            side = sides.pop()
                            roles = (((side, 'clean->C'), Jy, Jyx, Jx, y, x),) if side == 0 else (((side, 'clean->B'), Jx, Jxy, Jy, x, y),)
                        for tag, A, B, Jo, cc, oc in roles:
                            mA, mB, p, q = middles(A, B)
                            if kind == 'BC' and not bs_test(mA, mB): continue
                            if len(mA) < 2 or len(mB) < 1: continue
                            Mx0 = apply(M0, Jo)
                            W = set(M0) - set(Mx0); Xadd = set(Mx0) - set(M0)
                            fc = next((k for k, s in enumerate(A) if any(u in W for u in s[1] if u[0] != 'rh')), None)
                            st[(tag, 'first contact == p', fc == p, 'fc<p' if fc is not None and fc < p else '')] += 1
                            # created-by-Jo tiles adjacent: does any prefix step read (remove) an Jo-added tile? (impossible if prefix common, sanity)
                            st[(tag, 'offset clean-other (nest i,j)', tuple(a - b_ for a, b_ in zip(cc, oc)))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
