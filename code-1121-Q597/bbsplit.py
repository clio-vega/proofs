"""PROVE 11-21 control: B/B A-shaped read tiles (same recipe as g1aparity.py), split by ROLE of the clean journey
(BBy: second corner, BBx: first corner), by (A1) [-alpha edge vector of Jo's R_0], (A2) [even #b of Jo before removal],
and whether u-alpha is removed by Jo (the WS-b failure).  Lesion for (A1): R_0 replaced by the clean journey's R_a.
usage: python3 bbsplit.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
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
                    if not x < y: continue
                    c1, Kx = mvl(b, x, 'B', DV); c2, Ky = mvl(b, y, 'B', DV)
                    c12, Kyx = mvl(c1, y, 'B', DV); c21, Kxy = mvl(c2, x, 'B', DV)
                    if tiles_of(c12) == tiles_of(c21): continue
                    for nm, A_, B_, Jo_ in (('BBy', Ky, Kyx, Kx), ('BBx', Kx, Kxy, Ky)):
                        mA_, mB_, p_, q_ = middles(A_, B_)
                        if len(mA_) < 2: continue
                        MP_ = apply(M0, A_[:p_])
                        d_ = finddir(MP_, mA_)
                        if d_ is None: continue
                        al = disp(mA_[0], mA_[0]); mal = tuple(-v for v in al)
                        A1 = edgevec(rhof(Jo_[0]), mal); A1les = edgevec(rhof(mA_[0]), mal)
                        st[(nm, 'journeys', 'A1', A1, 'lesion A1 (R_a)', A1les)] += 1
                        seen = set()
                        for i in range(len(mA_) - 1):
                            Mi, whoM = made_by(MP_, list(mA_[:i+1]))
                            for _, S, _ in candidates(Mi, rhof(mA_[i+1]), d_):
                                for u in S:
                                    if u[0] == 'rh' or u in whoM or u in seen: continue
                                    seen.add(u)
                                    if not edgevec(u, mal): continue
                                    u1 = next(iter(T([u], mal)))
                                    if u1 not in M0: st[(nm, 'A-shaped, u-alpha not in M0')] += 1; continue
                                    jr = next((j for j, s in enumerate(Jo_) if u in s[1]), None)
                                    A2 = None if jr is None else sum(1 for s in Jo_[:jr] if s[0] == 'b') % 2 == 0
                                    rem = any(u1 in s[1] for s in Jo_)
                                    st[(nm, 'A-shaped tile', 'A1', A1, 'u removed by Jo', jr is not None, 'A2 even', A2, 'u-alpha removed by Jo', rem)] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
