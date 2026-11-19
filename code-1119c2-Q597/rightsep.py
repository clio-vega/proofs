"""PROVE 11-19 c2: is global survival explained by the engine's rightness filter?
For every distinct B/C journey (<=5), migration direction d (finddir on the full journey from the start board),
every non-rh tile t made at step j and every later step j' (j' > j):
  (S1) max x_d(t) <= min x_d(H_j')   [separation by the migration coordinate]
  (S2) max x_d(t) <= xmin(R_j')      [the tile lies left of the rhombus's leftmost vertex]
  (V)  t removed by step j' (violation of survival)."""
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
def work(w):
    st = Counter(); ns = 0; seenJ = set()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1; M0 = bare(b)
            for x in corners(set(b.nest[0]['pos'])):
                for TT in ('B', 'C'):
                    _, J = mvl(b, x, TT, DV)
                    if not J or tuple(J) in seenJ: continue
                    seenJ.add(tuple(J))
                    d = finddir(M0, J)
                    if d is None: st['no dir'] += 1; continue
                    rot = d.conjugate(); X = lambda P: (G.phi(P)*rot).real
                    st['journeys'] += 1
                    for j, (_, g, n) in enumerate(J):
                        for t in n:
                            if t[0] == 'rh': continue
                            mx = max(X(P) for P in t[1])
                            for jj in range(j+1, len(J)):
                                _, g2, _ = J[jj]
                                Hmin = min(X(P) for u in g2 for P in u[1])
                                Rmin = min(X(P) for P in rhof(J[jj])[1])
                                st[('S1', mx <= Hmin + 1e-9)] += 1
                                st[('S2', mx <= Rmin + 1e-9)] += 1
                                if t in g2: st['V removed later'] += 1
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
