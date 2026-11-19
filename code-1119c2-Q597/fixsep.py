"""PROVE 11-19 c2 (fixsep): does ONE linear functional per journey separate (scan 720 directions)? Also per-step functional (per j') ? Originally: is global survival explained by the engine's rightness filter?
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
                    import cmath, math
                    TH = [cmath.exp(1j*math.pi*k/360) for k in range(720)]
                    pts = lambda S: [G.phi(P) for u in S for P in u[1]]
                    ok = [True]*720; okstep = []
                    for jj in range(1, len(J)):
                        Hp = pts(J[jj][1])
                        prev = [t for j in range(jj) for t in J[j][2] if t[0] != 'rh' and t not in J[j+1:jj+1] ]
                        made = set(t for j in range(jj) for t in J[j][2] if t[0] != 'rh')
                        Tp = [G.phi(P) for t in made for P in t[1]]
                        if not Tp: continue
                        stepok = False
                        for k, z in enumerate(TH):
                            f = lambda w: (w*z.conjugate()).real
                            good = max(map(f, Tp)) <= min(map(f, Hp)) + 1e-9
                            if not good: ok[k] = False
                            else: stepok = True
                        okstep.append(stepok)
                    st['journeys'] += 1
                    st[('one functional separates all made tiles from all later hexagons', any(ok))] += 1
                    st[('per step: some functional separates made-so-far from H_jj (all steps)', all(okstep))] += 1
                    st[('journey shape', len(J) <= 3)] += 0
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
