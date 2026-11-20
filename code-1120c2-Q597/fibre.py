"""PROVE 11-20 c2 retreat: forgetting FIBRE test (connection 2026-11-20 falsifier 1), and the Hecke/braid pattern count.
For every non-commuting B/B pair (x,y) at a census state s, with common closing label m (c12.m = c21.m):
  D = xy(s) symmetric-difference yx(s)  (tiles).
  Prediction (a): D within wake(x) ∪ wake(y)  [first-journey wakes s△x(s), s△y(s)]; (a') within all four journey wakes.
  (b) m's wake covers D:  union-form  D ⊆ (c12△c12.m) ∪ (c21△c21.m)  is a TAUTOLOGY (instrument check, must be 100%);
      the content is the per-side form  D ⊆ c12△c12.m  and  D ⊆ c21△c21.m.
Split: 'nonvacuous' = nest A keeps >=2 rhombi after x,y (the 162).   usage: python3 fibre.py 5 NW NCH"""
import sys
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
import os; os.chdir('/home/clio/projects/proofs/code-1120c2-Q597')
def work(w):
    st = Counter(); ns = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1; S = frozenset(bare(b))
            cs = sorted(corners(set(b.nest[0]['pos'])))
            for a in range(len(cs)):
                for c in range(a+1, len(cs)):
                    x, y = cs[a], cs[c]
                    c1, _ = mvl(b, x, 'B', DV); c2, _ = mvl(b, y, 'B', DV)
                    c12, _ = mvl(c1, y, 'B', DV); c21, _ = mvl(c2, x, 'B', DV)
                    B12, B21 = frozenset(bare(c12)), frozenset(bare(c21))
                    if B12 == B21: continue
                    S1, S2 = frozenset(bare(c1)), frozenset(bare(c2))
                    D = B12 ^ B21
                    nA = len(c12.nest[0]['pos']); cls = 'nonvacuous' if nA >= 2 else 'vacuous (nest A emptied by m)'
                    wx, wy = S ^ S1, S ^ S2; w4 = wx | wy | (S1 ^ B12) | (S2 ^ B21)
                    st[(cls, 'pairs')] += 1
                    st[(cls, '(a) D within wake(x) u wake(y)', D <= (wx | wy))] += 1
                    st[(cls, "(a') D within the four journey wakes", D <= w4)] += 1
                    ms = []
                    for z in corners(set(c12.nest[0]['pos'])):
                        for T in 'BC':
                            e12, _ = mvl(c12, z, T, DV); e21, _ = mvl(c21, z, T, DV)
                            F12, F21 = frozenset(bare(e12)), frozenset(bare(e21))
                            if F12 == F21: ms.append((F12, F21, z, T))
                    st[(cls, 'number of common closing labels m', len(ms))] += 1
                    for F12, F21, z, T in ms:
                        m12, m21 = B12 ^ F12, B21 ^ F21
                        st[(cls, '(b) union form [tautology]', D <= (m12 | m21))] += 1
                        st[(cls, '(b) per-side: D within wake of m on c12', D <= m12, 'on c21', D <= m21)] += 1
                        st[(cls, 'm target', T, '|D|', len(D), '|D & m12|', len(D & m12), '|D & m21|', len(D & m21))] += 1
                    # Hecke/braid pattern: x.y.x = y.x.y  (third step = first label again on each side)
                    e1, _ = mvl(c12, x, 'B', DV) if x in corners(set(c12.nest[0]['pos'])) else (None, None)
                    e2, _ = mvl(c21, y, 'B', DV) if y in corners(set(c21.nest[0]['pos'])) else (None, None)
                    st[(cls, 'braid x.y.x = y.x.y', e1 is not None and e2 is not None and frozenset(bare(e1)) == frozenset(bare(e2)))] += 1
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
