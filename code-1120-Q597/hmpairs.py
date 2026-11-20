"""PROVE 11-20 retreat, Q598(d'): classify every corner pair (x,y) of nest A (<=5 census sampler) by
 - the pair type: B/B (x->B,y->B) and B/C (x->B,y->C);
 - commutation: the two orders give equal tilings;
 - the corner geometry: di=|i_x-i_y|, dj=|j_x-j_y| in nest coordinates (cells (i,j) = X + i u + j w), and whether the
   corners are 'adjacent' in the Hage-Malbos sense (consecutive corners of the ideal: no corner strictly between them).
usage: python3 hmpairs.py 5 NW NCH"""
import sys
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
import os; os.chdir('/home/clio/projects/proofs/code-1120-Q597')
def work(w):
    st = Counter(); ns = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1
            cs = sorted(corners(set(b.nest[0]['pos'])))
            for a in range(len(cs)):
                for c in range(len(cs)):
                    if a == c: continue
                    x, y = cs[a], cs[c]
                    adj = abs(a - c) == 1
                    geo = (abs(x[0]-y[0]), abs(x[1]-y[1]))
                    for TT in (('B', 'B'), ('B', 'C')):
                        if TT == ('B', 'B') and a > c: continue
                        c1, _ = mvl(b, x, TT[0], DV); c2, _ = mvl(b, y, TT[1], DV)
                        c12, _ = mvl(c1, y, TT[1], DV); c21, _ = mvl(c2, x, TT[0], DV)
                        com = bare(c12) == bare(c21)
                        st[('/'.join(TT), 'commute' if com else 'NON-commute', 'consecutive corners' if adj else 'non-consecutive',
                            'min(di,dj)', min(geo), 'max', max(geo) if max(geo) < 3 else '>=3')] += 1
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
