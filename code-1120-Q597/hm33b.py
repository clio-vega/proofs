"""PROVE 11-20 retreat Q598(d'): do the non-commuting B/B pairs close as 3+3 hexagons (HM 2012.15649 l.2033 closes the
critical branching in <=3+3 column steps)?  1-cells = single migrations (any current nest-A corner -> B or C).
For each non-commuting pair (x,y): S12 = {tiling of c12 . m}, S21 = {tiling of c21 . m'} over all 1-cells m, m';
closes_3+3 iff S12 & S21 nonempty.  Control: commuting B/B pairs (closing 2+2, should also close 3+3 trivially, m=m').
usage: python3 hm33.py 5 NW NCH"""
import sys
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
import os; os.chdir('/home/clio/projects/proofs/code-1120-Q597')
def one_more(c, DV):
    out = {}
    for z in corners(set(c.nest[0]['pos'])):
        for T in 'BC':
            cz, _ = mvl(c, z, T, DV); out[(z, T)] = frozenset(bare(cz))
    return out
def work(w):
    st = Counter(); ns = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1
            cs = sorted(corners(set(b.nest[0]['pos'])))
            for a in range(len(cs)):
                for c in range(a+1, len(cs)):
                    x, y = cs[a], cs[c]
                    c1, _ = mvl(b, x, 'B', DV); c2, _ = mvl(b, y, 'B', DV)
                    c12, _ = mvl(c1, y, 'B', DV); c21, _ = mvl(c2, x, 'B', DV)
                    com = bare(c12) == bare(c21)
                    A, Bm = one_more(c12, DV), one_more(c21, DV)
                    inter = set(A.values()) & set(Bm.values())
                    same_label = [k for k in A if k in Bm and A[k] == Bm[k]]
                    nA = len(c12.nest[0]['pos'])
                    diff_label = [(k, k2) for k in A for k2 in Bm if k != k2 and A[k] == Bm[k2]]
                    st[('B/B', 'commute' if com else 'NON-commute', 'closes 3+3 same label', bool(same_label),
                        'A-rhombi left after x,y', nA if nA < 3 else '>=3', 'NEG CONTROL: closes with different labels', bool(diff_label))] += 1
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
