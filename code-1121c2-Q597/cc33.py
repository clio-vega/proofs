"""PROVE 11-21 c2 retreat, Q598: the C/C census (the control absent from every table, DREAM 11-21).
From code-1120-Q597/hm33b.py, unchanged machinery.  For every pair of distinct nest-A corners x<y and kind in {BB, CC}:
commute?  For non-commuting pairs: closes 3+3 with a common third migration (same label (z,T) on both sides)?
Which targets T occur among the closing labels (B only / C only / both)?  Labels searched in BOTH orders (B,C) and (C,B)
as the search-order lesion.  Non-emptiness statistic: number of C/C pairs, number non-commuting.
usage: python3 cc33.py 5 NW NCH"""
import sys
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
def one_more(c, DV, order):
    out = {}
    for z in corners(set(c.nest[0]['pos'])):
        for T in order:
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
                    for K in ('B', 'C'):
                        c1, _ = mvl(b, x, K, DV); c2, _ = mvl(b, y, K, DV)
                        c12, _ = mvl(c1, y, K, DV); c21, _ = mvl(c2, x, K, DV)
                        com = bare(c12) == bare(c21)
                        st[(K + '/' + K, 'pairs')] += 1
                        st[(K + '/' + K, 'commute', com)] += 1
                        if com: continue
                        for order in ('BC', 'CB'):
                            A, Bm = one_more(c12, DV, order), one_more(c21, DV, order)
                            same = sorted({k[1] for k in A if k in Bm and A[k] == Bm[k]})
                            diff = any(A[k] == Bm[k2] for k in A for k2 in Bm if k != k2)
                            nA = len(c12.nest[0]['pos'])
                            st[(K + '/' + K, 'NON-commute', 'order', order, 'closes 3+3 same label', bool(same),
                                'third-step targets', ''.join(same), 'A-rhombi left', nA if nA < 3 else '>=3',
                                'NEG CONTROL different labels close', diff)] += 1
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
