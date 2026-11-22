"""PROVE 11-22: check of Lemma side (translation invariance of the pq/dl readout) on the census steps.
For every step s of every clean journey's middle, and v in {w_A = (0,0;0,1) in Lambda_1, u_A = (0,-1;0,0) in Lambda_0}:
   dl-dict of s (pure-cell moves, as bs6.sdelta) vs that of tr(s, v):
   v in Lambda_1 => side-0 entries identical, side-1 keys shifted;  v in Lambda_0 => side-1 identical, side-0 keys shifted.
Control: the opposite side must differ whenever the step moves a triangle of that side.   usage: python3 sidelemma.py"""
import os
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def work(w):
    st = Counter(); nstates = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            nstates += 1
            cs = corners(set(b.nest[0]['pos']))
            for x in cs[:2]:
                for T in 'BC':
                    b1, J = mvl(b, x, T, DV)
                    for s in J:
                        d0 = sdelta(s)
                        if d0 is None: continue
                        for nm, v, keep in (('w_A in L1', (0, 0, 0, 1), 0), ('u_A in L0', (0, -1, 0, 0), 1)):
                            d1 = sdelta(tr(s, v))
                            same = {k: d0[k] for k in d0 if k[0] == keep} == {k: d1[k] for k in d1 if k[0] == keep}
                            other = [k for k in d0 if k[0] != keep]
                            diff = {k: d0[k] for k in d0 if k[0] != keep} != {k: d1[k] for k in d1 if k[0] != keep}
                            st[(nm, 'side', keep, 'identical', same)] += 1
                            if other: st[(nm, 'CONTROL other side differs', diff)] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
