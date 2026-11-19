"""PROVE 11-19 c2: verify the run-monotonicity lemma on every type-a run of every B and C journey (<=5).
(a) for both edge classes e of the run's rhombus, cross(phi(e), phi(Delta R_j)) has constant nonzero sign along the run;
(b) survival within the run: no tile created by step j is removed by a later step j' of the same run;
(c) survival over the rest of the journey (stronger, not claimed by the lemma)."""
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585')
from mosaic import Geo
GG = Geo(-1)
def ph(v): return GG.phi(v) - GG.phi((0,0,0,0))
def cross(z, w): return (z.conjugate()*w).imag
def work(w):
    st = Counter(); ns = 0; seenJ = set()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1
            for x in corners(set(b.nest[0]['pos'])):
                for TT in ('B', 'C'):
                    _, J = mvl(b, x, TT, DV)
                    sig = tuple(J)
                    if sig in seenJ: continue
                    seenJ.add(sig)
                    # maximal type-a runs
                    i = 0
                    while i < len(J):
                        if J[i][0] != 'a': i += 1; continue
                        k = i
                        while k < len(J) and J[k][0] == 'a': k += 1
                        run = J[i:k]
                        R0 = rhof(run[0]); edges = [vec(e) for e in tedges(R0)]
                        cl = []
                        for e in edges:
                            if not any(par(e, c) for c in cl): cl.append(e)
                        for c in cl:
                            sgns = set()
                            for s in run:
                                dR = disp(s, s); z = cross(ph(c), ph(dR))
                                sgns.add(0 if abs(z) < 1e-9 else (1 if z > 0 else -1))
                            st[('(a) sign set of cross(e, dR) along run', tuple(sorted(sgns)))] += 1
                        made = {}
                        for j, (_, g, n) in enumerate(run):
                            for t in g:
                                if t in made and t[0] != 'rh': st['(b) VIOLATION within run (non-rh tile made earlier in run is removed)'] += 1
                            for t in n: made[t] = j
                        st[('(b) runs checked', len(run) > 1)] += 1
                        # (c) rest of the journey
                        madeall = set()
                        for j, (_, g, n) in enumerate(J):
                            if j >= i and j < k: madeall |= set(n) - {rhof((None, n))} if False else set(t for t in n if t[0] != 'rh')
                            if j >= k and (set(g) & madeall): st['(c) a run-made tile removed later in the journey'] += 1; break
                        else: st['(c) run-made tiles all survive the journey'] += 1
                        i = k
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
