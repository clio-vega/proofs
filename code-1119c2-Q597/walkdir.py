"""PROVE 11-19 c2 (walkdir): set of midpoint displacements of pi_i(R) per journey type; and per-cell crossing counts. Originally: product-picture reformulation of no-repeat.  For every step of every distinct B/C journey (<=5):
(i)  some 0-triangle moves  <=>  pi_0(R) changes;  same for 1-triangles and pi_1(R);
(ii) when 0-triangles move, pi_0(R) before and after are edges of the union U of their pi_0-cells (lozenge or one triangle), and differ;
(iii) the moved 0-cells' pi_0 sets are unchanged as a SET of cells (cells are fixed; only q-values move)."""
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
def cells(S, i):
    out = []
    for k, vs in S:
        if k != 'tri': continue
        a = frozenset(v[:2] for v in vs); b = frozenset(v[2:] for v in vs)
        if i == 0 and len(b) == 1: out.append(a)
        if i == 1 and len(a) == 1: out.append(b)
    return out
def proj(R, i): return frozenset((v[:2] if i == 0 else v[2:]) for v in R[1])
def work(w):
    st = Counter(); ns = 0; seenJ = set()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1
            for x in corners(set(b.nest[0]['pos'])):
                for TT in ('B', 'C'):
                    _, J = mvl(b, x, TT, DV)
                    if not J or tuple(J) in seenJ: continue
                    seenJ.add(tuple(J))
                    for typ, g, n in J:
                        Rb = next(t for t in g if t[0] == 'rh'); Ra = next(t for t in n if t[0] == 'rh')
                        for i in (0, 1):
                            mb = [sum(v[k] for v in proj(Rb, i)) for k in (0, 1)]; ma_ = [sum(v[k] for v in proj(Ra, i)) for k in (0, 1)]
                            dv = (ma_[0]-mb[0], ma_[1]-mb[1])
                            if dv != (0, 0): st[(TT, 'lattice', i, 'step', typ, '2*midpoint displacement', str(dv))] += 1
                            eb = sorted(proj(Rb, i)); ev = (eb[1][0]-eb[0][0], eb[1][1]-eb[0][1]); ev = max(ev, (-ev[0], -ev[1]))
                            st[(TT, 'lattice', i, 'edge class of pi(R)', str(ev))] += 1
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
