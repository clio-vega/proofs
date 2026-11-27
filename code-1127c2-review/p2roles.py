"""REVIEW 2026-11-27 c2: Lemma p2 (A1 tex l.129) with both role assignments, on every diverging B/C pair.
For a diverging pair (J_x -> B, J_y -> C on M0) and each ordered choice (clean, other):
  W    = M0 \\ (final board of other)
  fc   = least k with clean step k removing a tile of W   (None if never)
  p    = contact of clean (first index where clean and its re-run on M0.other differ)
  P2   : R_k (clean) is not a rhombus of other (incl. final), for all k <= fc
  base : R_0(clean) != R_0(other)        (hidden premise of the k=0 step)
  Claim: every M0 tile of U_fc (clean's start rhombus + hexagons 0..fc-1) is outside W  (proof's Claim, directly)
  CTRL : same P2 test over the whole clean journey (must fire, as in samepos.py)
Role tags: L = role at (1,a_T) at its contact with re-run (1,b) (the WB filter, stale_e.py); E the other.
usage: python3 p2roles.py SIZE NW NCH"""
import sys
exec(open("/home/clio/projects/proofs/code-1126-Q597/g3red.py").read().split("\ndef analyse(")[0])
def tilesof(bb): return {v[0:1] + (frozenset(v[1]),) for v in bb.tiles.values()}
def rhs(J): return [frozenset(rh(s[1])) for s in J] + [frozenset(rh(J[-1][2]))]
def check(clean, cleanR, other, Wo, M0, st, tag):
    mA, mB, p, q = middles(clean, cleanR)
    fc = next((k for k, s in enumerate(clean) if set(s[1]) & Wo), None)
    st[(tag, 'fc defined', fc is not None)] += 1
    if fc is None: return
    st[(tag, 'fc >= p', fc >= p)] += 1
    st[(tag, 'fc == p', fc == p)] += 1
    RC, RO = rhs(clean), set(rhs(other))
    st[(tag, 'base R0 != R0^o', RC[0] != rhs(other)[0])] += 1
    st[(tag, 'P2 k<=fc', not any(r in RO for r in RC[:fc+1]))] += 1
    st[(tag, 'P2 k<=p', not any(r in RO for r in RC[:p+1]))] += 1
    st[(tag, 'CTRL whole journey hits', any(r in RO for r in RC))] += 1
    st[(tag, 'CTRL k<=fc+1 hits', any(r in RO for r in RC[:fc+2]))] += 1
    U = {(('rh', frozenset(rh(clean[0][1]))))}
    for s in clean[:fc]: U |= set(s[1])
    st[(tag, 'Claim: U_fc M0-tiles outside W', not ((U & M0) & Wo))] += 1
    st[(tag, 'Claim: U_fc tiles all M0 or created by clean',
        all(t in M0 or any(t in s[2] for s in clean[:fc]) for t in U))] += 1
def work(w):
    st = Counter()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            M0 = tilesof(b)
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    b1, Jx = mvl(b, x, 'B', DV); b2, Jy = mvl(b, y, 'C', DV)
                    b12, Jyx = mvl(b1, y, 'C', DV); b21, Jxy = mvl(b2, x, 'B', DV)
                    if Jx == Jxy and Jy == Jyx: continue
                    Wx, Wy = M0 - tilesof(b1), M0 - tilesof(b2)
                    st[('pairs',)] += 1
                    # role tags as in stale_e.analyse2
                    info = []
                    for T, A, Bq in (('B', Jx, Jxy), ('C', Jy, Jyx)):
                        mA, mB, p, q = middles(A, Bq)
                        info.append(dict(T=T, p=p, tA=typed(A), tB=typed(Bq), n=len(A)))
                    Ls = [r for r in (0, 1) if info[r]['p'] < info[r]['n'] and info[r]['tA'][info[r]['p']] == (1, aT[info[r]['T']])]
                    geo = len(Ls) == 1 and info[Ls[0]]['tB'][info[Ls[0]]['p']] == (1, 'b')
                    J = {0: (Jx, Jxy, Wx), 1: (Jy, Jyx, Wy)}
                    for c in (0, 1):
                        o = 1 - c
                        if geo: role = 'clean=L,other=E' if c == Ls[0] else 'clean=E,other=L'
                        else: role = 'no-role(%d L)' % len(Ls)
                        check(J[c][0], J[c][1], J[o][0], J[o][2], M0, st, (role, 'clean->' + 'BC'[c]))
    return st
if __name__ == '__main__':
    tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for st in pool.imap_unordered(work, range(NCH)): tot.update(st)
    print('SIZE', SIZE)
    for k in sorted(tot, key=str): print(tot[k], k)
