"""Q585 G1 Step 0: R_ii = R_A(+,N) vs R_iii = R_A(-,W) on EVERY KT puzzle of size n (all boundaries), plus lesion kappa = R_A(-,S).
Prints: agreement counts, lesion firing, step-type histogram (a/b moves), tie counts, dst nests."""
import sys, copy, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from flock import *
import migrate
from readouts import unique_flock, norm
from collections import Counter

def readout2(board, X, sgn, dirn):
    b = copy.deepcopy(board); u, w = b.nest[X]['u'], b.nest[X]['w']
    E_ = tuple(sgn*a for a in u); N_ = tuple(sgn*a for a in w)
    fl = unique_flock(b, X, E_, N_)
    if not fl: return None, {}, Counter()
    dst, out, fr_, boxes, st = migrate_flock(b, X, fl, E_, N_, dirn)
    return dst, boxes, st

def strings(n):
    for k in range(n+1):
        for S in itertools.combinations(range(n), k):
            yield tuple(1 if i in S else 0 for i in range(n))

if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    st = Counter(); fails = []
    for n in range(2, nmax+1):
        left, right, bottom = side_edges(n)
        for sl in strings(n):
            for sr in strings(n):
                if sum(sl) != sum(sr): pass
                fixed = {}
                for e, l in zip(left, sl): fixed[E(*e)] = l
                for e, l in zip(right, sr): fixed[E(*e)] = l
                for P in puzzles(n, fixed):
                    M = build_mosaic(n, P, SIGN, RULE); b = Board(M)
                    if not b.nest[0]['pos']: st['alpha empty'] += 1; continue
                    t0 = migrate.Global_ties[0]
                    try:
                        d2, R2, s2 = readout2(b, 0, 1, 'N')
                        d3, R3, s3 = readout2(b, 0, -1, 'W')
                        dk, K, sk = readout2(b, 0, -1, 'S')
                    except Exception as ex:
                        st['ERR '+type(ex).__name__+' '+str(ex)[:40]] += 1; continue
                    st['n'] += 1; st[f'dst ii={d2} iii={d3} k={dk}'] += 1
                    st['ties met'] += migrate.Global_ties[0]-t0
                    st['steps a (ii)'] += s2['a']; st['steps b (ii)'] += s2['b']
                    st['steps a (iii)'] += s3['a']; st['steps b (iii)'] += s3['b']
                    if s2['a'] + s3['a'] > 0: st['elements with a-moves'] += 1
                    if norm(R2) == norm(R3): st['Rii=Riii'] += 1
                    else: fails.append((n, sl, sr, norm(R2), norm(R3)))
                    if norm(K) != norm(R2): st['lesion kappa fires (K!=Rii)'] += 1
                    if len(set(R2.values())) >= 2: st['nontrivial content (>=2 labels)'] += 1
                    if norm(R2) == norm(R3) and len(set(R2.values())) >= 2 and norm(K) != norm(R2): st['agree & lesion fires & >=2 labels'] += 1
        print(n, dict(st), flush=True)
    print('FAILS', len(fails)); [print(f) for f in fails[:10]]
