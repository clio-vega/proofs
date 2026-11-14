"""Restrict to c>=2 fibres (same full boundary, >=2 puzzles). Compare R_ii with R_iii and with same-shape controls:
   C1 = R_A(+,E) [orientation (E_A,N_A), migrate toward C]; C2 = R_A(-,S) = kappa; C3 = R_ii composed with a fibre cycle."""
import sys, copy, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from flock import *
import migrate
from census import readout2, strings
from readouts import norm
from collections import Counter, defaultdict
nmax = int(sys.argv[1]); st = Counter()
for n in range(2, nmax+1):
    left, right, bottom = side_edges(n)
    for sl in strings(n):
        for sr in strings(n):
            fixed = {}
            for e, l in zip(left, sl): fixed[E(*e)] = l
            for e, l in zip(right, sr): fixed[E(*e)] = l
            byb = defaultdict(list)
            for P in puzzles(n, fixed): byb[tuple(P[E(*e)] for e in bottom)].append(P)
            for bot, Ps in byb.items():
                if len(Ps) < 2: continue
                rows = []
                for P in Ps:
                    b = Board(build_mosaic(n, P, SIGN, RULE))
                    if not b.nest[0]['pos']: continue
                    t0 = migrate.Global_ties[0]
                    _, R2, s2 = readout2(b, 0, 1, 'N'); _, R3, s3 = readout2(b, 0, -1, 'W')
                    _, C1, _ = readout2(b, 0, 1, 'E'); _, C2, _ = readout2(b, 0, -1, 'S')
                    rows.append((norm(R2), norm(R3), norm(C1), norm(C2), migrate.Global_ties[0]-t0, s2, s3))
                if not rows: continue
                st['fibres'] += 1; st['elements'] += len(rows); st[f'fibre size {len(rows)}'] += 1
                st['Rii injective on fibre'] += len(set(r[0] for r in rows)) == len(rows)
                for i, r in enumerate(rows):
                    st['Rii=Riii'] += r[0] == r[1]
                    st['C1 same cells as Rii'] += set(dict(r[2])) == set(dict(r[0])); st['C1 = Rii'] += r[2] == r[0]
                    st['C2 same cells as Rii'] += set(dict(r[3])) == set(dict(r[0])); st['C2 = Rii'] += r[3] == r[0]
                    st['C3 (fibre-shifted Rii) = Riii'] += rows[(i+1) % len(rows)][0] == r[1]
                    st['ties met'] += r[4]; st['a-moves'] += r[5]['a'] + r[6]['a']
    print(n, dict(st), flush=True)
