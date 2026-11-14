"""Smallest c>=2 fibre (n=6): draw both mosaics; label alpha rhombi '(ii)k/(iii)k' processing indices; print arrivals."""
import sys, copy
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from readouts import *
from draw import draw
from census import strings
from lemma3 import arrivals
from collections import defaultdict
n = 6; left, right, bottom = side_edges(n); k = 0
for sl in strings(n):
    for sr in strings(n):
        fixed = {E(*e): l for e, l in zip(left, sl)}; fixed.update({E(*e): l for e, l in zip(right, sr)})
        byb = defaultdict(list)
        for P in puzzles(n, fixed): byb[tuple(P[E(*e)] for e in bottom)].append(P)
        for bot, Ps in byb.items():
            if len(Ps) < 2: continue
            for P in Ps:
                b = Board(build_mosaic(n, P, SIGN, RULE)); A = b.nest[0]
                if not A['pos']: continue
                u, w = A['u'], A['w']
                lab = {}
                for sgn, dirn, nm in ((1, 'N', 'ii'), (-1, 'W', 'iii')):
                    E_ = tuple(sgn*a for a in u); N_ = tuple(sgn*a for a in w); fl = unique_flock(b, 0, E_, N_)
                    key = {t: (v, box_of(b, 0, b.tiles[t][1], E_, N_)[0]) for t, v in fl.items()}
                    order = sorted(fl, key=lambda t: key[t], reverse=(dirn in 'NE'))
                    for i, t in enumerate(order): lab[t] = lab.get(t, '') + f'{nm}{i+1}:{fl[t]} '
                print('mosaic', k, 'boundary', sl, sr, bot, 'alpha cells', sorted(A['pos']))
                for t, s in lab.items(): print('   ', box_of(b, 0, b.tiles[t][1], u, w), s)
                print('   arrivals ii', arrivals(copy.deepcopy(b), 1, 'N')); print('   arrivals iii', arrivals(copy.deepcopy(b), -1, 'W'))
                draw(b, f'/tmp/small{k}.png', labels={t: s.split()[0][2:] for t, s in lab.items()}); k += 1
