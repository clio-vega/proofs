"""Per-rhombus comparison of R_ii and R_iii: start box (in + frame of A) -> (value, final box, step string)."""
import sys, copy, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from flock import *
import migrate
from readouts import unique_flock, norm
from census import strings
from collections import Counter

def trace(board, sgn, dirn):
    b = copy.deepcopy(board); X = 0; u, w = b.nest[X]['u'], b.nest[X]['w']
    E_ = tuple(sgn*a for a in u); N_ = tuple(sgn*a for a in w)
    fl = unique_flock(b, X, E_, N_)
    startbox = {t: box_of(b, 0, b.tiles[t][1], u, w) for t in fl}   # in + frame always
    dvec = {'E': E_, 'N': N_, 'W': tuple(-a for a in E_), 'S': tuple(-a for a in N_)}[dirn]
    d = unit_dir(dvec)
    key = {t: (v, box_of(b, 0, b.tiles[t][1], E_, N_)[0]) for t, v in fl.items()}
    order = sorted(fl, key=lambda t: key[t], reverse=(dirn in 'NE'))
    inv = {t: ij for ij, t in b.nest[0]['pos'].items()}
    recs = []; finals = {}
    for t in order:
        del b.nest[0]['pos'][inv[t]]
        R = t; seq = ''
        while True:
            vs = b.tiles[R][1]; done = False
            for Xd, info in enumerate(b.nest):
                if Xd == 0: continue
                ij = b.packed_at(info, vs)
                if ij is not None and ij not in info['pos'] and b.is_ideal(set(info['pos']) | {ij}):
                    info['pos'][ij] = R; done = True; dst = Xd; break
            if done: break
            R = b.step(R, d); seq += b.last
        finals[R] = (startbox[t], fl[t], seq)
    E2, N2, ang = induced(b, dst, E_, N_)
    out = []
    for R, (sb, v, seq) in finals.items():
        out.append((sb, v, box_of(b, dst, b.tiles[R][1], E2, N2), seq))
    return dst, out, [startbox[t] for t in order]

if __name__ == '__main__':
    nmax = int(sys.argv[1]); st = Counter(); shown = 0
    for n in range(2, nmax+1):
        left, right, bottom = side_edges(n)
        for sl in strings(n):
            for sr in strings(n):
                fixed = {}
                for e, l in zip(left, sl): fixed[E(*e)] = l
                for e, l in zip(right, sr): fixed[E(*e)] = l
                for P in puzzles(n, fixed):
                    b = Board(build_mosaic(n, P, SIGN, RULE))
                    if not b.nest[0]['pos']: continue
                    d2, o2, ord2 = trace(b, 1, 'N'); d3, o3, ord3 = trace(b, -1, 'W')
                    st['n'] += 1
                    m2 = {sb: (v, fb, s) for sb, v, fb, s in o2}; m3 = {sb: (v, fb, s) for sb, v, fb, s in o3}
                    st['same start->final box map'] += all(m2[k][1] == m3[k][1] for k in m2)
                    st['same start->stepstring'] += all(m2[k][2] == m3[k][2] for k in m2)
                    st['same values per start box'] += all(m2[k][0] == m3[k][0] for k in m2)
                    st['same processing order'] += ord2 == ord3
                    st['same multiset of stepstrings'] += sorted(s for *_, s in o2) == sorted(s for *_, s in o3)
                    if shown < 3 and len(o2) >= 3:
                        shown += 1; print('ii ', sorted(o2)); print('iii', sorted(o3)); print(ord2, ord3)
        print(n, dict(st), flush=True)
