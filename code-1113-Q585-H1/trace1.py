"""Trace Cor 3.4 journeys: per rhombus, source box (frame of its flock), value, destination box (induced frame)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from flock import *
from switch import sigma
from compare import beta_cm, trim
from commutor import LR, parts

def traced(board, src, flock, E, N, direction):
    dvec = {'E': E, 'N': N, 'W': tuple(-a for a in E), 'S': tuple(-a for a in N)}[direction]
    d = unit_dir(dvec)
    key = {t: (v, box_of(board, src, board.tiles[t][1], E, N)[0]) for t, v in flock.items()}
    srcbox = {t: box_of(board, src, board.tiles[t][1], E, N) for t in flock}
    order = sorted(flock, key=lambda t: key[t], reverse=(direction in 'NE'))
    inv = {t: ij for ij, t in board.nest[src]['pos'].items()}
    recs = []
    for t in order:
        del board.nest[src]['pos'][inv[t]]
        R, X, ij, st = board.journey(t, src, d)
        recs.append([srcbox[t], flock[t], R, X, dict(st)])
    E2, N2, ang = induced(board, X, E, N)
    for r in recs: r.append(box_of(board, X, board.tiles[r[2]][1], E2, N2))
    return recs

def run(T, nu):
    k0 = max(nu[0],1)+1; k1 = len(nu)+1
    b = canonical(nu, k0, k1); B = b.nest[1]; EB, NB = B['u'], B['w']; neg = lambda v: tuple(-a for a in v)
    F = {B['pos'][(c, r)]: v for (c, r), v in T.items()}
    mu_pos = [ij for ij in B['pos'] if ij not in T]
    boxesG = {box_of(b, 1, b.tiles[B['pos'][ij]][1], neg(EB), neg(NB)): ij for ij in mu_pos}
    mu_shape = tuple(sorted(Counter(r for _, r in mu_pos).values(), reverse=True))
    Ls = lr_tableaux(set(boxesG), mu_shape) if mu_pos else [{}]
    Gfl = {B['pos'][boxesG[bx]]: v for bx, v in Ls[0].items()}
    rF = traced(b, 1, F, EB, NB, 'E')
    rG = traced(b, 1, Gfl, neg(EB), neg(NB), 'W')
    return rF, rG, b

if __name__ == '__main__':
    lam, mu = (2,1), (2,1)
    for nu, bs in LR(lam, mu, 4).items():
        if len(bs) < 2: continue
        nt = trim(nu)
        for bb in bs:
            T = beta_cm(bb, lam, mu, nt)
            rF, rG, _ = run(T, nt)
            print('nu', nt, 'T', sorted(T.items()), 'Sigma', sorted(sigma(T, nt).items()))
            for r in rF: print('  F', r[0], r[1], '->', r[5], r[4])
            for r in rG: print('  G', r[0], r[1], '->', r[5], r[4])
