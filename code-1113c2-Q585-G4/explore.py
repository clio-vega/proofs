"""Look at the actual readouts: Tao(P), kappa(P)=T, Phi_G(P)=Sigma(T), with shapes and contents, for the 20-element census."""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0, '/home/clio/projects/proofs/code-1113-Q585-H1')
from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts
from switch import sigma
from flock import skew_of, content, is_lr

def mosaic_of(T, nt):
    k0 = max(nt[0], 1) + 1; k1 = len(nt) + 1
    b = canonical(nt, k0, k1); B = b.nest[1]
    F = {B['pos'][(c, r)]: v for (c, r), v in T.items()}
    migrate_flock(b, 1, F, B['u'], B['w'], 'E')
    return b, k0, k1

def elements(amax=4, bmax=4, cmin=2):
    for a in range(1, amax + 1):
        for c_ in range(1, bmax + 1):
            for lam in parts(a, 4):
                for mu in parts(c_, 4):
                    for nu, bs in LR(lam, mu, 4).items():
                        if len(bs) < cmin: continue
                        nt = trim(nu)
                        for bb in bs: yield lam, mu, nt, beta_cm(bb, lam, mu, nt)

def show(T):
    return ' / '.join(''.join(str(v) for c, v in sorted((c, v) for (c, rr), v in dict(norm(T)).items() if rr == r))
                      for r in sorted({r for (c, r), v in norm(T)}))

if __name__ == '__main__':
    for lam, mu, nt, T in elements():
        b, k0, k1 = mosaic_of(T, nt)
        _, tao = readout(b, 0, 1, 'N'); _, PG = readout(b, 1, -1, 'W')
        print('lam', trim(lam), 'mu', trim(mu), 'nu', nt, 'box', k1, 'x', k0)
        print('   T  ', skew_of(set(T)), content(T), show(T))
        print('   tao', skew_of(set(dict(norm(tao)))), content(tao), is_lr(dict(norm(tao))), show(tao))
        print('   PG ', skew_of(set(dict(norm(PG)))), content(PG), show(PG), ' Sigma ok', norm(PG) == norm(sigma(T, nt)))
