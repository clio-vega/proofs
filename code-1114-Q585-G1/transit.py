"""Transitivity test: R_ii lands flock F' in B (induced frame (E',N')). Migrate F' onward B->C (try each compass dir).
Does the result equal M_iii (same tiles) and the same tableau as R_iii?"""
import sys, copy
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from readouts import *
import migrate
from compare import beta_cm, trim
from commutor import LR, parts
from collections import Counter
def tiles_of(b): return sorted((k, tuple(sorted(vs))) for k, vs in b.tiles.values())
def go(b, sgn, dirn):
    u, w = b.nest[0]['u'], b.nest[0]['w']; E_ = tuple(sgn*a for a in u); N_ = tuple(sgn*a for a in w)
    fl = unique_flock(b, 0, E_, N_)
    return migrate_flock(b, 0, fl, E_, N_, dirn)
if __name__ == '__main__':
    st = Counter()
    for a in range(1, 6):
     for c_ in range(1, 6):
      for lam in parts(a, 4):
       for mu in parts(c_, 4):
        for nu, bs in LR(lam, mu, 4).items():
          if len(bs) < 2: continue
          nt = trim(nu)
          for bb in bs:
            T = beta_cm(bb, lam, mu, nt)
            b0 = canonical(nt, max(nt[0], 1)+1, len(nt)+1); B = b0.nest[1]
            migrate_flock(b0, 1, {B['pos'][(c, r)]: v for (c, r), v in T.items()}, B['u'], B['w'], 'E')
            b2 = copy.deepcopy(b0); d2, F2, (E2, N2), box2, _ = go(b2, 1, 'N')
            b3 = copy.deepcopy(b0); d3, F3, (E3, N3), box3, _ = go(b3, -1, 'W')
            st['n'] += 1
            for dirn in 'NESW':
                bb2 = copy.deepcopy(b2)
                try:
                    d4, F4, fr4, box4, _ = migrate_flock(bb2, 1, F2, E2, N2, dirn)
                except Exception as ex:
                    st[f'{dirn}: ERR {str(ex)[:30]}'] += 1; continue
                st[f'{dirn}: dst={d4}'] += 1
                if d4 == 2:
                    st[f'{dirn}: same mosaic as M_iii'] += tiles_of(bb2) == tiles_of(b3)
                    st[f'{dirn}: same tableau as R_iii'] += norm(box4) == norm(box3)
                    st[f'{dirn}: same frame as R_iii'] += fr4[:2] == (E3, N3)
    print(dict(st))
