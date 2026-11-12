"""Lesion of Cor 3.4: G = unique flock on the remainder in the STANDARD orientation (E_B,N_B), migrated EAST
(a legal migration, not Purbhoo's recipe).  If the comparison has power, this must differ from C somewhere."""
from compare import *
from betaq import beta_q
import flock
from collections import Counter
neg=lambda v: tuple(-a for a in v)
def cor34_lesion(f, lam, k0, k1):
    b = flock.canonical(lam, k0, k1); B = b.nest[1]; EB, NB = B['u'], B['w']
    F = {B['pos'][(c, r)]: v for (c, r), v in f.items()}
    mu_pos = [ij for ij in B['pos'] if ij not in f]
    boxesG = {flock.box_of(b, 1, b.tiles[B['pos'][ij]][1], EB, NB): ij for ij in mu_pos}
    mu_shape = tuple(sorted(Counter(r for _, r in mu_pos).values(), reverse=True))
    Ls = flock.lr_tableaux(set(boxesG), mu_shape); assert len(Ls) == 1
    Gfl = {B['pos'][boxesG[bx]]: v for bx, v in Ls[0].items()}
    flock.migrate_flock(b, 1, F, EB, NB, 'E')
    d2, G2, (E2,N2), boxG, st = flock.migrate_flock(b, 1, Gfl, EB, NB, 'E')
    # read the output in the frame (-E2,-N2) so that it sits in the positive quadrant like Purbhoo's G'
    return {flock.box_of(b, d2, b.tiles[R][1], neg(E2), neg(N2)): v for R, v in G2.items()}
S=Counter()
for a in range(1,6):
  for c_ in range(1,6):
    for lam in parts(a,NC):
      for mu in parts(c_,NC):
        for nu,bs in LR(lam,mu,NC).items():
          if len(bs)<2: continue
          nu_t=trim(nu); S['cells']+=1; diff=False
          for b in bs:
            T=beta_q(b,lam,mu); cb=C1(b,lam,mu,NC)[0]
            try: out=cor34_lesion(T,nu_t,nu_t[0]+1,len(nu_t)+1)
            except Exception as ex: S['lesion_crash']+=1; diff=True; continue
            if out!=beta_q(cb,mu,lam): diff=True
            if not is_lr(out): S['lesion_nonLR']+=1
          S['lesion_differs_from_C']+=diff
print(dict(S))
