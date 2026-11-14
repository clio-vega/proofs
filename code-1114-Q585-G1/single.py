"""Single-rhombus coherence: for each mosaic M (c>=2 census) and each outer corner rho of alpha:
 path X: rho A->C (dir E_A).   path Y: rho A->B (dir N_A), then the landed rhombus from B in each compass dir of B's frame.
 Does some Y land in C with the same mosaic as X?"""
import sys, copy
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts
from transit import tiles_of
from collections import Counter
def corners(S): return [(i, j) for (i, j) in S if (i+1, j) not in S and (i, j+1) not in S]
def mig1(b, X, ij, dvec):
    t = b.nest[X]['pos'].pop(ij); R, Y, ij2, s = b.journey(t, X, unit_dir(dvec)); return Y, ij2
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
        A = b0.nest[0]; uA, wA = A['u'], A['w']; uB, wB = b0.nest[1]['u'], b0.nest[1]['w']
        for rho in corners(set(A['pos'])):
            st['pairs'] += 1
            try:
                bx = copy.deepcopy(b0); dX = mig1(bx, 0, rho, uA)
                by = copy.deepcopy(b0); dY = mig1(by, 0, rho, wA)
            except Exception as ex: st['ERR1 '+str(ex)[:30]] += 1; continue
            st[f'X dst {dX[0]}  Y dst {dY[0]}'] += 1
            ok = False
            for nm, dv in (('+u', uB), ('-u', tuple(-x for x in uB)), ('+w', wB), ('-w', tuple(-x for x in wB))):
                bz = copy.deepcopy(by)
                try: dZ = mig1(bz, 1, dY[1], dv)
                except Exception as ex: st[f'Y2 {nm} ERR'] += 1; continue
                st[f'Y2 {nm} dst {dZ[0]}'] += 1
                if dZ[0] == 2 and tiles_of(bz) == tiles_of(bx): st[f'Y2 {nm} == X'] += 1; ok = True
            st['some Y2 == X'] += ok
print(dict(st))
