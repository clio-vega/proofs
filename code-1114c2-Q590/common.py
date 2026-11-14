"""Shared census for Q590 (c>=2 cells, |lambda|,|mu|<=5, gl4) -- identical generator to code-1114-Q585-G1/lemma3.py."""
import sys, copy
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
sys.path.insert(0,'/home/clio/projects/proofs/code-1114-Q585-G1')
from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts
from transit import tiles_of
def census():
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
            yield (lam, mu, nu, bb), b0
def corners(S): return [(i, j) for (i, j) in S if (i+1, j) not in S and (i, j+1) not in S]
def mig1(b, X, ij, dvec):
    t = b.nest[X]['pos'].pop(ij); R, Y, ij2, s = b.journey(t, X, unit_dir(dvec)); return Y, ij2, s
