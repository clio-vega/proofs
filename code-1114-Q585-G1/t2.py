"""T2 per rhombus: in the canonical mosaic M_ii=(empty, beta', gamma) carrying F', migrate F' B->C (dir S).
For each rhombus: start box in B (frame E2,N2) vs final box in C (induced frame). Also step strings and per-rhombus 'shape'."""
import sys, copy
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from readouts import *
import migrate
from compare import beta_cm, trim
from commutor import LR, parts
from collections import Counter
from transit import go
st = Counter(); shown = 0
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
        b = copy.deepcopy(b0); d2, F2, (E2, N2), box2, _ = go(b, 1, 'N')
        start = {t: box_of(b, 1, b.tiles[t][1], E2, N2) for t in F2}
        key = {t: (v, start[t][0]) for t, v in F2.items()}
        order = sorted(F2, key=lambda t: key[t])          # S migration: min first
        inv = {t: ij for ij, t in b.nest[1]['pos'].items()}
        d = unit_dir(tuple(-x for x in N2)); land = {}; seqs = []
        for t in order:
            del b.nest[1]['pos'][inv[t]]
            R, X, ij, s = b.journey(t, 1, d); land[t] = R; seqs.append(dict(s))
        E3, N3, _ = induced(b, 2, E2, N2)
        fin = {t: box_of(b, 2, b.tiles[land[t]][1], E3, N3) for t in order}
        st['n'] += 1
        # normalise both box sets by the same translation? compare raw and up to a common translation
        dx = {(fin[t][0]-start[t][0], fin[t][1]-start[t][1]) for t in order}
        st['per-rhombus constant translation'] += len(dx) == 1
        st['per-rhombus identical box'] += dx == {(0, 0)}
        if shown < 2: shown += 1; print([(start[t], F2[t], fin[t], seqs[i]) for i, t in enumerate(order)])
print(dict(st))
