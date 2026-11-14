"""c>=2 census (|lam|,|mu|<=5, <=4 rows; mosaics built by inverting Cor 3.2 from each LR tableau T).
R_ii vs R_iii, with same-fibre controls C1 = R_A(+,E), C2 = kappa, C3 = fibre-shifted R_ii; tie and move statistics."""
import sys, time
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from readouts import *
import migrate
from compare import beta_cm, trim
from commutor import LR, parts
from census import readout2
from collections import Counter
st = Counter(); t0 = time.time(); N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
for a in range(1, N+1):
 for c_ in range(1, N+1):
  for lam in parts(a, 4):
   for mu in parts(c_, 4):
    for nu, bs in LR(lam, mu, 4).items():
      if len(bs) < 2: continue
      nt = trim(nu); st['cells'] += 1; rows = []
      for bb in bs:
        T = beta_cm(bb, lam, mu, nt)
        b = canonical(nt, max(nt[0], 1)+1, len(nt)+1); B = b.nest[1]
        migrate_flock(b, 1, {B['pos'][(c, r)]: v for (c, r), v in T.items()}, B['u'], B['w'], 'E')
        t1 = migrate.Global_ties[0]
        _, R2, s2 = readout2(b, 0, 1, 'N'); _, R3, s3 = readout2(b, 0, -1, 'W')
        ties = migrate.Global_ties[0]-t1
        _, C1, _ = readout2(b, 0, 1, 'E'); _, C2, _ = readout2(b, 0, -1, 'S')
        rows.append((norm(R2), norm(R3), norm(C1), norm(C2))); st['ties in ii+iii'] += ties
        st['a-moves'] += s2['a']+s3['a']; st['b-moves'] += s2['b']+s3['b']; st['elements with ties'] += ties > 0
      st['n'] += len(rows); st['Rii injective on fibre'] += len(set(r[0] for r in rows)) == len(rows)
      for i, r in enumerate(rows):
        st['Rii=Riii'] += r[0] == r[1]
        st['C1 same cells'] += set(dict(r[2])) == set(dict(r[0])); st['C1 fires'] += r[2] != r[1]
        st['C2 same cells'] += set(dict(r[3])) == set(dict(r[0])); st['C2 fires'] += r[3] != r[1]
        st['C3 fires'] += rows[(i+1) % len(rows)][0] != r[1]
print(dict(st), '%.0fs' % (time.time()-t0))
