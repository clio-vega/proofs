"""Is the tie-break load-bearing for G1? Re-run the c>=2 census with the tie filter INVERTED (keep hexagons NOT containing
the tile across the horizontal edge), via a patched copy of Board.step. Report R_ii=R_iii and whether readouts changed."""
import sys, types, inspect, textwrap
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
import migrate
src = inspect.getsource(migrate.Board.step)
old = "cands = [zz for zz in cands if across & zz[2]]"; assert src.count(old) == 1
new = "cands = [zz for zz in cands if not (across & zz[2])]" if sys.argv[1] == 'invert' else old
ns = {}; exec(textwrap.dedent(src.replace(old, new)), migrate.__dict__, ns)
patched = ns['step']; orig = migrate.Board.step
from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts
from census import readout2
from collections import Counter
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
        migrate.Board.step = orig
        b = canonical(nt, max(nt[0], 1)+1, len(nt)+1); B = b.nest[1]
        migrate_flock(b, 1, {B['pos'][(c, r)]: v for (c, r), v in T.items()}, B['u'], B['w'], 'E')
        _, R2o, _ = readout2(b, 0, 1, 'N')
        migrate.Board.step = patched
        try:
            _, R2, _ = readout2(b, 0, 1, 'N'); _, R3, _ = readout2(b, 0, -1, 'W')
        except Exception as ex:
            st['ERR '+str(ex)[:50]] += 1; continue
        st['n'] += 1; st['Rii=Riii (patched)'] += norm(R2) == norm(R3); st['patched Rii == original Rii'] += norm(R2) == norm(R2o)
print(sys.argv[1], dict(st))
