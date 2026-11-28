"""PROVE 2026-11-28: ns_census.py on FRESH size-6 chunks (census_sized order of code-1118c2-Q597/bs6.py, idx % 24).
Chunks 0-18 and 20-23 were used by earlier sessions; 19 is fresh.  Uncapped closure (IBC.closure).
usage: python3 ns_census6.py NW CHUNK[,CHUNK...]"""
import sys, copy
sys.path.insert(0, '/home/clio/projects/proofs/code-1128-Q597')
src = open('/home/clio/projects/proofs/code-1128-Q597/ns_census.py').read().split("if __name__")[0]
exec(src)
from common import parts, LR, trim, beta_cm, canonical, migrate_flock
NCH = 24
def census_sized():
    for a in range(1, 7):
     for c_ in range(1, 7):
      if a <= 5 and c_ <= 5: continue
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
def work6(w):
    ST.clear()
    for i, (key, b0) in enumerate(census_sized()):
        if i % NCH != w: continue
        CTX['on'] = False
        try: cl = list(closure(b0))
        except RuntimeError as ex: ST[('closure exploded', str(ex))] += 1; continue
        CTX['on'] = True
        for b, dirs in cl:
            CTX['Tn'] = tn_tris(b); ST['states'] += 1
            for x in corners(set(b.nest[0]['pos'])):
                for T in 'BC':
                    nb = copy.deepcopy(b); t = nb.nest[0]['pos'].pop(x); CTX['k'] = 0; CTX['last_created'] = None
                    nb.journey(t, 0, unit_dir(dirs[T])); ST[('journeys', T)] += 1
        CTX['on'] = False
    return Counter(ST)
if __name__ == '__main__':
    NW = int(sys.argv[1]); chunks = [int(c) for c in sys.argv[2].split(',')]
    tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(work6, chunks): tot.update(o)
    print('SIZE 6 chunks', chunks, 'of', NCH)
    for k in sorted(tot, key=str): print(tot[k], k)
