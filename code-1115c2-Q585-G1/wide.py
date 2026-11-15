"""Wider test of (MC) and (SW) beyond the 60-state DFS cap: cells with |lam|,|mu| <= 6, up to NR rows, c>=2.
States are sampled by seeded RANDOM mixed walks (a different sampler from the c1/c2 DFS), K states per cell.
Tests at every sampled state: MC (all ordered corner pairs, x->B,y->C), SW column and row, and the B/B control.
Prints a running tally; writes the tally to wide.out every cell.  Engine unchanged."""
import sys, copy, random, time, json
sys.path.insert(0,'/home/clio/projects/proofs/code-1114c2-Q590')
from common import *
from collections import Counter
NR, MAXS, K, TLIM = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
rng = random.Random(20261115); st = Counter(); bad = []; t0 = time.time()
def cells():
    for a in range(1, MAXS+1):
     for c_ in range(1, MAXS+1):
      for lam in parts(a, NR):
       for mu in parts(c_, NR):
        for nu, bs in LR(lam, mu, NR).items():
          if len(bs) < 2: continue
          nt = trim(nu)
          for bb in bs:
            if a <= 5 and c_ <= 5 and NR <= 4: continue      # skip the old census: new cells only
            T = beta_cm(bb, lam, mu, nt)
            b0 = canonical(nt, max(nt[0], 1)+1, len(nt)+1); B = b0.nest[1]
            migrate_flock(b0, 1, {B['pos'][(c, r)]: v for (c, r), v in T.items()}, B['u'], B['w'], 'E')
            yield (lam, mu, nu, bb), b0
def same(p, q): return tiles_of(p) == tiles_of(q)
for key, b0 in cells():
    if time.time() - t0 > TLIM: st['STOPPED at time limit'] = 1; break
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    for _ in range(K):
        b = copy.deepcopy(b0); steps = rng.randrange(0, len(A['pos']))
        for _ in range(steps):
            cs = corners(set(b.nest[0]['pos']))
            if not cs: break
            mig1(b, 0, rng.choice(cs), DV[rng.choice('BC')])
        S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            for y in cs:
                if x == y: continue
                p = copy.deepcopy(b); _, l1, _ = mig1(p, 0, x, DV['B']); _, l2, _ = mig1(p, 0, y, DV['C'])
                q = copy.deepcopy(b); _, m2, _ = mig1(q, 0, y, DV['C']); _, m1, _ = mig1(q, 0, x, DV['B'])
                ok = same(p, q) and (l1, l2) == (m1, m2); st['MC tested'] += 1; st['MC holds'] += ok
                if not ok and len(bad) < 5: bad.append(('MC', key, sorted(S), x, y))
                if x < y:
                    p = copy.deepcopy(b); _, l1, _ = mig1(p, 0, x, DV['B']); _, l2, _ = mig1(p, 0, y, DV['B'])
                    q = copy.deepcopy(b); _, m2, _ = mig1(q, 0, y, DV['B']); _, m1, _ = mig1(q, 0, x, DV['B'])
                    st['BB tested'] += 1; st['BB commutes (control, expect <)'] += same(p, q) and (l1, l2) == (m1, m2)
        for h in cs:
            for nm, g in (('col', (h[0], h[1]-1)), ('row', (h[0]-1, h[1]))):
                if g not in S or g not in corners(S - {h}): continue
                q = copy.deepcopy(b); _, Xh, _ = mig1(q, 0, h, DV['B']); _, Lg, _ = mig1(q, 0, g, DV['C'])
                r = copy.deepcopy(b); _, Lh, _ = mig1(r, 0, h, DV['C']); _, Xg, _ = mig1(r, 0, g, DV['B'])
                ok = same(q, r) and Xh == Xg and Lh == Lg; st['SW-'+nm+' tested'] += 1; st['SW-'+nm+' holds'] += ok
                if not ok and len(bad) < 5: bad.append(('SW-'+nm, key, sorted(S), h, g))
    st['cells'] += 1; st['max |alpha|'] = max(st['max |alpha|'], len(A['pos']))
    json.dump(dict(st=dict(st), bad=[repr(x) for x in bad], secs=time.time()-t0), open('wide.out', 'w'))
print(dict(st)); print(bad)
