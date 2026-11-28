"""PROVE 2026-11-28 c2: (SW) on the UNCAPPED closure (IBC.closure), plus row swap, merge statistic, vertex-diagonal lesion.
usage: python3 swfull.py 5            (size-5 census, parallel over seeds)
       python3 swfull.py 6 CH[,CH]    (census_sized chunks idx%24)"""
import sys, copy
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, '/home/clio/projects/reviews/code-1125c2-BS')
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590')
import migrate
from common import corners, tiles_of, census, parts, LR, trim, beta_cm, canonical, migrate_flock
from flock import unit_dir
from indep_bs_census import closure
import indep_bs_census as IBC

LOG = []
_step = migrate.Board.step
def logged(self, R, d):
    before = set(self.tiles)
    R2 = _step(self, R, d)
    if R2 is not None:
        LOG.append(frozenset(self.tiles[t][1] for t in self.tiles if t not in before))  # created tiles' vertex tuples
    return R2
migrate.Board.step = logged

def mv(b, ij, dvec):
    LOG.clear(); IBC.STEPLOG.clear(); t = b.nest[0]['pos'].pop(ij); R, Y, ij2, s = b.journey(t, 0, unit_dir(dvec)); return (Y, ij2), list(LOG)

def suffix_len(a, b):
    k = 0
    while k < min(len(a), len(b)) and a[-1-k] == b[-1-k]: k += 1
    return k

def test_state(b, dirs, st):
    S = set(b.nest[0]['pos']); cs = corners(S)
    for x in cs:
        i, j = x
        for kind, y in (('col', (i, j-1)), ('row', (i-1, j)), ('diag', (i+1, j-1)), ('diag2', (i-1, j+1))):
            if y not in S: continue
            if kind.startswith('diag'):
                if y not in cs: continue
            elif y not in corners(S - {x}): st[kind+' skip'] += 1; continue
            b1 = copy.deepcopy(b); lxB, Jx1 = mv(b1, x, dirs['B']); lyC, Jy1 = mv(b1, y, dirs['C'])
            b2 = copy.deepcopy(b); lxC, Jx2 = mv(b2, x, dirs['C']); lyB, Jy2 = mv(b2, y, dirs['B'])
            st[kind+' n'] += 1
            same = tiles_of(b1) == tiles_of(b2)
            st[kind+' same tiling'] += same
            st[kind+' C-box equal'] += lyC == lxC
            st[kind+' B-box equal'] += lxB == lyB
            if kind == 'col':
                # merge: y's 2nd journey vs x's solo journey to the same nest (solo = first move on the other side)
                st['col C-side y->C ends with suffix of solo x->C'] += suffix_len(Jy1, Jx2) > 0
                st['col B-side y->B ends with suffix of solo x->B'] += suffix_len(Jy2, Jx1) > 0
                st['col C-side y->C identical to solo x->C'] += Jy1 == Jx2
                st['col B-side y->B identical to solo x->B'] += Jy2 == Jx1
                if not same: st['FAIL'] += 1

def work5(seed_idx):
    st = Counter()
    for k, (key, b0) in enumerate(census()):
        if k != seed_idx: continue
        for b, dirs in closure(b0):
            st['states'] += 1; test_state(b, dirs, st)
    return st

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
def work6(arg):
    w, sub, nsub = arg; st = Counter(); m = 0
    for i, (key, b0) in enumerate(census_sized()):
        if i % NCH != w: continue
        m += 1
        if m % nsub != sub: continue
        try: cl = list(closure(b0))
        except RuntimeError: st['closure exploded'] += 1; continue
        for b, dirs in cl: st['states'] += 1; test_state(b, dirs, st)
    return st

if __name__ == '__main__':
    tot = Counter()
    if sys.argv[1] == '5':
        nseeds = sum(1 for _ in census())
        with Pool(3, maxtasksperchild=1) as p:
            for o in p.imap_unordered(work5, range(nseeds)): tot.update(o)
        print('SIZE 5 full closure, seeds', nseeds)
    else:
        chs = [int(c) for c in sys.argv[2].split(',')]; nsub = 4
        with Pool(3, maxtasksperchild=1) as p:
            for o in p.imap_unordered(work6, [(c, s, nsub) for c in chs for s in range(nsub)]): tot.update(o)
        print('SIZE 6 chunks', chs, 'of', NCH)
    for k in sorted(tot): print(tot[k], k)
